'''
GWS Format
----------

A parser for a subset of the GlycoWorkbench sequence format.
'''
import re
import warnings
from collections import deque

from glypy.structure import (
    Monosaccharide, Glycan, Link, AmbiguousLink,
    Substituent, constants, named_structures, UnknownPosition)

from glypy.composition.structure_composition import substituent_compositions

from glypy.utils import invert_dict

from glypy.io import format_constants_map
from glypy.io.nomenclature import identity
from glypy.io import iupac, file_utils

monosaccharide_reference = {k: v for k,
                            v in named_structures.monosaccharides.items()}

special_base_types = {
    # "Neu5Ac", "Neu5Gc", "Neu",
    # "Kdn", "Kdo",
    "Oli", "Tyv",
    "Psi", "Fru", "Sor", "Tag",
    "Xul", "Sed"
}

special_base_types = {
    s: monosaccharide_reference[s]
    for s in special_base_types
}

special_base_type_resolver = identity.MonosaccharideIdentifier(
    special_base_types)

anomer_map_from = dict(format_constants_map.anomer_map)
anomer_map_from['?'] = anomer_map_from.pop('x')
anomer_map_to = invert_dict(anomer_map_from)

Stem = constants.Stem
Configuration = constants.Configuration
Modification = constants.Modification
SuperClass = constants.SuperClass


deserializer = iupac.MonosaccharideDeserializer()

reducing_or_free = re.compile(r"(redEnd|freeEnd)")

residue_pattern = re.compile(r"""^(?:(?P<monosaccharide>(?:(?P<configuration>[DL?])-)?
(?P<base_type>(?:[A-Z][a-z]{2}?|(?:[a-z]{3}[A-Z][a-z]{2})))
(?P<bound_substituent>Ac|NAc|N|S|P|A)?,
(?P<ring_type>[fpo?]))|
(?P<substituent>Ac|NAc|N|S|P|A)|
(?P<modification>deoxy))
""", re.VERBOSE)

linkage_pattern = re.compile(
    r"^--(?P<parent_position>\?|\d+)(?:(?P<anomer>[ab?])(?P<child_position>\?|\d+)?)?")


def _make_substituent_name(name):
    return ''.join(t.title() for t in name.split("_")).replace("(", "").replace(")", "")


class GWSError(file_utils.ParserError):
    pass


substituents_map_to = {
    name: _make_substituent_name(name) for name in substituent_compositions
}

# Special Cases
substituents_map_to['n_acetyl'] = "NAc"
substituents_map_to['n_glycolyl'] = "NGc"
substituents_map_to['n_sulfate'] = "NS"
substituents_map_to['sulfate'] = "S"
substituents_map_to["methyl"] = "Me"
substituents_map_to["acetyl"] = "Ac"
substituents_map_to["glycolyl"] = "Gc"
substituents_map_to["fluoro"] = "F"
substituents_map_to["amino"] = "N"
substituents_map_to['phosphate'] = 'P'
substituents_map_to['phospho_ethanolamine'] = 'PEtn'
substituents_map_to['ethanolamine'] = 'Etn'

substituents_map_from = invert_dict(substituents_map_to)
substituents_map_from['Phosphate'] = 'phosphate'

_modification_map_to = {
    'deoxy': 'd',
}


_substituent_replacement_rules = {
    'NeuAc': [
        ('n_acetyl', 'acetyl')
    ],
    'NeuGc': [
        ('n_glycolyl', 'glycolyl')
    ],
    'Neu': [
        ('amino', None)
    ]
}


def set_ring_bounds(residue, ring_type):
    if residue.ring_start == UnknownPosition:
        residue.ring_end = UnknownPosition
    elif ring_type == 'p':
        residue.ring_end = residue.ring_start + 4
    elif ring_type == 'f':
        residue.ring_end = residue.ring_start + 3
    elif ring_type == 'o':
        residue.ring_end = residue.ring_start = 0
    else:
        residue.ring_end = residue.ring_start = UnknownPosition


def build_residue(residue_dict):
    if residue_dict['monosaccharide']:
        residue_dict['substituent'] = residue_dict.pop(
            'bound_substituent', '') or ''
        residue_dict['configuration'] = residue_dict['configuration'] or '?'
        residue, _ = deserializer.build_residue(residue_dict)
        return residue, 'monosaccharide'
    elif residue_dict['substituent']:
        substs = [subst for pos, subst in deserializer.substituent_parser(
            residue_dict['substituent'])]
        return substs[0], 'substituent'
    elif residue_dict['modification']:
        return Modification[residue_dict['modification']], 'modification'
    else:
        raise ValueError("Could not determine node type from %r" %
                         (residue_dict, ))


def apply_linkage(parent, child, linkage_spec, child_type):
    if linkage_spec['anomer'] and child_type == 'monosaccharide':
        child.anomer = linkage_spec['anomer']
    if parent is None:
        return
    child_position = linkage_spec['child_position']
    if child_position[:1].isdigit():
        child_position = int(child_position)
    else:
        child_position = UnknownPosition
    parent_position = linkage_spec['parent_position']
    if parent_position[:1].isdigit():
        parent_position = int(parent_position)
    else:
        parent_position = UnknownPosition
    if child_type == 'monosaccharide':
        parent.add_monosaccharide(child, parent_position, child_position)
    elif child_type == 'substituent':
        parent.add_substituent(child, parent_position)
    elif child_type == "modification":
        parent.add_modification(child, parent_position)
    else:
        raise GWSError(child_type)


def loads(text):
    '''Parse a single GWS glycan sequence plus metadata.

    Parameters
    ----------
    text : str
        The sequence to parse

    Returns
    -------
    structure : :class:`~.Glycan`
        The parsed glycan structure
    metadata : :class:`tuple`
        Unstructured metadata associated with the sequence.
    '''
    reducing_end = reducing_or_free.search(text)
    if reducing_end:
        text = text[reducing_end.end():]
    branch_queue = deque()
    root = None
    last_residue = root
    i = 0
    while text:
        if text[0] == "(":
            while text[0] == "(":
                text = text[1:]
                branch_queue.append((last_residue, i))
        elif text[0] == ')':
            text = text[1:]
            last_residue, last_id = branch_queue.pop()
        elif text[0] == '$':
            metadata = text[1:].split(",")
            break
        elif text[0] == '}':
            warnings.warn("Unknown symbol }, ignoring")
            text = text[1:]
        i += 1
        linkage = linkage_pattern.search(text)
        if linkage:
            text = text[linkage.end():]
        else:
            raise GWSError(
                "Failed to parse linkage from \"%s\"..." % (text[:20], ))
        residue = residue_pattern.search(text)
        if residue:
            text = text[residue.end():]
        else:
            raise GWSError(
                "Failed to parse residue from \"%s\"..." % (text[:20], ))
        residue, node_type = build_residue(residue.groupdict())
        apply_linkage(last_residue, residue, linkage.groupdict(), node_type)
        if root is None:
            root = residue
        if node_type == 'monosaccharide':
            # refuse to handle bridge substituents
            last_residue = residue
    return Glycan(root).reindex(), metadata


parse_gws = loads


test_data = '''redEnd--?a1D-GalNAc,p--3b1D-Gal,p--3a2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p--3b1D-Gal,p--3a2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p--3b1D-Gal,p--3a2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--??1Hex,p--??1Hex,p--??1Hex,p$MONO,Und,-H,0,redEnd
redEnd--??1Hex,p--??1Hex,p--??1Hex,p$MONO,Und,-H,0,redEnd

redEnd--??1Hex,p--??1Hex,p--??1Hex,p--??1Hex,p$MONO,Und,-H,0,redEnd
redEnd--??1Hex,p--??1Hex,p--??1Hex,p$MONO,Und,-H,0,redEnd
redEnd--??1Hex,p--??1Hex,p--??1Hex,p--??1Hex,p$MONO,Und,-H,0,redEnd
redEnd--??1Hex,p--??1Hex,p(--??1Hex,p)--??1Hex,p$MONO,Und,-H,0,redEnd
redEnd--??1Hex,p--??1Hex,p(--??1Hex,p)--??1Hex,p$MONO,Und,-H,0,redEnd

redEnd--??1Hex,p--??1Hex,p--??1Hex,p(--??1Hex,p)--??1Hex,p$MONO,Und,-H,0,redEnd

redEnd--??1Hex,p--??1Hex,p--??1Hex,p(--??1Hex,p)--??1Hex,p$MONO,Und,-H,0,redEnd
redEnd--??1Hex,p--??1Hex,p--??1Hex,p(--??1Hex,p)--??1D-GlcA,p$MONO,Und,-H,0,redEnd
redEnd--??1Hex,p--??1Hex,p--??1Hex,p(--??1Hex,p)--??1D-GlcA,p$MONO,Und,-H,0,redEnd
redEnd--??1Hex,p--??1Hex,p(--??1Hex,p--??1Hex,p--??1Hex,p)--??1D-GlcA,p$MONO,Und,-H,0,redEnd

redEnd--??1Hex,p--??1Hex,p(--??1Hex,p--??1Hex,p--??1Hex,p)--??1D-GlcA,p$MONO,Und,-H,0,redEnd
redEnd--??1Hex,p--??1Hex,p--??1Hex,p--??1Hex,p--??1Hex,p}--??1D-GlcA,p$MONO,Und,-H,0,redEnd
redEnd--??1Hex,p--??1Hex,p--??1Hex,p--??1Hex,p--??1Hex,p}--??1D-GlcA,p$MONO,Und,-H,0,redEnd

redEnd--??1Hex,p--??1Hex,p--??1Hex,p--??1Hex,p--??1Hex,p(--??1Hex,p)--??1D-GlcA,p$MONO,Und,-H,0,redEnd
redEnd--??1Hex,p--??1Hex,p--??1Hex,p--??1Hex,p--??1Hex,p(--??1Hex,p)--??1D-GlcA,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p--3b1D-Gal,p--4a1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p--3b1D-Gal,p--4a1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p)--3b1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p)--3b1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-GlcNAc,p--4b1D-Gal,p)--6a2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-GlcNAc,p--4b1D-Gal,p)--6a2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--3a2D-NeuAc,p)--6a2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--3a2D-NeuAc,p)--6a2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p--6b1D-GlcNAc,p--4b1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p--6b1D-GlcNAc,p--4b1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--4a1D-GlcNAc,p)--6b1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--4a1D-GlcNAc,p)--6b1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p)--3b1D-Gal,p--??1D-GlcNAc,p((--??1D-Gal,p)--??1L-Fuc,p)--??1S$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p)--3b1D-Gal,p--??1D-GlcNAc,p((--??1D-Gal,p)--??1L-Fuc,p)--??1S$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p((--??1D-Gal,p)--??1L-Fuc,p)--??1S$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p((--??1D-Gal,p)--??1L-Fuc,p)--??1S$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--6b1D-Gal,p)--3b1D-GlcNAc,p(--3a1L-Fuc,p)--4b1D-Gal,p)--6b1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--6b1D-Gal,p)--3b1D-GlcNAc,p(--3a1L-Fuc,p)--4b1D-Gal,p)--6b1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--3a2D-NeuAc,p)--6b1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--3a2D-NeuAc,p)--6b1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p)--3b1D-Gal,p--??1D-GlcNAc,p(--??1L-Fuc,p)--??1S$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p(--??1S)--??1L-Fuc,p)--3b1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p)--3b1D-Gal,p--??1D-GlcNAc,p(--??1L-Fuc,p)--??1S$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p(--??1S)--??1L-Fuc,p)--3b1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p--3b1D-Gal,p--??1D-GlcNAc,p((--??1S)--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p--3b1D-Gal,p--??1D-GlcNAc,p((--??1S)--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p(--??1S)--??1L-Fuc,p)--3b1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p(--??1S)--??1L-Fuc,p)--3b1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd

redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--4b1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--4b1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1L-Fuc,p)--??1HexNAc,p$MONO,Und,-H,0,redEnd

redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1L-Fuc,p)--??1HexNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd

redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??2D-NeuAc,p)--??1D-GlcNAc,p--??1D-GalNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??2D-NeuAc,p)--??1D-GlcNAc,p--??1D-GalNAc,p$MONO,Und,-H,0,redEnd

redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--3b1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--3b1D-GlcNAc,p$MONO,Und,-H,0,redEnd

redEnd--?a1D-GalNAc,p(--3b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??2D-NeuAc,p)--6b1D-GlcNAc,p--??1D-GalNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??2D-NeuAc,p)--6b1D-GlcNAc,p--??1D-GalNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??2D-NeuAc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??2D-NeuAc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p)--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p)--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd

redEnd--?a1D-GalNAc,p(--3b1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--??1D-GalNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--??1D-GalNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1L-Fuc,p)--??1D-GalNAc,p)--6b1D-GlcNAc,p$MONO,Und,-H,0,redEnd

redEnd--?a1D-GalNAc,p(--3b1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-H,0,redEnd

redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p}--??2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??2D-NeuAc,p)--6b1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p--??2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--??1D-GalNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p)--??1D-GlcNAc,p--??1D-Gal,p--??2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p)--??1D-GlcNAc,p--??1D-Gal,p--??2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??2D-NeuAc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??2D-NeuAc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p(--??1S)--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p--??1S--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p(--??1S)--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p--??1S--??1D-Gal,p$MONO,Und,-H,0,redEnd

redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GalNAc,p)--3b1D-GlcNAc,p}--??1L-Fuc,p$MONO,Und,-H,0,redEnd

redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GalNAc,p)--3b1D-GlcNAc,p}--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p--??2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??2D-NeuAc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p--??2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??2D-NeuAc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p(--??1S)--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p(--??1S)--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-H,0,redEnd

redEnd--?a1D-GalNAc,p(--3b1D-GlcNAc,p--??1D-Gal,p--??2D-NeuAc,p)--6b1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-GlcNAc,p--??1D-Gal,p--??2D-NeuAc,p)--6b1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1L-Fuc,p)--??2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1L-Fuc,p)--??2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p--3b1D-Gal,p(--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-H,0,redEnd

redEnd--?a1D-GalNAc,p--3b1D-Gal,p(--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1L-Fuc,p)--??1D-GalNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1L-Fuc,p)--??1D-GalNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1L-Fuc,p)--??1D-GalNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1L-Fuc,p)--??1D-GalNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p)--??1D-Gal,p(--??1D-GlcNAc,p)--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p)--??1D-Gal,p(--??1D-GlcNAc,p)--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p(--??1S)--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p(--??1S)--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p--??2D-NeuAc,p$MONO,Und,-H,0,redEnd

redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p--??2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p--3b1D-Gal,p--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p--3b1D-Gal,p--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-H,0,redEnd

redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd

redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p)--??1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p)--??1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p)--??1D-Gal,p(--??1D-GlcNAc,p)--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p)--??1D-Gal,p(--??1D-GlcNAc,p)--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-H,0,redEnd

redEnd--?a1D-GalNAc,p(--3b1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p$MONO,Und,-H,0,redEnd

redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p)--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p)--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p(--??1S)--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p(--??1S)--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd

redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p)--??1D-GlcNAc,p--??1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p)--??1D-GlcNAc,p--??1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-Gal,p)--??1L-Fuc,p)--6b1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-Gal,p)--??1L-Fuc,p)--6b1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p)--??1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p)--??1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p)--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p)--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p)--??1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p)--??1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p)--??1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p)--??1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??2D-NeuAc,p)--6b1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p--??2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??2D-NeuAc,p)--6b1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p--??2D-NeuAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p(--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p(--??1L-Fuc,p)--??1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p(--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p(--??1L-Fuc,p)--??1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p)--6b1D-GlcNAc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p)--6b1D-GlcNAc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p--??1D-Gal,p(--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p--??1D-Gal,p(--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p)--6b1D-GlcNAc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p)--6b1D-GlcNAc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd

redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-H,0,redEnd

redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1L-Fuc,p)--??1D-GalNAc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1L-Fuc,p)--??1D-GalNAc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p(--??1L-Fuc,p)--??1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p(--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p(--??1L-Fuc,p)--??1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p(--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p(--??1D-Gal,p--??1L-Fuc,p)--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-2H,0,redEnd

redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p(--??1D-Gal,p--??1L-Fuc,p)--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-2H,0,redEnd



redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-2H,0,redEnd


redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-2H,0,redEnd


redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--6b1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p$MONO,Und,-H,0,redEnd

redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--6b1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p$MONO,Und,-H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p)--6b1D-GlcNAc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p(--??1L-Fuc,p)--??1D-Gal,p)--6b1D-GlcNAc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p(--??1L-Fuc,p)--??1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p(--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p(--??1L-Fuc,p)--??1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p(--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p(--??1D-GalNAc,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd

redEnd--?a1D-GalNAc,p(--3b1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p--??1D-GlcNAc,p--??1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p(--??1D-Gal,p--??1L-Fuc,p)--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p(--??1D-Gal,p--??1L-Fuc,p)--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p(--??1D-Gal,p--??1L-Fuc,p)--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p((--??1D-Gal,p)--??1L-Fuc,p)--??1L-Fuc,p)--??1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--6b1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p--??1D-Gal,p--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p--??1D-GlcNAc,p)--6b1D-GlcNAc,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-2H,0,redEnd
redEnd--?a1D-GalNAc,p(--3b1D-Gal,p(--??1D-GlcNAc,p(--??1D-Gal,p)--??1L-Fuc,p)--??1D-GlcNAc,p--??1D-Gal,p)--6b1D-GlcNAc,p--??1D-Gal,p$MONO,Und,-2H,0,redEnd
'''
