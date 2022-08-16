'''
CFG Format
----------

An experimental parser for the `Consortium for Functional Glycomics <http://www.functionalglycomics.org/>`_
(CFG) glycan line format.

'''
import re
import warnings

from collections import deque, namedtuple
from functools import partial

from glypy.structure import (
    Monosaccharide, Glycan, Link, AmbiguousLink,
    Substituent, constants, named_structures,
    UnknownPosition, SuperClass, Modification)
from glypy.composition import Composition
from glypy.composition.structure_composition import substituent_compositions
from glypy.composition.composition_transform import has_derivatization, derivatize
from glypy.io import format_constants_map
from glypy.io.nomenclature import identity
from glypy.utils import invert_dict

from glypy.io.file_utils import ParserInterface, ParserError


# A static copy of monosaccharide names to structures for copy-free comparison
monosaccharide_reference = {k: v for k, v in named_structures.monosaccharides.items()}

anomer_map_from = dict(format_constants_map.anomer_map)
anomer_map_from['?'] = anomer_map_from.pop('x')
anomer_map_to = invert_dict(anomer_map_from)
anomer_map_from['beta'] = anomer_map_from['b']
anomer_map_from['alpha'] = anomer_map_from['a']
anomer_map_from[u"\u03B1"] = anomer_map_from['a']
anomer_map_from[u"\u03B2"] = anomer_map_from['b']

def _make_substituent_name(name):
    return ''.join(t.title() for t in name.split("_")).replace("(", "").replace(")", "")


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

substituents_map_from = invert_dict(substituents_map_to)
substituents_map_from['Phosphate'] = 'phosphate'


class CFGError(ParserError):
    pass


LinkageSpecification = namedtuple("LinkageSpecification", ("child_position", "parent_position", "has_ambiguity"))

class LinkageDeserializer(object):
    pattern = re.compile(r"(?P<child_linkage>[0-9?/]+)->?(?P<parent_linkage>[0-9?/]+)?")

    def parse(self, linkage_string):
        if linkage_string is None:
            return None
        match = self.pattern.search(linkage_string)
        if match is None:
            raise CFGError(linkage_string)
        has_ambiguity = "/" in linkage_string
        match_groups = match.groupdict()
        child_linkage = match_groups['child_linkage']
        if child_linkage is not None:
            child_linkage = self.parse_position(child_linkage)
        parent_linkage = match_groups['parent_linkage']
        if parent_linkage is not None:
            parent_linkage = self.parse_position(parent_linkage)
        return LinkageSpecification(child_linkage, parent_linkage, has_ambiguity)

    def parse_position(self, position_string):
        if position_string == '?':
            return UnknownPosition
        elif '/' in position_string:
            return list(map(int, position_string.split('/')))
        else:
            return int(position_string)

    def __call__(self, linkage_string):
        return self.parse(linkage_string)


class SubstituentDeserializer(object):
    def __init__(self, error_on_missing=True):
        self.error_on_missing = error_on_missing

    def substituent_from_cfg(self, tokens):
        for position, name in tokens:
            if position is None:
                position = UnknownPosition
            else:
                try:
                    position = int(position)
                except (ValueError, TypeError):
                    warnings.warn("Unable to interpret substituent position %r" % (position))
                    position = UnknownPosition

            try:
                name = (substituents_map_from[name])
            except KeyError:
                # Acidic special case:
                # Often, acidic monosaccharides are written with a trailing A like a substituent while
                # GlycoCT treats acidic groups as modifications. If an A appears in the substituent suffix
                # it will fail to be cast as a substituent, but pass it along raw and it will be handled
                # downstream by :func:`monosaccharide_from_iupac`.
                if name == "A":
                    pass
                else:  # pragma: no cover
                    if not self.error_on_missing:
                        warnings.warn("No translation rule found to convert %s into a Substituent" % name)
                        continue
                    else:
                        raise CFGError("No translation rule found to convert %s into a Substituent" % name)
            yield int(position), name

    def symbol_to_name(self, symbol):
        name = substituents_map_from[symbol]
        return name

    def __call__(self, substituents):
        return self.substituent_from_cfg(substituents)


def aminate_substituent(substituent):
    if substituent.name.startswith("n_"):
        # already aminated
        return substituent
    aminated = Substituent("n_" + substituent.name)
    if aminated.composition == {}:
        raise ValueError("Could not aminate substituent")
    return aminated


class MonosaccharideDeserializer(object):
    _pattern = r"""
        (:?\((?P<substituent_prefix_position>[0-9\?]+?)
             (?P<substituent_prefix_name>[A-Za-z]+?)\))?
        (?P<base_type>(?:[A-Z][a-z]{2}?|(:?[a-z]{3}[A-Z][a-z]{2})))
        (?:(?P<substituent_position>\d+?)?(?P<substituent_name>[A-Za-z]+?))?
        (?P<anomer>a|b|\?|alpha|beta|\u03B1|\u03B2)
        (?P<linkage>[0-9?/]+(:?->?|,)[0-9?/]+?)?
        $"""
    try:
        # convert to unicode for Py2
        _pattern = _pattern.decode("raw_unicode_escape")
    except AttributeError:
        pass
    pattern = re.compile(_pattern, re.VERBOSE | re.UNICODE)

    def __init__(self, substituent_deserializer=None):
        if substituent_deserializer is None:
            substituent_deserializer = SubstituentDeserializer()
        self.substituent_deserializer = substituent_deserializer
        self.linkage_parser = LinkageDeserializer()

    def has_pattern(self, string):
        return self.pattern.search(string)

    def extract_pattern(self, monosaccharide_str):
        match = self.pattern.search(monosaccharide_str)
        if match is None:
            raise CFGError("Cannot find monosaccharide pattern in {}".format(monosaccharide_str))
        match_dict = match.groupdict()
        return match_dict

    def build_residue(self, match_dict):
        try:
            anomer = anomer_map_from[match_dict['anomer']]
        except KeyError:
            anomer = anomer_map_from['?']
        base_type = match_dict["base_type"]
        linkage = match_dict.get("linkage")
        original_base_type = base_type
        try:
            residue = named_structures.monosaccharides[base_type]
        except KeyError:
            raise CFGError("Unknown Residue Base-type %r" % (original_base_type,))

        base_is_modified = len(residue.substituent_links) + len(residue.modifications) > 0
        residue.anomer = anomer
        substituents = []
        if match_dict['substituent_name'] is not None:
            substituents.append((match_dict['substituent_position'],
                                 match_dict['substituent_name']))
        if match_dict['substituent_prefix_name'] is not None:
            substituents.append((
                match_dict['substituent_prefix_position'],
                match_dict['substituent_prefix_name']))
        self.set_substituents(residue, substituents, base_is_modified, base_type)
        return residue, linkage

    def set_substituents(self, residue, substituent_set, base_is_modified, base_type):
        i = 0
        for position, substituent in self.substituent_deserializer(substituent_set):
            i += 1
            if position == UnknownPosition and base_is_modified:
                # Guess at what the user might mean using base_type
                if base_type == "Neu" and substituent in ["acetyl", "glycolyl"] and i == 1:
                    position = 5
                else:
                    raise ValueError(
                        "Cannot have ambiguous location of substituents on a base type which"
                        " has default modifications or substituents. {} {}".format(
                            residue, (position, substituent)))
            # Often, acidic monosaccharides will be suffixed "A" instead of prefixed "a".
            # Handle this here.
            if substituent == "A":
                residue.add_modification(Modification.a, position)
                continue

            substituent = Substituent(substituent)
            try:
                residue.add_substituent(
                    substituent, position,
                    parent_loss=substituent.attachment_composition_loss(), child_loss='H')
            except ValueError:
                # Highly modified large bases have a degenerate encoding, where additional qualifications following
                # base name *replace* an existing substituent. This behavior may not be expected in other more
                # common cases.
                if base_type in {"Neu", "Kdo"}:
                    occupancy = 0
                    try:
                        unplaced = residue.substituent_links[position][0].child
                        residue.drop_substituent(position)
                        if unplaced.name == "amino":
                            try:
                                substituent = aminate_substituent(substituent)
                            except ValueError:
                                pass
                    except ValueError:
                        # The site contains a modification which can be present alongside the substituent
                        occupancy = 1
                    except IndexError:
                        occupancy = 1
                    try:
                        residue.add_substituent(
                            substituent, position, occupancy,
                            parent_loss=substituent.attachment_composition_loss(), child_loss='H')
                    except ValueError:
                        raise CFGError("Can't resolve %s" % substituent)
                else:
                    raise

    def add_monosaccharide_bond(self, residue, parent, linkage):
        if parent is not None and linkage != ():
            parent.add_monosaccharide(residue, position=linkage[1], child_position=linkage[0])

    def monosaccharide_from_cfg(self, monosaccharide_str, parent=None):
        match_dict = self.extract_pattern(monosaccharide_str)
        residue, linkage = self.build_residue(match_dict)
        linkage = self.linkage_parser(linkage)

        self.add_monosaccharide_bond(residue, parent, linkage)
        return residue, linkage

    def __call__(self, cfg_str, parent=None):
        return self.monosaccharide_from_cfg(cfg_str, parent=parent)

    def finalize(self, glycan):
        pass


class SpacerDeserializer(object):
    pattern = re.compile(r"""
    -(?P<spacer_name>[A-Za-z]+)
    (?P<linkage>\d+)$
    """, re.VERBOSE | re.UNICODE)

    def has_pattern(self, string):
        return self.pattern.search(string)

    def remove_pattern(self, cfg_str):
        return self.pattern.sub("", cfg_str)

    def extract_pattern(self, cfg_str):
        match = self.pattern.search(cfg_str)
        if match is None:
            raise CFGError("Cannot find spacer pattern in {}".format(cfg_str))
        match_dict = match.groupdict()
        return match_dict

    def spacer_from_cfg(self, cfg_str):
        return


class GlycanDeserializer(object):
    def __init__(self, monosaccharide_deserializer=None, set_default_positions=True):
        if monosaccharide_deserializer is None:
            monosaccharide_deserializer = MonosaccharideDeserializer()
        self.monosaccharide_deserializer = monosaccharide_deserializer
        self.set_default_positions = set_default_positions
        self.spacer_deserializer = SpacerDeserializer()

    new_branch_open = re.compile(r"(\))$")

    def add_monosaccharide(self, parent_node, child_node, linkage):
        self.monosaccharide_deserializer.add_monosaccharide_bond(
            child_node, parent_node, linkage)

    def glycan_from_cfg(self, text, structure_class=Glycan, **kwargs):
        last_outedge = None
        root = None
        last_residue = None
        branch_stack = []

        # Remove the base
        if self.spacer_deserializer.has_pattern(text):

            text = self.spacer_deserializer.remove_pattern(text)

        while len(text) > 0:
            # If starting a new branch
            match = self.new_branch_open.search(text)
            if match is not None:
                step = match.end(1) - match.start(1)
                text = text[:-step]
                branch_stack.append((last_residue, root, last_outedge))
                root = None
                last_residue = None
                last_outedge = None
            # If ending a branch
            elif text[-1] == '(':
                try:
                    branch_parent, old_root, old_last_outedge = branch_stack.pop()
                    # child_position, parent_position = last_outedge
                    self.add_monosaccharide(branch_parent, root, last_outedge)
                    root = old_root
                    last_residue = branch_parent
                    last_outedge = old_last_outedge
                    text = text[:-1]
                except IndexError:
                    raise CFGError("Bad branching at {}".format(len(text)))
            # Parsing a residue
            else:
                match = self.monosaccharide_deserializer.has_pattern(text)
                if match:
                    next_residue, outedge = self.monosaccharide_deserializer(
                        text[match.start(): match.end()], last_residue)
                    if root is None:
                        last_outedge = outedge
                        root = next_residue
                    last_residue = next_residue
                    text = text[:match.start()]
                else:
                    raise CFGError("Could not identify residue '...{}' at {}".format(text[-30:], len(text)))

        res = structure_class(root=root)
        self.monosaccharide_deserializer.finalize(res)
        res.reindex()
        if self.set_default_positions:
            self.set_default_positions_for_common_cases(res)
        return res

    def __call__(self, text, **kwargs):
        return self.glycan_from_cfg(text, **kwargs)

    def set_default_positions_for_common_cases(self, glycan):
        for node in glycan:
            candidates = []
            for pos, link in list(node.substituent_links.items()):
                if pos == UnknownPosition:
                    candidates.append(link)
            for candidate in candidates:
                substituent = candidate.to(node)
                position = None
                if substituent.name == 'n_acetyl':
                    if node.superclass == SuperClass.hex:
                        position = 2
                    elif node.superclass == SuperClass.non:
                        position = 5
                if position is None:
                    continue
                if not node.is_occupied(position):
                    candidate.break_link(refund=True)
                    candidate.parent_position = position
                    candidate.apply()
        return glycan


glycan_parser = GlycanDeserializer()


def loads(text):
    '''Parse a single CFG glycan sequence.

    .. note::
        The spacer, if any, is ignored.

    Parameters
    ----------
    text : str
        The sequence to parse

    Returns
    -------
    structure : :class:`~.Glycan`
        The parsed glycan structure
    '''
    return glycan_parser(text)
