import re

from collections import deque, namedtuple

from glypy.structure import (
    Monosaccharide, Glycan, Link, AmbiguousLink,
    Substituent, constants, named_structures, UnknownPosition)

from glypy.composition import Composition
from glypy.composition.structure_composition import substituent_compositions
from glypy.composition.composition_transform import has_derivatization, derivatize
from glypy.io import format_constants_map
from glypy.io.nomenclature import identity
from glypy.utils import invert_dict


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

reducing_or_free = re.compile(r"(redEnd|freeEnd)")

residue_pattern = re.compile(r"""^(?P<monosaccharide>(?P<configuration>[DL?])-
(?P<base_type>(?:[A-Z][a-z]{2}?|(?:[a-z]{3}[A-Z][a-z]{2})))
(?P<bound_substituent>Ac|NAc|N|S|P)?,
(?P<ring_type>[fpo?]))|
(?P<substituent>Ac|NAc|N|S|P)|
(?P<modification>deoxy)
""", re.VERBOSE)

linkage_pattern = re.compile(
    r"^--(?P<parent_position>\?|\d+)(?:(?P<anomer>[ab?])(?P<child_position>\?|\d+)?)?")


def parse_gws(text):
    reducing_end = reducing_or_free.search(text)
    if reducing_end:
        text = text[reducing_end.end():]
    while text:
        if text[0] == "(":
            raise ValueError("No support for branches yet")
        elif text[0] == '$':
            print("Hit configuration, quitting")
            break
        print(len(text))
        linkage = linkage_pattern.search(text)
        if linkage:
            text = text[linkage.end():]
        else:
            raise ValueError(
                "Failed to parse linkage from \"%s\"..." % (text[:10], ))
        residue = residue_pattern.search(text)
        if residue:
            text = text[residue.end():]
        else:
            raise ValueError(
                "Failed to parse residue from \"%s\"..." % (text[:20], ))
        print(linkage.groupdict())
        print(residue.groupdict())


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
