'''
Glyconnect
----------

A simple dialect of the Glyconnect/GlycoMod glycan composition notation.
'''
import re

from glypy.structure import glycan_composition
from glypy.structure.glycan_composition import (
    FrozenGlycanComposition,
    FrozenMonosaccharideResidue,
    SubstituentResidue)
from glypy.structure.glycan import Glycan
from glypy.utils import invert_dict


try:
    import requests
except ImportError:
    requests = None

#: The set of defined symbols and their mappings.
defined_symbols = {
    "Hex": FrozenMonosaccharideResidue.from_iupac_lite("Hex"),
    "HexNAc": FrozenMonosaccharideResidue.from_iupac_lite('HexNAc'),
    "dHex": FrozenMonosaccharideResidue.from_iupac_lite('dHex'),
    "NeuAc": FrozenMonosaccharideResidue.from_iupac_lite("NeuAc"),
    "NeuGc": FrozenMonosaccharideResidue.from_iupac_lite("NeuGc"),
    "S": SubstituentResidue("sulfate"),
    "P": SubstituentResidue("phosphate"),
    "Xyl": FrozenMonosaccharideResidue.from_iupac_lite("Xyl"),
    "HexA": FrozenMonosaccharideResidue.from_iupac_lite("HexA"),
    "Pent": FrozenMonosaccharideResidue.from_iupac_lite("Pen"),
    "Kdn": FrozenMonosaccharideResidue.from_iupac_lite("Kdn"),
    "Su": SubstituentResidue("sulfate"),
}


monosaccharide_to_symbol = invert_dict(defined_symbols)


tokenizer = re.compile(r"([^:\s]+):(\d+)")


def loads(string):
    '''Parse a GlyConnect glycan composition into a :class:`~.FrozenGlycanComposition`

    Parameters
    ----------
    string: str
        The string to parse

    Returns
    -------
    :class:`~.FrozenGlycanComposition`

    Raises
    ------
    :class:`KeyError`: Raised if a key isn't defined by the GlyConnect dialect
    '''
    tokens = tokenizer.findall(string)
    gc = FrozenGlycanComposition()
    for mono, count in tokens:
        mono = defined_symbols[mono]
        count = int(count)
        gc[mono] += count
    return gc


def dumps(composition):
    '''Encode :class:`~.GlycanComposition` or :class:`~.Glycan` into the GlyConnect
    glycan composition text format.

    Parameters
    ----------
    composition: :class:`~.GlycanComposition` or :class:`~.Glycan`
        The structure to format

    Returns
    -------
    :class:`str`

    Raises
    ------
    :class:`KeyError`: Raised if a key isn't defined by the GlyConnect Compozitor dialect
    '''
    if isinstance(composition, Glycan):
        composition = FrozenGlycanComposition.from_glycan(composition)
    tokens = []
    for key, value in composition.items():
        key = monosaccharide_to_symbol[key]
        tokens.append("%s:%d" % (key, value))
    return ' '.join(tokens)


API_SERVER = "https://glyconnect.expasy.org/api"


def from_glytoucan_id(glytoucan_id):
    response = requests.post(
        "{server}/structures/search/glytoucan".format(server=API_SERVER),
        {"glytoucan_id": glytoucan_id})
    response.raise_for_status()
    data = response.json()
    return data

