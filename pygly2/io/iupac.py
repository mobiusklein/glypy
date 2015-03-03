import pygly2
from pygly2.structure import monosaccharide, constants
from pygly2.io import format_constants_map
from pygly2.utils import invert_dict


anomer_map_from = dict(format_constants_map.anomer_map)
anomer_map_to = invert_dict(anomer_map_from)


Stem = constants.Stem
Modification = constants.Modification

def monosaccharide_to_iupac(residue):
    template = "{anomer}-{configuration}-{modification}{base_type}{ring_type}{substituent}"
    anomer = anomer_map_to[residue.anomer]
    configuration = residue.configuration[0].name.upper()
    modification = ""
    if Modification.d in residue.modifications.values():
        modification = 'd'
    base_type = ""
    if len(residue.stem) == 1:
        base_type = residue.stem[0].name.title()
    else:
        base_type = resolve_special_base_type(residue)
    ring_type = residue.ring_type.name[0]
    substituent = resolve_substituent(residue)
    return template.format(
        anomer=anomer,
        configuration=configuration,
        modification=modification,
        base_type=base_type,
        ring_type=ring_type,
        substituent=substituent
        )


def resolve_substituent(residue):
    substituent = ""
    names = {sub.name for p, sub in residue.substituents()}
    if 'n_acetyl' in names:
        substituent = "NAc"
    elif 'n_glycolyl' in names:
        substituent = 'NGc'
    try:
        if resolve_special_base_type(residue) == 'Neu' and substituent[0] == 'N':
            substituent = substituent[1:]
    except ValueError:
        pass
    return substituent


def resolve_special_base_type(residue):
    if residue.superclass == "NON":
        if residue.stem == (Stem.gro, Stem.gal):
            return "Neu"
    else:
        raise ValueError("Could not resolve special base_type for {}".format(residue))


def to_iupac(structure):
    if isinstance(structure, Monosaccharide):
        return monosaccharides_to_iupac(structure)
    else:
        raise NotImplementedError("Non-monosaccharide translations not supported yet.")