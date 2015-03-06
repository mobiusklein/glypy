from collections import deque

import pygly2
from pygly2.structure import Monosaccharide, Glycan, constants, named_structures
from pygly2.io import format_constants_map
from pygly2.io.nomenclature import identity
from pygly2.utils import invert_dict


# A static copy of monosaccharide names to structures for copy-free comparison
monosaccharide_reference = {k: v for k, v in named_structures.monosaccharides.items()}


anomer_map_from = dict(format_constants_map.anomer_map)
anomer_map_from['?'] = anomer_map_from.pop('x')
anomer_map_to = invert_dict(anomer_map_from)


Stem = constants.Stem
Modification = constants.Modification
SuperClass = constants.SuperClass


def monosaccharide_to_iupac(residue):
    template = "{anomer}-{configuration}-{modification}{base_type}{ring_type}{substituent}"
    anomer = anomer_map_to[residue.anomer]
    configuration = residue.configuration[0].name.upper()
    modification = ""
    base_type = resolve_special_base_type(residue)
    if base_type is None:
        base_type = residue.stem[0].name.title()
    if Modification.d in residue.modifications.values() and "Neu" not in base_type:
        modification = 'd'
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
    multi = False
    for name, pos in get_relevant_substituents(residue):
        part = ""
        if 'n_acetyl' == name:
            part = "NAc"
        elif 'n_glycolyl' == name:
            part = "NGc"
        elif 'sulfate' == name:
            part = 'S'

        # If there is a substituent after the first, successive ones are placed in parentheses
        if multi:
            substituent += "({}{})".format(pos, part)
        else:
            substituent += "{}{}".format(pos, part)
            multi = True
    return substituent


def get_relevant_substituents(residue):
    '''
    Retrieve the set of substituents not implicitly included
    in the base type's symbol name.
    '''
    positions = [p for p, sub in residue.substituents()]
    substituents = [sub.name for p, sub in residue.substituents()]
    if identity.is_a(residue, monosaccharide_reference["HexNAc"], exact=False) or\
       identity.is_a(residue, monosaccharide_reference["NeuAc"], exact=False):
        i = substituents.index("n_acetyl")
        substituents.pop(i)
        positions.pop(i)
    elif identity.is_a(residue, monosaccharide_reference["NeuGc"], exact=False):
        i = substituents.index("n_glycolyl")
        substituents.pop(i)
        positions.pop(i)
    return zip(substituents, positions)


def resolve_special_base_type(residue):
    if residue.superclass == SuperClass.non:
        if residue.stem == (Stem.gro, Stem.gal):
            substituents = [sub.name for p, sub in residue.substituents()]
            if 'n_acetyl' in substituents:
                return "NeuAc"
            elif 'n_glycolyl' in substituents:
                return "NeuGc"
            return "Neu"
    elif residue.superclass == SuperClass.hex:
        if identity.is_a(residue, monosaccharide_reference["HexNAc"], exact=False):
            return residue.stem[0].name.title() + "NAc"

    return None


def glycan_to_iupac(structure=None, attach=''):
    '''
    Translate a |Glycan| structure into IUPAC Three Letter Code.
    Recursively operates on branches.

    Parameters
    ----------
    structure: Glycan or Monosaccharide
        The glycan to be translated. Translation starts from `glycan.root` if `structure`
        is a |Glycan|.
    max_tolerance: int
        The maximum amount of deviance to allow when translating |Monosaccharide| objects
        into nomenclature symbols

    Returns
    -------
    deque
    '''
    base = structure.root if isinstance(structure, Glycan) else structure
    stack = [(attach, base)]
    outstack = deque()
    while(len(stack) > 0):
        outedge_pos, node = stack.pop()
        if outedge_pos in {None, -1}:
            outedge_pos = ''
        link = "({outedge_pos}-{attach})".format(outedge_pos=outedge_pos, attach=attach)
        if attach == "":
            link = ""
        outstack.appendleft('{node}{link}'.format(node=monosaccharide_to_iupac(node), link=link))
        attach = 1
        children = list(node.children())
        if len(children) > 1:
            for pos, child in children[:-1]:
                branch = '({branch})'.format(
                    branch=''.join(glycan_to_iupac(child, pos)),
                    attach_pos=pos
                )
                outstack.appendleft(branch)
            pos, child = children[-1]
            stack.append((pos, child))
        elif len(children) == 1:
            pos, child = children[0]
            stack.append((pos, child))
    return outstack


def to_iupac(structure):
    if isinstance(structure, Monosaccharide):
        return monosaccharide_to_iupac(structure)
    else:
        return ''.join(list(glycan_to_iupac(structure)))
