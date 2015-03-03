import re
from collections import OrderedDict, deque

from pygly2.io import format_constants_map
from pygly2.io.nomenclature import identity
from pygly2.structure import constants, named_structures, Monosaccharide, Glycan, Substituent
from pygly2.utils import invert_dict

Stem = constants.Stem
Configuration = constants.Configuration
monosaccharide_reference = {k: v for k, v in named_structures.monosaccharides.items()}

monosaccharides_to = OrderedDict((
    ("Glc", 'G'),
    ("Gal", 'A'),
    ("GlcNAc", 'GN'),
    ("GalNAc", 'AN'),
    ("Man", "M"),
    ("Neu", "N"),
    ("NeuAc", "NN"),
    ("NeuGc", "NJ"),
    ("KDN", "K"),
    ("Kdo", "W"),
    ("GalA", "L"),
    ("IdoA", "I"),
    ("Rha", "H"),
    ("Fuc", "F"),
    ("Xyl", "X"),
    ("Rib", "B"),
    ("Ara", "R"),
    ("GlcA", "U"),
    ("All", 'O'),
    ("Api", 'P'),
    ("Fru", "E")
))

monosaccharides_from = invert_dict(dict(monosaccharides_to))

substituents_to = {
    'amino': 'Q',
    'ethanolominephosphate': 'PE',
    'inositol': "IN",
    'methyl': "ME",
    'n_acetyl': 'N',
    'o_acetyl': 'T',
    'phosphate': 'P',
    'phosphocholine': "PC",
    'pyruvate': 'PYR',
    'sulfate': 'S',
    'sulfide': 'SH',
    '2-aminoethylphosphonic acid': 'EP'
}

substituents_from = invert_dict(substituents_to)


anomer_map_from = dict(format_constants_map.anomer_map)
anomer_map_to = invert_dict(anomer_map_from)


def get_relevant_substituents(residue):
    positions = [p for p, sub in residue.substituents()]
    substituents = [sub.name for p, sub in residue.substituents()]
    if identity.is_a(residue, monosaccharide_reference["HexNAc"]) or\
       identity.is_a(residue, monosaccharide_reference["NeuAc"]):
        i = substituents.index("n_acetyl")
        substituents.pop(i)
        positions.pop(i)
    elif identity.is_a(residue, "NeuGc"):
        i = substituents.index("n_glycolyl")
        substituents.pop(i)
        positions.pop(i)
    return zip(substituents, positions)

def substituent_to_linearcode(substituent, position=None):
    symbol = None
    try:

        symbol = substituents_to[substituent.name]
    except:
        symbol = substituents_to[substituent]
    if position is not None:
        if position == -1:
            position = ''
        return "{}{}".format(position, symbol)
    else:
        return symbol

def monosaccharide_to_linearcode(monosaccharide, max_tolerance=3):
    tolerance = 0
    while tolerance <= max_tolerance:
        for k, v in monosaccharides_to.items():
            if k not in monosaccharide_reference:
                continue
            if identity.is_a(monosaccharide, monosaccharide_reference[k], tolerance=tolerance):
                residue_sym = v
                substituents_sym =  [(substituent_to_linearcode(*s)) for s in get_relevant_substituents(monosaccharide)]
                if len(substituents_sym) > 0:
                    residue_sym = residue_sym + '[{}]'.format(', '.join(substituents_sym))
                residue_sym = residue_sym + anomer_map_to[monosaccharide.anomer]
                return residue_sym

        tolerance += 1
    raise ValueError("Cannot map {} to LinearCode".format(monosaccharide))


def priority(sym):
    i = 0
    for key, value in monosaccharides_to.items():
        i += 1
        if key == sym or value == sym:
            return i
    return -1


def glycan_to_linearcode(glycan=None, monosaccharide=None, max_tolerance=3):
    stack = [(1, glycan.root if glycan is not None else monosaccharide)]
    outstack = deque()
    while(len(stack) > 0):
        outedge_pos, node = stack.pop()
        if outedge_pos in {1, -1}:
            outedge_pos = ''
        outstack.appendleft(monosaccharide_to_linearcode(node, max_tolerance=max_tolerance) + str(outedge_pos))
        children = []
        for pos, child in node.children():
            rank = priority(child)
            children.append((pos, child, rank))
        if len(children) > 1:
            ordered_children = sorted(children, key=lambda x: x[2])
            for pos, child, rank in ordered_children[:-1]:
                branch = '({branch}{attach_pos})'.format(
                    branch=''.join(glycan_to_linearcode(monosaccharide=child, max_tolerance=max_tolerance)),
                    attach_pos=pos
                    )
                outstack.appendleft(branch)
            pos, child, rank = children[-1]
            stack.append((pos, child))
        elif len(children) == 1:
            pos, child, rank = children[0]
            stack.append((pos, child))
    return outstack



def to_linearcode(structure):
    if isinstance(structure, Monosaccharide):
        return monosaccharide_to_linearcode(structure)
    else:
        return ''.join(list(glycan_to_linearcode(glycan=structure)))


class LinearCodeException(Exception):
    pass


def monosaccharide_from_linearcode(residue_str, parent=None):
    base_type, anomer, outedge, substituents = re.search(r"([A-Z]+)(.)(.)?(\[.*\])?", residue_str).groups()
    base = named_structures.monosaccharides[monosaccharides_from[base_type]]
    base.anomer = anomer_map_from[anomer]
    if substituents is not None:
        for subst_str in substituents[1:-1]:
            pos, name = re.search(r'(\d*)(.+)', subst_str).groups()
            subst_object = Substituent(name)
            try:
                pos = int(pos)
            except:
                pos = -1
            base.add_substituent(subst_object, position=pos)
    try:
        outedge = int(outedge)
    except:
        outedge = -1

    if parent is not None:
        parent.add_monosaccharide(base, position=outedge, child_position=1)

    return base, outedge


def tokenize_linearcode(text):

    last_outedge = None
    root = None
    last_residue = None
    branch_stack = []
    while len(text) > 0:
        # If starting a new branch
        if text[-1] == ')':
            branch_stack.append((last_residue, root, last_outedge))
            root = None
            last_residue = None
            last_outedge = None
            text = text[:-1]
        # If ending a branch
        elif text[-1] == '(':
            try:
                branch_parent, old_root, old_last_outedge = branch_stack.pop()
                branch_parent.add_monosaccharide(root, position=last_outedge, child_position=1)
                root = old_root
                last_residue = branch_parent
                last_outedge = old_last_outedge
                text = text[:-1]
            except IndexError:
                raise LinearCodeException("Bad branching at {}".format(len(text)))
        # Parsing a residue
        else:
            match = re.search(r"([A-Z]+)(.)(.)?(\[.*\])?$", text)
            if match:
                next_residue, outedge = monosaccharide_from_linearcode(text[match.start() : match.end()], last_residue)
                if root is None:
                    last_outedge = outedge
                    root = next_residue
                last_residue = next_residue
                text = text[:match.start()]
            else:
                raise LinearCodeException("Could not identify residue at {}".format(len(text)))

    return Glycan(root)


