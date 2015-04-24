import re
from collections import deque

from pygly2.structure import Monosaccharide, Glycan, constants, named_structures
from pygly2.composition.structure_composition import substituent_compositions
from pygly2.io import format_constants_map
from pygly2.io.nomenclature import identity
from pygly2.utils import invert_dict


class IUPACException(Exception):
    pass


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
    if Modification.d in residue.modifications.values() and\
       ("Neu" not in base_type) and ("Kd" not in base_type):
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


substituents_map_to = {
    name: ''.join(t.title()[:2] for t in name.split("_")) for name in substituent_compositions
}

# Special Cases
substituents_map_to['n_acetyl'] = "NAc"
substituents_map_to['n_glycolyl'] = "NGc"
substituents_map_to['sulfate'] = "S"

substituents_map_from = invert_dict(substituents_map_to)


def resolve_substituent(residue):
    substituent = ""
    multi = False
    for name, pos in get_relevant_substituents(residue):
        if name in substituents_map_to:
            part = substituents_map_to[name]
        else:
            part = [t.title()[:2] for t in name.split("_")]
            substituents_map_to[name] = part
            substituents_map_from[part] = name
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
    if identity.is_a(residue, monosaccharide_reference["NeuAc"], exact=False):
        i = substituents.index("n_acetyl")
        substituents.pop(i)
        positions.pop(i)
    elif identity.is_a(residue, monosaccharide_reference["NeuGc"], exact=False):
        i = substituents.index("n_glycolyl")
        substituents.pop(i)
        positions.pop(i)
    elif identity.is_a(residue, monosaccharide_reference["Neu"], exact=False):
        i = substituents.index("amino")
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
    return None


def glycan_to_iupac(structure=None, attach=None, open_edge='-(', close_edge=')-', open_branch='[', close_branch=']'):
    '''
    Translate a |Glycan| structure into IUPAC Three Letter Code.
    Recursively operates on branches.

    Parameters
    ----------
    structure: Glycan or Monosaccharide
        The glycan to be translated. Translation starts from `glycan.root` if `structure`
        is a |Glycan|.
    attach: int
        The point from the structure tree is attached to its parent. Used for recursively
        handling branches. Defaults to |None|.
    open_edge: str
        The token at the start of an edge. Defaults to the more verbose IUPAC '-('
    close_edge: str
        The token at the end of an edge. Defaults to the more verbose IUPAC ')-'
    open_branch: str
        The token at the start of a branch. Defaults to '['
    close_branch: str
        The token at the end of a branch. Defaults to ']'

    Returns
    -------
    deque
    '''
    base = structure.root if isinstance(structure, Glycan) else structure
    stack = [(1, base)]
    outstack = deque()
    while(len(stack) > 0):
        outedge_pos, node = stack.pop()
        if outedge_pos in {None, -1}:
            outedge_pos = ''
        link = "{oe}{outedge_pos}-{attach}{ce}".format(
            outedge_pos=outedge_pos, attach=attach, oe=open_edge, ce=close_edge)
        if attach is None:
            link = ""
        elif attach != 1 and link[-1] == '-':
            link = link[:-1]
        outstack.appendleft('{node}{link}'.format(node=monosaccharide_to_iupac(node), link=link))
        attach = 1
        children = list(node.children())
        if len(children) > 1:
            for pos, child in children[:-1]:
                branch = '{ob}{branch}{cb}'.format(
                    branch=''.join(glycan_to_iupac(child, pos,
                                                   open_edge=open_edge, close_edge=close_edge,
                                                   open_branch=open_branch, close_branch=close_branch)),
                    attach_pos=pos,
                    ob=open_branch,
                    cb=close_branch
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


monosaccharide_parser = re.compile(r'''(?P<anomer>[abo?])-
                                       (?P<configuration>[LDl])-
                                       (?P<modification>[a-z0-9_]*)
                                       (?P<base_type>[^-]+?)
                                       (?P<ring_type>[pfo?])
                                       (?P<substituent>[^-]*)
                                       (?P<linkage>-\([0-9?]-[0-9?]\)-)?$''', re.VERBOSE)


def monosaccharide_from_iupac(monosaccharide_str, parent=None):
    match = monosaccharide_parser.search(monosaccharide_str)
    if match is None:
        raise IUPACException("Cannot find monosaccharide pattern in {}".format(monosaccharide_str))
    match_dict = match.groupdict()
    anomer = anomer_map_from[match_dict['anomer']]
    base_type = match_dict["base_type"]
    configuration = match_dict["configuration"].lower()
    ring_type = match_dict['ring_type']
    modification = match_dict['modification']
    linkage = [d for d in match_dict.get('linkage') or "" if d.isdigit() or d == "?"]

    residue = named_structures.monosaccharides[base_type]

    if len(residue.configuration) == 1:
        residue.configuration = (configuration,)
    residue.anomer = anomer

    if ring_type == 'p':
        residue.ring_end = residue.ring_start + 4
    elif ring_type == 'f':
        residue.ring_end = residue.ring_start + 3
    elif ring_type == 'o':
        residue.ring_end = residue.ring_start = 0
    else:
        residue.ring_end = residue.ring_start = None

    if modification != "":
        residue.add_modification(modification, 3)
    for position, substituent in substituent_from_iupac(match_dict["substituent"]):
        residue.add_substituent(substituent, position)

    def tryint(i):
        try:
            return int(i)
        except:
            return -1

    linkage = map(tryint, linkage)

    if parent is not None and linkage != ():
        parent.add_monosaccharide(residue, position=linkage[1], child_position=linkage[0])
    return residue, linkage


def substituent_from_iupac(substituents):
    parts = re.split(r"\(|\)", substituents)
    for part in parts:
        if part == "":
            continue
        split_part = re.split(r"(\d+)?", part)
        if len(split_part) == 3:
            _, position, name = split_part
        else:
            position = -1
            name = split_part[0]
        name = (substituents_map_from[name])
        yield int(position), name


def glycan_from_iupac(text):
    last_outedge = None
    root = None
    last_residue = None
    branch_stack = []

    while len(text) > 0:
        print text
        # If starting a new branch
        if text[-1] == ']':
            branch_stack.append((last_residue, root, last_outedge))
            root = None
            last_residue = None
            last_outedge = None
            text = text[:-1]
        # If ending a branch
        elif text[-1] == '[':
            try:
                branch_parent, old_root, old_last_outedge = branch_stack.pop()
                branch_parent.add_monosaccharide(root, position=last_outedge, child_position=1)
                root = old_root
                last_residue = branch_parent
                last_outedge = old_last_outedge
                text = text[:-1]
            except IndexError:
                raise IUPACException("Bad branching at {}".format(len(text)))
        # Parsing a residue
        else:
            match = monosaccharide_parser.search(text)
            if match:
                next_residue, outedge = monosaccharide_from_iupac(text[match.start(): match.end()], last_residue)
                if root is None:
                    last_outedge = outedge
                    root = next_residue
                last_residue = next_residue
                text = text[:match.start()]
            else:
                raise IUPACException("Could not identify residue '...{}' at {}".format(text[-10:], len(text)))

    res = Glycan(root)
    if len(res) > 1:
        return res
    else:
        return res.root
