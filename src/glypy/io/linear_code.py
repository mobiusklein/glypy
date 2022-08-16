'''
A module for operating on GlycoMinds Linear Code

Assumes that the structure's root is the right-most residue, as shown in
:title-reference:`A Novel Linear Code Nomenclature for Complex Carbohydrates, Banin et al.`.

Currently does not handle the sigils indicating deviation from the common forms.
'''

import re
from collections import OrderedDict, deque

from six import string_types as basestring

from glypy.io import format_constants_map
from glypy.io.nomenclature import identity
from glypy.structure import constants, named_structures, Monosaccharide, Glycan, Substituent
from glypy.utils import invert_dict

from glypy.io.file_utils import ParserInterface, ParserError

Stem = constants.Stem
Configuration = constants.Configuration

# A static copy of monosaccharide names to structures for copy-free comparison
monosaccharide_reference = {k: v.clone() for k, v in named_structures.monosaccharides.items()}
# Unset the anomericity as this does not influence monosaccharide resolution
for k, v in monosaccharide_reference.items():
    v.anomer = None

#: A mapping from common monosaccharide names to their symbol, ordered by priority
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
    # ("Api", 'P'),
    ("Fru", "E")
))

#: A mapping from symbol to common monosaccharide name
monosaccharides_from = invert_dict(dict(monosaccharides_to))

#: A mapping from common substituent names to symbol
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

#: A mapping from symbol to common substituent names
substituents_from = invert_dict(substituents_to)


anomer_map_from = dict(format_constants_map.anomer_map)
anomer_map_from['?'] = anomer_map_from.pop('x')
anomer_map_to = invert_dict(anomer_map_from)


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


def substituent_to_linear_code(substituent, position=None):
    '''
    Translate a |Substituent| to Linear Code. Include's the substituent's
    position if it is known.

    Parameters
    ----------
    substituents: Substituent or str
        The structure or the name to translate
    position:
        The position of the structure to try including

    Returns
    -------
    str

    Raises
    ------
    KeyError:
        When an unknown symbol is encountered
    '''
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


def monosaccharide_to_linear_code(monosaccharide, max_tolerance=3):
    '''
    Perform iteratively permissive attempts to translate `monosaccharide` into
    a nomenclature symbol.

    .. note::
        Uses a multi-pass approach. Could alternatively do a single pass
        and keep the best match.

    Parameters
    ----------
    monosaccharide: Monosaccharide
        The residue to be translated
    max_tolerance: int
        The maximum error tolerance to allow while looking for a match

    Returns
    -------
    str

    Raises
    ------
    ValueError:
        When no suitable translation can be found
    KeyError:
        When an unknown symbol is encountered
    '''
    tolerance = 0
    if identity.is_generic_monosaccharide(monosaccharide):
        raise LinearCodeError("Linear Code does not support generic monosaccharide %s" % str(monosaccharide))
    while tolerance <= max_tolerance:
        for k, v in monosaccharides_to.items():
            if k not in monosaccharide_reference:
                continue
            if identity.is_a(monosaccharide, monosaccharide_reference[k], tolerance=tolerance):
                residue_sym = v
                substituents_sym = [(
                    substituent_to_linear_code(
                        *s)) for s in get_relevant_substituents(monosaccharide)]
                if len(substituents_sym) > 0:
                    residue_sym = residue_sym + '[{}]'.format(', '.join(substituents_sym))
                residue_sym = residue_sym + anomer_map_to[monosaccharide.anomer]
                return residue_sym

        tolerance += 1
    raise LinearCodeError("Cannot map {} to Linear Code".format(monosaccharide))


def priority(sym):
    '''
    Calculate the branching priority for a given symbol or |Monosaccharide|.
    Used when deciding when a residue is considered branch or backbone.

    Parameters
    ----------
    sym: str or Monosaccharide

    Returns
    -------
    int
    '''
    i = 0
    is_str = isinstance(sym, basestring)
    for key, value in monosaccharides_to.items():
        i += 1
        if (identity.is_a(sym, monosaccharide_reference[key])) if not is_str else value == sym:
            return i
    return -1


def glycan_to_linear_code(structure=None, max_tolerance=3):
    '''
    Translate a |Glycan| structure into Linear Code. Called from :func:`to_linear_code`.
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
    stack = [(1, base)]
    outstack = deque()
    while(len(stack) > 0):
        outedge_pos, node = stack.pop()
        if outedge_pos in {1, -1}:
            outedge_pos = ''
        outstack.appendleft(monosaccharide_to_linear_code(node, max_tolerance=max_tolerance) + str(outedge_pos))
        children = []
        for pos, child in node.children():
            rank = priority(child)
            children.append((pos, child, rank))
        if len(children) > 1:
            ordered_children = sorted(children, key=lambda x: x[2])
            for pos, child, rank in ordered_children[:-1]:
                branch = '({branch}{attach_pos})'.format(
                    branch=''.join(glycan_to_linear_code(child, max_tolerance=max_tolerance)),
                    attach_pos=pos
                )
                outstack.appendleft(branch)
            pos, child, rank = ordered_children[-1]
            stack.append((pos, child))
        elif len(children) == 1:
            pos, child, rank = children[0]
            stack.append((pos, child))
    return outstack


def to_linear_code(structure):
    '''
    Translates `structure` to Linear Code.

    Parameters
    ----------
    structure: Monosaccharide or Glycan

    Returns
    -------
    str
    '''
    if isinstance(structure, Monosaccharide):
        return monosaccharide_to_linear_code(structure)
    else:
        return ''.join(list(glycan_to_linear_code(structure)))


class LinearCodeError(ParserError):
    pass


def monosaccharide_from_linear_code(residue_str, parent=None):
    '''
    Helper function for :func:`parse_linear_code`. Given a residue string
    of the form "(base_type)(anomer)(outedge?)([substituents]?)", construct
    a |Monosaccharide| object. If `parent` is not |None|, connect the
    resulting |Monosaccharide| to `parent` at `outedge`.
    '''
    base_type, substituents, anomer, outedge = re.search(r"([A-Z]+)(\[.*?\])?([abo\?]?)(.)?", residue_str).groups()
    base = named_structures.monosaccharides[monosaccharides_from[base_type]]
    base.anomer = anomer_map_from.get(anomer, None)
    if substituents is not None:
        for subst_str in substituents[1:-1].split(','):
            pos, name = re.search(r'(\d*)(.+)', subst_str).groups()
            subst_object = Substituent(substituents_from[name])
            try:
                pos = int(pos)
            except (ValueError, TypeError):
                pos = -1
            base.add_substituent(subst_object, position=pos)
    try:
        outedge = int(outedge)
    except (ValueError, TypeError):
        outedge = -1

    if parent is not None:
        parent.add_monosaccharide(base, position=outedge, child_position=min(base.open_attachment_sites()[0]))

    return base, outedge


def parse_linear_code(text, structure_class=Glycan):
    '''
    Parse the character string `text`, extracting GlycoMinds Linear Code-format
    carbohydrate structures, converting them into a |Glycan| object.

    .. note:: Assumes that the structure's root is the right-most residue

    Supports only *concrete* structures.

    The resulting structure will be canonicalized.

    Parameters
    ----------
    text: str
        The string to be parsed

    Returns
    -------
    Glycan

    Raises
    ------
    LinearCodeError:
        When an error is encountered while parsing resulting from malformed syntax
        or structure

    ValueError:
        When a symbol is encountered for which no suitable translation could be found

    KeyError:
        When an unknown symbol is encountered
    '''
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
                branch_parent.add_monosaccharide(root, position=last_outedge,
                                                 child_position=min(root.open_attachment_sites()[0]))
                root = old_root
                last_residue = branch_parent
                last_outedge = old_last_outedge
                text = text[:-1]
            except IndexError:
                raise LinearCodeError("Bad branching at {}".format(len(text)))
        # Parsing a residue
        else:
            match = re.search(r"([A-Z]+?)(\[[^\]]*?\])?(.)(.)?$", text)
            if match:
                next_residue, outedge = monosaccharide_from_linear_code(text[match.start(): match.end()], last_residue)
                if root is None:
                    last_outedge = outedge
                    root = next_residue
                last_residue = next_residue
                text = text[:match.start()]
            else:
                raise LinearCodeError("Could not identify residue '...{}' at {}".format(text[-10:], len(text)))

    res = structure_class(root=root).reindex()
    res.canonicalize()
    if len(res) > 1:
        return res
    else:
        return res.root


#: Common alias for :func:`to_linear_code`
dumps = to_linear_code

#: Common alias for :func:`parse_linear_code`
loads = parse_linear_code


class LinearCodeParser(ParserInterface):
    def process_result(self, line):
        structure = loads(line)
        return structure


Monosaccharide.register_serializer('linear_code', dumps)
Glycan.register_serializer('linear_code', dumps)
