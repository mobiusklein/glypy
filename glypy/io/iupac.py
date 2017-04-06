import re
from collections import deque

from glypy.structure import Monosaccharide, Glycan, constants, named_structures, Substituent
from glypy.composition.structure_composition import substituent_compositions
from glypy.composition.composition_transform import has_derivatization, derivatize
from glypy.io import format_constants_map
from glypy.io.nomenclature import identity
from glypy.utils import invert_dict

from glypy.io.file_utils import ParserInterface


# A static copy of monosaccharide names to structures for copy-free comparison
monosaccharide_reference = {k: v for k, v in named_structures.monosaccharides.items()}


anomer_map_from = dict(format_constants_map.anomer_map)
anomer_map_from['?'] = anomer_map_from.pop('x')
anomer_map_to = invert_dict(anomer_map_from)
anomer_map_from['beta'] = anomer_map_from['b']
anomer_map_from['alpha'] = anomer_map_from['a']


Stem = constants.Stem
Configuration = constants.Configuration
Modification = constants.Modification
SuperClass = constants.SuperClass


def tryint(i):
    try:
        return int(i)
    except:
        return -1


class IUPACError(Exception):
    pass


def _make_substituent_name(name):
    return ''.join(t.title() for t in name.split("_")).replace("(", "")


substituents_map_to = {
    name: _make_substituent_name(name) for name in substituent_compositions
}

# Special Cases
substituents_map_to['n_acetyl'] = "NAc"
substituents_map_to['n_glycolyl'] = "NGc"
substituents_map_to['sulfate'] = "S"
substituents_map_to["methyl"] = "Me"
substituents_map_to["acetyl"] = "Ac"
substituents_map_to["glycolyl"] = "Gc"
substituents_map_to["fluoro"] = "F"
substituents_map_to["amino"] = "N"

substituents_map_from = invert_dict(substituents_map_to)


class SubstituentSerializer(object):
    def __init__(self, monosaccharides=None):
        if monosaccharides is None:
            monosaccharides = monosaccharide_reference
        self.monosaccharide_reference = monosaccharides

    def __call__(self, residue, monosaccharide_reference=None, **kwargs):
        return self.resolve_substituents(residue, monosaccharide_reference)

    def serialize_substituent(self, substituent):
        name = substituent.name
        if name in substituents_map_to:
            part = substituents_map_to[name]
        else:
            part = _make_substituent_name(name)
            substituents_map_to[name] = part
            substituents_map_from[part] = name
        return part

    def resolve_substituents(self, residue, monosaccharides=None):
        if monosaccharides is None:
            monosaccharides = self.monosaccharide_reference
        substituent = ""
        multi = False
        for name, pos in self.get_relevant_substituents(residue, monosaccharides):
            if pos in {-1, None}:
                pos = ""
            if name in substituents_map_to:
                part = substituents_map_to[name]
            else:
                part = _make_substituent_name(name)
                substituents_map_to[name] = part
                substituents_map_from[part] = name
            # If there is a substituent after the first, successive ones are placed in parentheses
            if multi:
                substituent += "({}{})".format(pos, part)
            else:
                substituent += "{}{}".format(pos, part)
                multi = True
        return substituent

    def get_relevant_substituents(self, residue, monosaccharides=None):
        '''
        Retrieve the set of substituents not implicitly included
        in the base type's symbol name.
        '''
        if monosaccharides is None:
            monosaccharides = self.monosaccharide_reference
        positions = [p for p, sub in residue.substituents() if not sub._derivatize]
        substituents = [sub.name for p, sub in residue.substituents() if not sub._derivatize]
        if identity.is_a(residue, monosaccharides["NeuAc"], exact=False, short_circuit=True):
            try:
                i = substituents.index("n_acetyl")
                substituents.pop(i)
                j = positions.pop(i)
                substituents.insert(i, "acetyl")
                positions.insert(i, j)
            except:  # pragma: no cover
                pass
        elif identity.is_a(residue, monosaccharides["NeuGc"], exact=False, short_circuit=True):
            # i = substituents.index("n_glycolyl")
            # substituents.pop(i)
            # positions.pop(i)
            pass
        elif identity.is_a(residue, monosaccharides["Neu"], exact=False, short_circuit=True):
            i = substituents.index("amino")
            substituents.pop(i)
            positions.pop(i)

        return zip(substituents, positions)


resolve_substituent = SubstituentSerializer()


class ModificationSerializer(object):
    def extract_modifications(self, modifications, base_type):
        buff = []
        template = '{position}-{name}'
        pos_mod_pairs = list(modifications.items())
        try:
            pos, mods = map(list, zip(*pos_mod_pairs))
        except ValueError:
            pos, mods = [], []
        if "Neu" in base_type or "Kd" in base_type:
            for mod in [Modification.d, Modification.keto, Modification.a]:
                try:
                    pop_ix = mods.index(mod)
                    pos.pop(pop_ix)
                    mods.pop(pop_ix)
                except:  # pragma: no cover
                    pass

        elif "Fuc" in base_type:
            for mod in [Modification.d]:
                pop_ix = mods.index(mod)
                pos.pop(pop_ix)
                mods.pop(pop_ix)

        pos_mod_pairs = zip(pos, mods)
        for pos, mod in pos_mod_pairs:
            if pos != -1:
                buff.append(template.format(position=pos, name=mod.name))
            else:
                buff.append(mod.name)
        return ','.join(buff)

    def __call__(self, modifications, base_type):
        return self.extract_modifications(modifications, base_type)


extract_modifications = ModificationSerializer()


class ModificationDeserializer(object):
    def parse_modifications(self, modification_string):
        buff = modification_string.split(",")
        pairs = []
        for token in buff:
            if token == '':
                continue
            try:
                pos, mod = token.split("-")
            except:
                pos = -1
                mod = token
            pairs.append((int(pos), Modification[mod]))
        return pairs

    def __call__(self, modification_string):
        return self.parse_modifications(modification_string)


parse_modifications = ModificationDeserializer()


class MonosaccharideSerializer(object):
    def __init__(self, monosaccharides=None, substituent_resolver=None, modification_extractor=None):
        if monosaccharides is None:
            monosaccharides = monosaccharide_reference
        self.monosaccharide_reference = monosaccharides
        if substituent_resolver is None:
            substituent_resolver = SubstituentSerializer(monosaccharides)
        self.substituent_resolver = substituent_resolver
        if modification_extractor is None:
            modification_extractor = ModificationSerializer()
        self.modification_extractor = modification_extractor

    def resolve_special_base_type(self, residue, monosaccharide_reference=None):
        if monosaccharide_reference is None:
            monosaccharide_reference = self.monosaccharide_reference
        if residue.superclass == SuperClass.non:
            if residue.stem == (Stem.gro, Stem.gal):
                substituents = [sub.name for p, sub in residue.substituents() if not sub._derivatize]
                modifications = [mod for p, mod in residue.modifications.items()]
                if Modification.a in modifications and\
                   Modification.keto in modifications and\
                   Modification.d in modifications:
                    if len(substituents) == 0:
                        return "Kdn"
                    elif "n_acetyl" in substituents:
                        return "Neu"  # Ac
                    elif "n_glycolyl" in substituents:
                        return "Neu"  # Gc
                    elif "amino" in substituents:
                        return "Neu"  # _

        elif residue.superclass == SuperClass.oct:
            if residue.stem == (Stem.man,):
                if Modification.a in residue.modifications[1] and\
                   Modification.keto in residue.modifications[2] and\
                   Modification.d in residue.modifications[3]:
                    return "Kdo"
        elif residue.stem == (Stem.gal,) and residue.configuration == (Configuration.l,):
            if Modification.d in residue.modifications.values():
                return "Fuc"

        return None

    def monosaccharide_to_iupac(self, residue):
        template = "{anomer}-{configuration}-{modification}{base_type}{ring_type}{substituent}"
        anomer = anomer_map_to[residue.anomer]
        if residue.configuration[0] is Configuration.Unknown:
            configuration = "?"
        else:
            configuration = residue.configuration[0].name.upper()
        modification = ""
        base_type = self.resolve_special_base_type(residue)
        if base_type is None:
            if len(residue.stem) == 1 and residue.stem[0] is not Stem.Unknown:
                base_type = residue.stem[0].name.title()
            else:
                base_type = residue.superclass.name.title()
        modification = self.modification_extractor(residue.modifications, base_type)
        ring_type = residue.ring_type.name[0]
        substituent = self.substituent_resolver(residue)
        return template.format(
            anomer=anomer,
            configuration=configuration,
            modification=modification,
            base_type=base_type,
            ring_type=ring_type,
            substituent=substituent
        )

    def __call__(self, residue):
        return self.monosaccharide_to_iupac(residue)


class DerivatizationAwareMonosaccharideSerializer(MonosaccharideSerializer):
    def monosaccharide_to_iupac(self, residue):
        string = super(DerivatizationAwareMonosaccharideSerializer, self).monosaccharide_to_iupac(residue)
        deriv = has_derivatization(residue)
        if deriv:
            string = "%s^%s" % (string, self.substituent_resolver.serialize_substituent(deriv))
        return string


monosaccharide_to_iupac = MonosaccharideSerializer()
resolve_special_base_type = monosaccharide_to_iupac.resolve_special_base_type


class GlycanSerializer(object):
    def __init__(self, monosaccharide_serializer=None, open_edge='-(', close_edge=')-',
                 open_branch='[', close_branch=']'):
        if monosaccharide_serializer is None:
            monosaccharide_serializer = MonosaccharideSerializer()
        self.monosaccharide_serializer = monosaccharide_serializer

        self.open_edge = open_edge
        self.close_edge = close_edge
        self.open_branch = open_branch
        self.close_branch = close_branch

    def glycan_to_iupac(self, structure=None, attach=None, is_branch=False):
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

        Returns
        -------
        deque
        '''
        base = structure.root if isinstance(structure, Glycan) else structure
        stack = [(attach, base)]
        outstack = deque()
        while(len(stack) > 0):
            outedge, node = stack.pop()
            link = ""
            if outedge is not None:
                link = "{oe}{attach}-{outedge_pos}{ce}".format(
                    outedge_pos=outedge.parent_position,
                    attach=outedge.child_position,
                    oe=self.open_edge, ce=self.close_edge)
            # Branch linkage does not start with leading dash
            if is_branch and link[-1] == '-':
                link = link[:-1]
            outstack.appendleft('{node}{link}'.format(node=self.monosaccharide_serializer(node), link=link))
            # Reset for next pass through the loop
            is_branch = False
            children = list((p, link) for p, link in node.links.items() if link.is_parent(node))
            if len(children) > 1:
                for pos, link in children[:-1]:
                    branch = '{ob}{branch}{cb}'.format(
                        branch=''.join(self.glycan_to_iupac(link.child, link, is_branch=True)),
                        ob=self.open_branch,
                        cb=self.close_branch
                    )
                    outstack.appendleft(branch)
                pos, link = children[-1]
                stack.append((link, link.child))
            elif len(children) == 1:
                pos, link = children[0]
                stack.append((link, link.child))
        return outstack

    def __call__(self, structure):
        return ''.join(self.glycan_to_iupac(structure))


glycan_to_iupac = GlycanSerializer()


def to_iupac(structure):
    """Translate `structure` into its textual representation using IUPAC Three Letter Code

    Parameters
    ----------
    structure : |Glycan| or |Monosaccharide|
        The structure to be translated

    Returns
    -------
    |str|
    """
    if isinstance(structure, Monosaccharide):
        return monosaccharide_to_iupac(structure)
    else:
        return glycan_to_iupac(structure)


def aminate_substituent(substituent):
    if substituent.name.startswith("n_"):
        # already aminated
        return substituent
    aminated = Substituent("n_" + substituent.name)
    if aminated.composition == {}:
        raise ValueError("Could not aminate substituent")
    return aminated


monosaccharide_parser = re.compile(r'''(?P<anomer>[abo?])-
                                       (?P<configuration>[LD?])-
                                       (?P<modification>[a-z0-9_\-,]*)
                                       (?P<base_type>[^-]{3}?)
                                       (?P<ring_type>[xpfo?])
                                       (?P<substituent>[^-]*?)
                                       (?P<linkage>-\([0-9?]-[0-9?]\)-?)?$''', re.VERBOSE)


class SubstituentDeserializer(object):
    def substituent_from_iupac(self, substituents):
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
                    import warnings
                    warnings.warn("No translation rule found to convert %s into a Substituent" % name)
                    continue
            yield int(position), name

    def symbol_to_name(self, symbol):
        name = substituents_map_from[symbol]
        return name

    def __call__(self, substituents):
        return self.substituent_from_iupac(substituents)


substituent_from_iupac = SubstituentDeserializer()


class MonosaccharideDeserializer(object):
    pattern = re.compile(r'''(?:(?P<anomer>[abo?]|alpha|beta)-)?
                         (?P<configuration>[LD?])-
                         (?P<modification>[a-z0-9_\-,]*)
                         (?P<base_type>[^-]{3}?)
                         (?P<ring_type>[xpfo?])?
                         (?P<substituent>[^-]*?)
                         (?P<linkage>-\([0-9?]->?[0-9?]\)-?)?
                         $''', re.VERBOSE)

    def __init__(self, modification_parser=None, substituent_parser=None):
        if modification_parser is None:
            modification_parser = ModificationDeserializer()
        self.modification_parser = modification_parser
        if substituent_parser is None:
            substituent_parser = SubstituentDeserializer()
        self.substituent_parser = substituent_parser

    def has_pattern(self, string):
        return self.pattern.search(string)

    def extract_pattern(self, monosaccharide_str):
        match = self.pattern.search(monosaccharide_str)
        if match is None:
            raise IUPACError("Cannot find monosaccharide pattern in {}".format(monosaccharide_str))
        match_dict = match.groupdict()
        return match_dict

    def ring_bounds(self, residue, ring_type):
        if ring_type == 'p':
            residue.ring_end = residue.ring_start + 4
        elif ring_type == 'f':
            residue.ring_end = residue.ring_start + 3
        elif ring_type == 'o':
            residue.ring_end = residue.ring_start = 0
        else:
            residue.ring_end = residue.ring_start = None

    def build_residue(self, match_dict):
        try:
            anomer = anomer_map_from[match_dict['anomer']]
        except KeyError:
            anomer = anomer_map_from['?']
        base_type = match_dict["base_type"]
        configuration = match_dict["configuration"].lower()
        ring_type = match_dict['ring_type']

        modification = match_dict['modification']

        linkage = [d for d in match_dict.get('linkage') or "" if d.isdigit() or d == "?"]

        residue = named_structures.monosaccharides[base_type]
        base_is_modified = len(residue.substituent_links) + len(residue.modifications) > 0

        if len(residue.configuration) == 1:
            residue.configuration = (configuration,)

        residue.anomer = anomer
        self.ring_bounds(residue, ring_type)
        self.set_modifications(residue, modification)
        self.set_substituents(residue, match_dict['substituent'], base_is_modified, base_type)
        return residue, linkage

    def set_substituents(self, residue, substituent_string, base_is_modified, base_type):
        i = 0
        for position, substituent in self.substituent_parser(substituent_string):
            i += 1
            if position == -1 and base_is_modified:
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
                        raise IUPACError("Can't resolve %s" % substituent)
                else:
                    raise

    def set_modifications(self, residue, modification_string):
        for pos, mod in self.modification_parser(modification_string):
            residue.add_modification(mod, pos)

    def monosaccharide_from_iupac(self, monosaccharide_str, parent=None):
        match_dict = self.extract_pattern(monosaccharide_str)
        residue, linkage = self.build_residue(match_dict)
        linkage = list(map(tryint, linkage))

        self.add_monosaccharide_bond(residue, parent, linkage)
        return residue, linkage

    def add_monosaccharide_bond(self, residue, parent, linkage):
        if parent is not None and linkage != ():
            parent.add_monosaccharide(
                residue, position=linkage[1], child_position=linkage[0])

    def __call__(self, monosaccharide_str, parent=None):
        return self.monosaccharide_from_iupac(monosaccharide_str, parent=parent)

    def finalize(self, glycan):
        pass


class DerivatizationAwareMonosaccharideDeserializer(MonosaccharideDeserializer):
    pattern = re.compile(r'''(?:(?P<anomer>[abo?]|alpha|beta)-)?
                             (?P<configuration>[LD?])-
                             (?P<modification>[a-z0-9_\-,]*)
                             (?P<base_type>[^-]{3}?)
                             (?P<ring_type>[xpfo?])?
                             (?P<substituent>[^-]*?)
                             (?P<derivatization>\^[^\s-]*?)?
                             (?P<linkage>-\([0-9?]->?[0-9?]\)-?)?$''', re.VERBOSE)

    def add_monosaccharide_bond(self, residue, parent, linkage):
        if parent is not None and linkage != ():
            try:
                parent.add_monosaccharide(residue, position=linkage[1], child_position=linkage[0])
            except ValueError:
                parent_substituent_links_at_site = parent.substituent_links[linkage[1]]
                if (parent_substituent_links_at_site and parent_substituent_links_at_site[0].child._derivatize):
                    parent.drop_substituent(linkage[1], parent_substituent_links_at_site[0].child)
                residue_substituent_links_at_site = residue.substituent_links[linkage[0]]
                if residue_substituent_links_at_site and residue_substituent_links_at_site[0].child._derivatize:
                    residue.drop_substituent(linkage[0], residue_substituent_links_at_site[0].child)

                parent.add_monosaccharide(residue, position=linkage[1], child_position=linkage[0])

    def apply_derivatization(self, residue, deriv):
        if deriv.startswith("^"):
            deriv = deriv[1:]
            deriv = self.substituent_parser.symbol_to_name(deriv)
            derivatize(residue, deriv)
        else:
            raise IUPACError("Derivatization Extension Must Start with '^'")

    def monosaccharide_from_iupac(self, monosaccharide_str, parent=None):
        match_dict = self.extract_pattern(monosaccharide_str)
        residue, linkage = self.build_residue(match_dict)
        linkage = list(map(tryint, linkage))

        self.add_monosaccharide_bond(residue, parent, linkage)

        deriv = match_dict.get("derivatization", '')
        if deriv is not None and deriv != "":
            self.apply_derivatization(residue, deriv)

        return residue, linkage

    def finalize(self, glycan):
        for node in glycan:
            neg_capacity = -node._remaining_capacity()
            if neg_capacity > 0:
                unknowns = node.substituent_links[-1]
                to_remove = []
                for unknown in unknowns:
                    if unknown.child.node_type is Substituent.node_type and unknown.child._derivatize:
                        if neg_capacity > 0:
                            to_remove.append(unknown)
                            neg_capacity -= 1
                        else:
                            break
                for link_to_remove in to_remove:
                    link_to_remove.break_link(refund=True)
                if neg_capacity > 0:
                    raise ValueError("Could not completely remove overload from %s" % (node,))


monosaccharide_from_iupac = MonosaccharideDeserializer()


class GlycanDeserializer(object):
    def __init__(self, monosaccharide_deserializer=None):
        if monosaccharide_deserializer is None:
            monosaccharide_deserializer = MonosaccharideDeserializer()
        self.monosaccharide_deserializer = monosaccharide_deserializer

    def add_monosaccharide(self, parent_node, child_node, parent_position, child_position):
        # parent_node.add_monosaccharide(
        #     child_node, position=parent_position, child_position=child_position)
        self.monosaccharide_deserializer.add_monosaccharide_bond(
            child_node, parent_node, (child_position, parent_position))

    def glycan_from_iupac(self, text, **kwargs):
        last_outedge = None
        root = None
        last_residue = None
        branch_stack = []

        # Remove the base
        text = re.sub(r"\(\?->?$", "", text)

        while len(text) > 0:

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
                    child_position, parent_position = last_outedge
                    self.add_monosaccharide(branch_parent, root, parent_position, child_position)
                    root = old_root
                    last_residue = branch_parent
                    last_outedge = old_last_outedge
                    text = text[:-1]
                except IndexError:
                    raise IUPACError("Bad branching at {}".format(len(text)))
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
                    raise IUPACError("Could not identify residue '...{}' at {}".format(text[-30:], len(text)))
        res = Glycan(root)
        self.monosaccharide_deserializer.finalize(res)
        return res

    def __call__(self, text, **kwargs):
        return self.glycan_from_iupac(text, **kwargs)


glycan_from_iupac = GlycanDeserializer()


def from_iupac(text, **kwargs):
    """Parse the given text into an instance of |Glycan|. If there is only a single monosaccharide
    in the output, just the Monosaccharide instance is returned.

    Parameters
    ----------
    text : |str|

    Returns
    -------
    |Glycan| or |Monosaccharide|
        If the resulting structure is just a single monosaccharide, the returned value is a Monosaccharide.
    """
    res = glycan_from_iupac(text, **kwargs)
    if len(res) > 1:
        return res
    else:
        return res.root


loads = from_iupac
dumps = to_iupac


class IUPACParser(ParserInterface):
    def process_result(self, line):
        structure = loads(line)
        return structure


Monosaccharide.register_serializer("iupac", dumps)
Glycan.register_serializer("iupac", dumps)
