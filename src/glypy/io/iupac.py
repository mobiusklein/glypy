import re
import warnings

from collections import deque, namedtuple
from functools import partial

from glypy.structure import (
    Monosaccharide, Glycan, Link, AmbiguousLink,
    Substituent, constants, named_structures, UnknownPosition)
from glypy.composition import Composition
from glypy.composition.structure_composition import substituent_compositions
from glypy.composition.composition_transform import has_derivatization, derivatize
from glypy.io import format_constants_map
from glypy.io.nomenclature import identity
from glypy.utils import invert_dict

from glypy.io.file_utils import ParserInterface, ParserError


# A static copy of monosaccharide names to structures for copy-free comparison
monosaccharide_reference = {k: v for k, v in named_structures.monosaccharides.items()}

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

special_base_type_resolver = identity.MonosaccharideIdentifier(special_base_types)

anomer_map_from = dict(format_constants_map.anomer_map)
anomer_map_from['?'] = anomer_map_from.pop('x')
anomer_map_to = invert_dict(anomer_map_from)
anomer_map_from['beta'] = anomer_map_from['b']
anomer_map_from['alpha'] = anomer_map_from['a']
anomer_map_from[u"\u03B1"] = anomer_map_from['a']
anomer_map_from[u"\u03B2"] = anomer_map_from['b']


Stem = constants.Stem
Configuration = constants.Configuration
Modification = constants.Modification
SuperClass = constants.SuperClass


def tryint(i):
    try:
        return int(i)
    except ValueError:
        return -1


class IUPACError(ParserError):
    pass


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


class SubstituentSerializer(object):
    """Build the textual encoding for the relevant substituents for
    a provided monosaccharide.

    Attributes
    ----------
    monosaccharide_reference : :class:`dict`
        Map base type to :class:`~.Monosaccharide`
    """

    def __init__(self, monosaccharides=None, substitution_rules=None, substituent_map=None):
        if monosaccharides is None:
            monosaccharides = monosaccharide_reference
        if substitution_rules is None:
            substitution_rules = _substituent_replacement_rules
        if substituent_map is None:
            substituent_map = substituents_map_to
        self.monosaccharide_reference = monosaccharides
        self.substitution_rules = substitution_rules
        self.substituent_map = substituent_map

    def __call__(self, residue, **kwargs):
        """Alias for :meth:`resolve_substituents`
        """
        return self.resolve_substituents(residue, **kwargs)

    def serialize_substituent(self, substituent):
        """Obtain an IUPAC-compatible name for ``substituent``

        Parameters
        ----------
        substituent : :class:`~.Substituent`
            The subsituent group to get the name for

        Returns
        -------
        :class:`str`
        """
        name = substituent.name
        if name in self.substituent_map:
            part = self.substituent_map[name]
        else:
            part = _make_substituent_name(name)
            warnings.warn("Registering IUPAC  name %r for %r" % (name, part))
            self.substituent_map[name] = part
            if part not in substituents_map_from:
                substituents_map_from[part] = name
        return part

    def resolve_substituents(self, residue, **kwargs):
        """Build a textual encoding of the substituent list for ``residue``.

        Parameters
        ----------
        residue : :class:`~.Monosaccharide`
            The residue to build the substituent list for

        Returns
        -------
        :class:`str`
        """
        substituent = ""
        multi = False
        for name, pos in self.get_relevant_substituents(residue):
            if pos in {UnknownPosition, None}:
                pos = ""
            if name in self.substituent_map:
                part = self.substituent_map[name]
            else:
                part = _make_substituent_name(name)
                warnings.warn("Registering IUPAC  name %r for %r" %
                              (name, part))
                self.substituent_map[name] = part
                if part not in substituents_map_from:
                    substituents_map_from[part] = name
            # If there is a substituent after the first, successive ones are placed in parentheses
            if multi:
                substituent += "({}{})".format(pos, part)
            else:
                substituent += "{}{}".format(pos, part)
                multi = True
        return substituent

    def _test_for_replacement(self, residue, reference, positions, substituents, replacements, exact=False):
        if identity.is_a(residue, reference, exact=exact, short_circuit=True):
            self._substituent_replacement(positions, substituents, replacements)

    def _substituent_replacement(self, positions, substituents, pairs):
        for target, replacement in pairs:
            try:
                i = substituents.index(target)
                substituents.pop(i)
                j = positions.pop(i)
                if replacement is not None:
                    substituents.insert(i, replacement)
                    positions.insert(i, j)
            except Exception:  # pragma: no cover
                pass

    def get_relevant_substituents(self, residue):
        '''
        Retrieve the set of substituents not implicitly included
        in the base type's symbol name.

        Certain base types have implied substituent groups or partial substituent
        groups. For example, from the perspective of :mod:`glypy`, "n-acetyl" is
        a discrete unit, but from a structural perspective it is a substituted amine
        group that was later acetylated. The "Neu" base type implies an amination of
        carbon 5. In "Neu5Ac", the amine at carbon 5 is acetylated, but because the
        amine is implied, the N of the "NAc" signifier is omitted.

        In IUPAC's trivial coding, substituents are listed following
        the base type, with each substituent encoded as an optional position
        specifier followed immediately by a shortened version of the substituent group's
        name. The first substituent immediately follows the base type, and subsequent
        substituent groups are enclosed in parentheses.
        '''
        monosaccharides = self.monosaccharide_reference
        positions = [p for p, sub in residue.substituents() if not sub._derivatize]
        substituents = [sub.name for p, sub in residue.substituents() if not sub._derivatize]

        for reference_name, replacements in self.substitution_rules.items():
            reference = monosaccharides[reference_name]
            self._test_for_replacement(residue, reference, positions, substituents, replacements)
        return zip(substituents, positions)


resolve_substituent = SubstituentSerializer()


class ModificationSerializer(object):
    """Build the textual encoding for the relevant modifications for
    a provided monosaccharide base type.
    """

    def extract_modifications(self, modifications, base_type):
        """Build a string representing the relevant modifications for.

        Certain base types imply a collection of modifications by default. These
        modifications should not be included in the textual encoding.

        In IUPAC's trivial coding, modifications are specified as a comma-separated
        list preceding the base definition with an optional position designation linked
        by a dash character to the modification name. For example "deoxy" and "?-deoxy"
        both signify a deoxidation at an unknown position, and "6-deoxy" indicates the
        deoxidation appears at carbon 6.

        Parameters
        ----------
        modifications : :class:`~.MultiMap`
            A mapping between position and modifications
        base_type : :class:`str`
            The monosaccharide base type to build the list for. Uses
            a hard-coded list of modifications

        Returns
        -------
        :class:`str`
        """
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
                except Exception:  # pragma: no cover
                    pass

        elif "Fuc" in base_type or "Qui" in base_type or "Rha" in base_type:
            for mod in [Modification.d]:
                pop_ix = mods.index(mod)
                pos.pop(pop_ix)
                mods.pop(pop_ix)
        elif "Fru" in base_type or "Psi" in base_type or "Sor" in base_type or "Tag" in base_type:
            for mod in [Modification.keto]:
                pop_ix = mods.index(mod)
                pos.pop(pop_ix)
                mods.pop(pop_ix)
        pos_mod_pairs = zip(pos, mods)
        for pos, mod in pos_mod_pairs:
            if pos != UnknownPosition:
                buff.append(template.format(position=pos, name=mod.name))
            else:
                buff.append(mod.name)
        out = ','.join(buff)
        if out:
            out += '-'
        return out

    def __call__(self, modifications, base_type):
        """An alias for :meth:`extract_modifications`
        """
        return self.extract_modifications(modifications, base_type)


extract_modifications = ModificationSerializer()


class ModificationDeserializer(object):
    """Parses modification signifiers from text into position, :class:`~.Modification` pairs

    Attributes
    ----------
    modification_map : :class:`dict`
        Mapping from text representation to :class:`~.Modification` to provide additional
        names for the existing modification name mapping.
    """

    def __init__(self, modification_map=None):
        if modification_map is None:
            modification_map = _modification_map_to.copy()
        else:
            t = _modification_map_to.copy()
            t.update(modification_map)
            modification_map = t
        self.modification_map = modification_map

    def parse_modifications(self, modification_string):
        """Parses the text for site-modification definitions.

        In IUPAC's trivial coding, modifications are specified as a comma-separated
        list preceding the base definition with an optional position designation linked
        by a dash character to the modification name. For example "deoxy" and "?-deoxy"
        both signify a deoxidation at an unknown position, and "6-deoxy" indicates the
        deoxidation appears at carbon 6.

        Parameters
        ----------
        modification_string : str

        Returns
        -------
        list:
            The list of (position, modification) pairs parsed from the string

        Raises
        ------
        IUPACError:
            If the modification signifier cannot be translated
        """
        buff = modification_string.split(",")
        pairs = []
        for token in buff:
            if token == '':
                continue
            try:
                pos, mod = token.split("-")
            except Exception:
                pos = UnknownPosition
                mod = token
            try:
                mod_t = self.modification_map.get(mod, mod)
                pos = int(pos)
                pairs.append((pos, Modification[mod_t]))
            except KeyError:
                raise IUPACError("Could not determine modification from %s" % modification_string)
        return pairs

    def __call__(self, modification_string):
        """An alias for :meth:`parse_modifications`
        """
        return self.parse_modifications(modification_string)


parse_modifications = ModificationDeserializer()


class MonosaccharideSerializer(object):
    """Serialize a :class:`~.Monosaccharide` object to IUPAC text

    Attributes
    ----------
    modification_extractor: :class:`ModificationSerializer`
        Convert modifications to a text list
    monosaccharide_reference : :class:`dict`
        Map base type to :class:`~.Monosaccharide`
    substituent_resolver : :class:`SubstituentSerializer`
        Convert substituents to a text list
    """

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

    def resolve_special_base_type(self, residue):
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
        elif residue.stem == (Stem.gal,):
            if Modification.d in residue.modifications.values():
                return "Fuc"
        elif residue.stem == (Stem.man,):
            if Modification.d in residue.modifications.values():
                return "Rha"
        elif residue.stem == (Stem.glc,):
            if Modification.d in residue.modifications.values():
                return "Qui"
        query = special_base_type_resolver.query(residue)
        if query:
            return special_base_type_resolver.name_map[query]
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
    """A derivatization aware version of :class:`MonosaccharideSerializer` which
    deviates from the standard IUPAC code to encode derivatization.

    If a :class:`~.Monosaccharide` object has a derivatizing substituent attached to
    it, as detected by :func:`~.has_derivatization`, those substituent groups will
    normally be ignored. With this subclass, a single entry will be appended to the
    monosaccharide encoding joined by an "^" character. For example a permethylated
    hexose would be written "Hex^Me".
    """

    def monosaccharide_to_iupac(self, residue):
        string = super(DerivatizationAwareMonosaccharideSerializer, self).monosaccharide_to_iupac(residue)
        deriv = has_derivatization(residue)
        if deriv:
            string = "%s^%s" % (string, self.substituent_resolver.serialize_substituent(deriv))
        return string


class SimpleMonosaccharideSerializer(DerivatizationAwareMonosaccharideSerializer):
    def monosaccharide_to_iupac(self, residue):
        """
        Encode a subset of traits of a :class:`Monosaccharide`-like object
        using a limited subset of the IUPAC three letter code.

        Parameters
        ----------
        residue: :class:`~Monosaccharide`
            The object to be encoded

        Returns
        -------
        str
        """
        template = "{modification}{base_type}{substituent}"
        modification = ""
        base_type = self.resolve_special_base_type(residue)
        if base_type is None:
            if len(residue.stem) == 1 and residue.stem[0] is not Stem.Unknown:
                base_type = residue.stem[0].name.title()
            else:
                base_type = residue.superclass.name.title()
        modification = self.modification_extractor(residue.modifications, base_type)
        substituent = self.substituent_resolver(residue)
        string = template.format(
            modification=modification,
            base_type=base_type,
            substituent=substituent
        )

        deriv = has_derivatization(residue)
        if deriv:
            string = "%s^%s" % (string, self.substituent_resolver.serialize_substituent(deriv))
        return string


monosaccharide_to_iupac = MonosaccharideSerializer()
resolve_special_base_type = monosaccharide_to_iupac.resolve_special_base_type


class LinkageSerializer(object):
    def __init__(self, open_edge='-(', close_edge=')-', open_branch='[', close_branch=']'):
        self.open_edge = open_edge
        self.close_edge = close_edge
        self.open_branch = open_branch
        self.close_branch = close_branch

    def format_linkage(self, linkage):
        text = "{oe}{attach}-{linkage_pos}{ce}".format(
            linkage_pos=linkage.parent_position if linkage.parent_position != UnknownPosition else "?",
            attach=linkage.child_position if linkage.child_position != UnknownPosition else "?",
            oe=self.open_edge, ce=self.close_edge)
        return text

    def format_branch(self, branch):
        branch = '{ob}{branch}{cb}'.format(
            branch=''.join(branch),
            ob=self.open_branch,
            cb=self.close_branch
        )
        return branch


class SimpleLinkageSerializer(LinkageSerializer):
    def __init__(self, open_edge="(", close_edge=")", open_branch="[", close_branch="]"):
        super(SimpleLinkageSerializer, self).__init__(open_edge, close_edge, open_branch, close_branch)

    def format_linkage(self, linkage):
        template = "{oe}{anomer}{attach}-{linkage_pos}{ce}"
        text = template.format(
            oe=self.open_edge, anomer=anomer_map_to.get(linkage.child.anomer, "?"),
            linkage_pos=linkage.parent_position if linkage.parent_position != UnknownPosition else "?",
            attach=linkage.child_position if linkage.child_position != UnknownPosition else "?",
            ce=self.close_edge)
        return text


class GlycanSerializer(object):
    """Converts a :class:`~.Glycan` structure to IUPAC format.

    Also works on individual :class:`~.Monosaccharide` objects, but
    will traverse any links they have to other nodes.

    Attributes
    ----------
    linkage_serializer : :class:`LinkageSerializer`
        An object that converts a :class:`~.Link` object into text
    monosaccharide_serializer : :class:`MonosaccharideSerializer`
        An object that converts a :class:`~.Monosaccharide` object into text
    """

    def __init__(self, monosaccharide_serializer=None, linkage_serializer=None):
        if monosaccharide_serializer is None:
            monosaccharide_serializer = MonosaccharideSerializer()
        if linkage_serializer is None:
            linkage_serializer = LinkageSerializer()
        self.monosaccharide_serializer = monosaccharide_serializer
        self.linkage_serializer = linkage_serializer

    def branch_to_iupac(self, structure=None, attach=None, is_branch=False):
        '''Translate a |Glycan| structure's branch into IUPAC Three Letter Code.
        Recursively operates on branches.

        Parameters
        ----------
        structure: :class:`~.Glycan` or :class:`~.Monosaccharide`
            The glycan to be translated. Translation starts from `glycan.root` if `structure`
            is a |Glycan|. May also be a :class:`~.Monosaccharide` which is the root of a
            branch of the overall structure.
        attach: int
            The point from the structure tree is attached to its parent. Used for recursively
            handling branches. Defaults to |None|.
        is_branch: :class:`bool`
            Whether this structure contains the root of the overall structure or a branch

        Returns
        -------
        :class:`collections.deque`
        '''
        base = structure.root if isinstance(structure, Glycan) else structure
        stack = [(attach, base)]
        outstack = deque()
        while(len(stack) > 0):
            outedge, node = stack.pop()
            link = ""
            if outedge is not None:
                link = self.linkage_serializer.format_linkage(outedge)
            # Branch linkage does not start with leading dash
            if is_branch and link[-1] == '-':
                link = link[:-1]
            outstack.appendleft('{node}{link}'.format(node=self.monosaccharide_serializer(node), link=link))
            # Reset for next pass through the loop
            is_branch = False
            children = list((p, link) for p, link in node.links.items() if link.is_parent(node))
            if len(children) > 1:
                for pos, link in children[:-1]:
                    branch = self.linkage_serializer.format_branch(
                        self.branch_to_iupac(link.child, link, is_branch=True))
                    outstack.appendleft(branch)
                pos, link = children[-1]
                stack.append((link, link.child))
            elif len(children) == 1:
                pos, link = children[0]
                stack.append((link, link.child))
        return outstack

    def glycan_to_iupac(self, structure, **kwargs):
        '''Translate a |Glycan| structure's branch into IUPAC Three Letter Code.structure

        Calls :meth:`branch_to_iupac`, a recursive function.

        Parameters
        ----------
        structure: Glycan or Monosaccharide
            The glycan to be translated. Translation starts from `glycan.root` if `structure`
            is a |Glycan|.

        Returns
        -------
        :class:`str`
        '''
        return ''.join(self.branch_to_iupac(structure))

    def __call__(self, structure):
        """An alias for :meth:`glycan_to_iupac`
        """
        return self.glycan_to_iupac(structure)


glycan_to_iupac = GlycanSerializer()

glycan_to_iupac_simple = GlycanSerializer(SimpleMonosaccharideSerializer(), SimpleLinkageSerializer())


def to_iupac(structure, dialect=None):
    """Translate `structure` into its textual representation using IUPAC Three Letter Code

    Parameters
    ----------
    structure : |Glycan| or |Monosaccharide|
        The structure to be translated
    dialect: :class:`str`
        One of "extended" or "simple", controlling whether the long-form linkage
        and monosaccharide notation is used, or the more compact simplified form
        is used. Defaults to "extended".
    Returns
    -------
    |str|

    See Also
    --------
    :class:`GlycanSerializer`
    """
    if dialect is None:
        dialect = 'extended'
    if isinstance(structure, Monosaccharide):
        return monosaccharide_to_iupac(structure)
    else:
        if dialect == 'simple':
            return glycan_to_iupac_simple(structure)
        return glycan_to_iupac(structure)


def aminate_substituent(substituent):
    if substituent.name.startswith("n_"):
        # already aminated
        return substituent
    aminated = Substituent("n_" + substituent.name)
    if aminated.composition == {}:
        raise ValueError("Could not aminate substituent")
    return aminated


class SubstituentDeserializer(object):
    def __init__(self, substituents_map=None, error_on_missing=True):
        if substituents_map is None:
            substituents_map = substituents_map_from
        self.error_on_missing = error_on_missing
        self.substituents_map = substituents_map

    def substituent_from_iupac(self, substituents):
        parts = re.split(r"\(|\)", substituents)
        for part in parts:
            if part == "":
                continue
            # split_part = re.split(r"(\d+)?", part)
            split_part = re.split(r"(\d+)", part)

            if len(split_part) == 3:
                _, position, name = split_part
            else:
                position = UnknownPosition
                name = split_part[0]
            try:
                name = self.substituents_map[name]
            except KeyError:
                # Acidic special case:
                # Often, acidic monosaccharides are written with a trailing A like a substituent while
                # GlycoCT treats acidic groups as modifications. If an A appears in the substituent suffix
                # it will fail to be cast as a substituent, but pass it along raw and it will be handled
                # downstream by :func:`monosaccharide_from_iupac`.
                if name == "A":
                    pass
                elif name.startswith("O"):
                    if name[1:] in self.substituents_map:
                        # Some dialects prefix non-amine substituents with O to differentiate them
                        name = self.substituents_map[name[1:]]
                else:  # pragma: no cover
                    if not self.error_on_missing:
                        warnings.warn("No translation rule found to convert %s into a Substituent" % name)
                        continue
                    else:
                        raise IUPACError("No translation rule found to convert %s into a Substituent" % name)
            yield int(position), name

    def symbol_to_name(self, symbol):
        name = self.substituents_map[symbol]
        return name

    def __call__(self, substituents):
        return self.substituent_from_iupac(substituents)


substituent_from_iupac = SubstituentDeserializer()


LinkageSpecification = namedtuple("LinkageSpecification", ("child_position", "parent_position", "has_ambiguity"))


class LinkageDeserializer(object):
    pattern = re.compile(r"\((?P<child_linkage>[0-9?/]+)->?(?P<parent_linkage>[0-9?/]+)?\)?")

    def parse(self, linkage_string):
        if linkage_string is None:
            return None
        match = self.pattern.search(linkage_string)
        if match is None:
            raise ValueError(linkage_string)
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


class SimpleLinkageDeserializer(LinkageDeserializer):
    pattern = re.compile(r"""\((?P<anomer>[abo?])
                             (?P<child_linkage>[0-9?/]+)->?
                             (?P<parent_linkage>[0-9?/]+)?\)?""",
                         re.VERBOSE)


parse_linkage_structure = LinkageDeserializer()


class MonosaccharideDeserializer(object):
    _pattern = r'''(?:(?P<anomer>[abo?]|alpha|beta|\u03B1|\u03B2)-)?
                   (?P<configuration>[LD?])-
                   (?P<modification>[a-z0-9_\-,]*?)
                   (?P<base_type>(?:[A-Z][a-z]{2}?|(?:[a-z]{3}[A-Z][a-z]{2})))
                   (?P<ring_type>[xpfo?])?
                   (?P<substituent>[^-]*?)
                   (?P<linkage>-\([0-9?/]+->?[0-9?/]+\)-?)?
                   $'''
    try:
        # convert to unicode for Py2
        _pattern = _pattern.decode("raw_unicode_escape")
    except AttributeError:
        pass
    pattern = re.compile(_pattern, re.VERBOSE | re.UNICODE)

    linkage_parser = parse_linkage_structure

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

    def build_residue(self, match_dict):
        try:
            anomer = anomer_map_from[match_dict['anomer']]
        except KeyError:
            anomer = anomer_map_from['?']
        base_type = match_dict["base_type"]
        configuration = match_dict["configuration"].lower()
        ring_type = match_dict['ring_type']

        modification = (match_dict['modification'] or '').rstrip("-")

        linkage = match_dict.get("linkage")
        original_base_type = base_type
        # alternate carbon backbone size encoded as stem{3}Superclass{3}
        # instead of Stem{3}
        if len(base_type) == 6:
            superclass_type = base_type[3:].lower()
            base_type = base_type[:3].title()
        else:
            superclass_type = None
        try:
            residue = named_structures.monosaccharides[base_type]
        except KeyError:
            raise IUPACError("Unknown Residue Base-type %r" % (original_base_type,))
        base_is_modified = len(residue.substituent_links) + len(residue.modifications) > 0
        if superclass_type is not None:
            residue.superclass = superclass_type
            residue.ring_start = UnknownPosition
            residue.ring_end = UnknownPosition

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
            if position == UnknownPosition and base_is_modified:
                # Guess at what the user might mean using base_type
                if base_type == "Neu" and substituent in ["acetyl", "glycolyl"] and i == 1:
                    position = 5
                # else:
                #     raise ValueError(
                #         "Cannot have ambiguous location of substituents on a base type which"
                #         " has default modifications or substituents. {} {}".format(
                #             residue, (position, substituent)))
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

    def parse_linkage_structure(self, linkage):
        return self.linkage_parser(linkage)

    def monosaccharide_from_iupac(self, monosaccharide_str, parent=None):
        match_dict = self.extract_pattern(monosaccharide_str)
        residue, linkage = self.build_residue(match_dict)
        linkage = self.parse_linkage_structure(linkage)

        self.add_monosaccharide_bond(residue, parent, linkage)
        return residue, linkage

    def add_monosaccharide_bond(self, residue, parent, linkage):
        if parent is not None and linkage != ():
            if linkage.has_ambiguity:
                bond = AmbiguousLink(
                    parent, residue, parent_position=linkage.parent_position,
                    child_position=linkage.child_position,
                    parent_loss=Composition("H"), child_loss=Composition("OH"))
                bond.find_open_position()
            else:
                parent.add_monosaccharide(residue, position=linkage[1], child_position=linkage[0])

    def __call__(self, monosaccharide_str, parent=None):
        return self.monosaccharide_from_iupac(monosaccharide_str, parent=parent)

    def finalize(self, glycan):
        pass


class DerivatizationAwareMonosaccharideDeserializer(MonosaccharideDeserializer):
    _pattern = r'''(?:(?P<anomer>[abo?]|alpha|beta|\u03B1|\u03B2)-)?
                   (?P<configuration>[LD?])-
                   (?P<modification>[a-z0-9_\-,]*)
                   (?P<base_type>[^-]{3}?)
                   (?P<ring_type>[xpfo?])?
                   (?P<substituent>[^-]*?)
                   (?P<derivatization>\^[^\s-]*?)?
                   (?P<linkage>-\([0-9?/]+->?[0-9?/]+\)-?)?$'''
    try:
        # convert to unicode for Py2
        _pattern = _pattern.decode("raw_unicode_escape")
    except AttributeError:
        pass
    pattern = re.compile(_pattern, re.VERBOSE | re.UNICODE)

    def add_monosaccharide_bond(self, residue, parent, linkage):
        if parent is not None and linkage != ():
            try:
                if linkage.has_ambiguity:
                    bond = AmbiguousLink(
                        parent, residue, parent_position=linkage.parent_position,
                        child_position=linkage.child_position,
                        parent_loss=Composition("H"), child_loss=Composition("OH"))
                    bond.find_open_position()
                else:
                    parent.add_monosaccharide(residue, position=linkage[1], child_position=linkage[0])
            except ValueError:
                parent_substituent_links_at_site = parent.substituent_links[linkage[1]]
                if (parent_substituent_links_at_site and parent_substituent_links_at_site[0].child._derivatize):
                    parent.drop_substituent(linkage[1], parent_substituent_links_at_site[0].child)
                residue_substituent_links_at_site = residue.substituent_links[linkage[0]]
                if residue_substituent_links_at_site and residue_substituent_links_at_site[0].child._derivatize:
                    residue.drop_substituent(linkage[0], residue_substituent_links_at_site[0].child)

                if linkage.has_ambiguity:
                    bond = AmbiguousLink(
                        parent, residue, parent_position=linkage.parent_position,
                        child_position=linkage.child_position,
                        parent_loss=Composition("H"), child_loss=Composition("OH"))
                    bond.find_open_position()
                else:
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
        linkage = self.parse_linkage_structure(linkage)

        self.add_monosaccharide_bond(residue, parent, linkage)

        deriv = match_dict.get("derivatization", '')
        if deriv is not None and deriv != "":
            self.apply_derivatization(residue, deriv)

        return residue, linkage

    def finalize(self, glycan):
        for node in glycan:
            neg_capacity = -node._remaining_capacity()
            if neg_capacity > 0:
                unknowns = node.substituent_links[UnknownPosition]
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

    def strip_derivatization(self, residue_str, **kwargs):
        base = residue_str.rsplit("^")[0]
        return self(base, **kwargs)


class SimpleMonosaccharideDeserializer(DerivatizationAwareMonosaccharideDeserializer):
    pattern = re.compile(
        r'''(?P<modification>[a-z0-9_\-,]*)
            (?P<base_type>(?:[A-Z][a-z]{2}?|(?:[a-z]{3}[A-Z][a-z]{2})))
            (?P<ring_type>[pfox])?
            (?P<substituent>[^-]*?)
            (?P<derivatization>\^[^\s-]*?)?
            (?P<linkage>-?\((?P<anomer>[ab?o]?)[0-9?/]+->?[0-9?/]+\)-?)?$''', re.VERBOSE)

    linkage_parser = SimpleLinkageDeserializer()

    def parse_linkage_structure(self, linkage):
        return self.linkage_parser(linkage)

    def build_residue(self, match_dict):
        base_type = match_dict["base_type"]
        modification = (match_dict['modification'] or '').rstrip("-")

        try:
            residue = named_structures.monosaccharides[base_type]
        except KeyError:
            raise IUPACError("Unknown Residue Base-type %r" % (base_type,))
        base_is_modified = len(residue.substituent_links) + len(residue.modifications) > 0
        ring_type = match_dict.get('ring_type')
        if ring_type is not None:
            self.ring_bounds(residue, ring_type)
        self.set_modifications(residue, modification)
        self.set_substituents(residue, match_dict['substituent'], base_is_modified, base_type)
        linkage = match_dict.get("linkage")
        return residue, linkage


monosaccharide_from_iupac = MonosaccharideDeserializer()


class GlycanDeserializer(object):
    def __init__(self, monosaccharide_deserializer=None, set_default_positions=True):
        if monosaccharide_deserializer is None:
            monosaccharide_deserializer = MonosaccharideDeserializer()
        self.monosaccharide_deserializer = monosaccharide_deserializer
        self.set_default_positions = set_default_positions

    new_branch_open = re.compile(r"(\]-?)$")

    def add_monosaccharide(self, parent_node, child_node, linkage):
        # parent_node.add_monosaccharide(
        #     child_node, position=parent_position, child_position=child_position)
        self.monosaccharide_deserializer.add_monosaccharide_bond(
            child_node, parent_node, linkage)

    def glycan_from_iupac(self, text, structure_class=Glycan, **kwargs):
        last_outedge = None
        root = None
        last_residue = None
        branch_stack = []

        # Remove the base
        text = re.sub(r"\((\d*|\?)->?$", "", text)

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
            elif text[-1] == '[':
                try:
                    branch_parent, old_root, old_last_outedge = branch_stack.pop()
                    # child_position, parent_position = last_outedge
                    self.add_monosaccharide(branch_parent, root, last_outedge)
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
        res = structure_class(root=root)
        self.monosaccharide_deserializer.finalize(res)
        res.reindex()
        if self.set_default_positions:
            self.set_default_positions_for_common_cases(res)
        return res

    def __call__(self, text, **kwargs):
        return self.glycan_from_iupac(text, **kwargs)

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


def set_default_positions(glycan):  # pragma: nocover
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


glycan_from_iupac = GlycanDeserializer()

glycan_from_iupac_simple = GlycanDeserializer(SimpleMonosaccharideDeserializer())


def from_iupac(text, structure_class=Glycan, resolve_default_positions=True, dialect=None, **kwargs):
    """Parse the given text into an instance of |Glycan|. If there is only a single monosaccharide
    in the output, just the Monosaccharide instance is returned.

    Parameters
    ----------
    text : |str|
        The text to parser
    resolve_default_positions: :class:`bool`
        Whether to assume default positions for common monosaccharide modifiers
        that are omitted for brevity, such as the postion of n-acetyl on HexNAc.
    dialect: :class:`str`
        One of "extended" or "simple", controlling whether the long-form linkage
        and monosaccharide notation is used, or the more compact simplified form
        is used. Defaults to "extended".
    **kwargs:
        Forwarded to :func:`glycan_from_iupac`

    Returns
    -------
    |Glycan| or |Monosaccharide|
        If the resulting structure is just a single monosaccharide, the returned value is a Monosaccharide.
    """
    if dialect is None:
        dialect = 'extended'
    if dialect != 'simple':
        res = glycan_from_iupac(
            text, structure_class=structure_class,
            set_default_positions=resolve_default_positions,
            **kwargs)
    else:
        res = glycan_from_iupac_simple(
            text, structure_class=structure_class,
            set_default_positions=resolve_default_positions,
            **kwargs)
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


load = IUPACParser.load


Monosaccharide.register_serializer("iupac", dumps)
Glycan.register_serializer("iupac", dumps)

_dumps_simple = partial(dumps, dialect='simple')
_dumps_extended = partial(dumps, dialect='extended')

Monosaccharide.register_serializer("iupac_simple", _dumps_simple)
Glycan.register_serializer("iupac_simple", _dumps_simple)

Monosaccharide.register_serializer("iupac_extended", _dumps_extended)
Glycan.register_serializer("iupac_extended", _dumps_extended)
