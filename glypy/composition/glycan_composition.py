import re
from glypy import Composition, Monosaccharide, Glycan, Substituent
from glypy.utils import tree
from glypy.utils.multimap import OrderedMultiMap

from glypy.structure.base import SaccharideCollection
from glypy.io.iupac import (
    parse_modifications, named_structures, Modification,
    substituent_from_iupac, Stem, extract_modifications,
    resolve_substituent, SuperClass, Configuration, aminate_substituent,
    resolve_special_base_type as _resolve_special_base_type)


from .composition_transform import (
    derivatize, has_derivatization, strip_derivatization,
    _derivatize_reducing_end, _strip_derivatization_reducing_end)

monosaccharide_parser_lite = re.compile(r'''(?P<modification>[a-z0-9_\-,]*)
                                       (?P<base_type>[A-Z][a-z]+)
                                       (?P<substituent>[^-]*?)$''', re.VERBOSE)


def from_iupac_lite(monosaccharide_str):
    """
    Parse a string in a limited subset of IUPAC three letter code into
    an instance of :class:`MonosaccharideResidue`.

    Parameters
    ----------
    monosaccharide_str: str
        The string to be parsed

    Returns
    -------
    MonosaccharideResidue
    """
    match = monosaccharide_parser_lite.search(monosaccharide_str)
    if match is None:
        raise ValueError("Cannot find monosaccharide pattern in {}".format(monosaccharide_str))
    match_dict = match.groupdict()
    base_type = match_dict["base_type"]

    modification = match_dict['modification']

    residue = named_structures.monosaccharides[base_type]
    base_is_modified = len(residue.substituent_links) + len(residue.modifications) > 0

    residue.ring_end = residue.ring_start = None

    for pos, mod in parse_modifications(modification):
        residue.add_modification(mod, pos)
    i = 0
    strict = False
    for position, substituent in substituent_from_iupac(match_dict["substituent"]):
        i += 1
        if position == -1 and base_is_modified:
            # Guess at what the user might mean using base_type
            if base_type == "Neu" and substituent in ["acetyl", "glycolyl"] and i == 1:
                position = 5
            elif strict:
                raise ValueError(
                    "Cannot have ambiguous location of substituents on a base type which"
                    " has default modifications or substituents. {} {}".format(
                        residue, (position, substituent)))

        substituent = Substituent(substituent)
        try:
            residue.add_substituent(
                substituent, position,
                parent_loss=substituent.attachment_composition_loss(),
                child_loss='H')
        except ValueError:
            # Highly modified large bases have a degenerate encoding, where additional qualifications following
            # base name *replace* an existing substituent. This behavior may not be expected in other more
            # common cases.
            occupancy = 0
            if base_type in {"Neu", "Kdo"}:
                try:
                    unplaced = residue.substituent_links[position][0].child
                    residue.drop_substituent(position)
                    if unplaced.name == "amino":
                        try:
                            substituent = aminate_substituent(substituent)
                        except ValueError:
                            pass
                except:
                    # The site contains a modification which can be present alongside the substituent
                    occupancy = 1
                try:
                    residue.add_substituent(
                        substituent, position, occupancy,
                        parent_loss=substituent.attachment_composition_loss(),
                        child_loss='H')
                except:
                    residue.add_substituent(
                        substituent, -1, occupancy,
                        parent_loss=substituent.attachment_composition_loss(),
                        child_loss='H')
            else:
                residue.add_substituent(
                    substituent, -1, occupancy,
                    parent_loss=substituent.attachment_composition_loss(),
                    child_loss='H')

    return MonosaccharideResidue.from_monosaccharide(residue)


def to_iupac_lite(residue):
    """
    Encode a subset of traits of a :class:`Monosaccharide`-like object
    using a limited subset of the IUPAC three letter code. The information
    present is sufficient to reconstruct a :class:`MonosaccharideResidue` instance
    reflecting the base type and its native substituents and modificats.

    .. note::
        This function is not suitable for use on whole |Glycan| objects. Instead,
        see :meth:`GlycanComposition.from_glycan` and :meth:`GlycanComposition.serialize`

    Parameters
    ----------
    residue: Monosaccharide
        The object to be encoded

    Returns
    -------
    str

    See Also
    --------
    :func:`from_iupac_lite`
    """
    template = "{modification}{base_type}{substituent}"
    modification = ""
    base_type = _resolve_special_base_type(residue)
    if base_type is None:
        if residue.stem[0] is not Stem.Unknown:
            base_type = residue.stem[0].name.title()
        else:
            base_type = residue.superclass.name.title()

    # Omit unknown coordinates on modifications and substituents
    modification = extract_modifications(residue.modifications, base_type).replace("-1-", "")
    substituent = resolve_substituent(residue).replace("-1", "")
    return template.format(
        modification=modification,
        base_type=base_type,
        substituent=substituent
        )


def drop_stem(residue):
    if _resolve_special_base_type(residue) is None:
        residue.stem = (None,)


def drop_positions(residue):
    if _resolve_special_base_type(residue) is None:
        modifications = OrderedMultiMap()
        for k, v in residue.modifications.items():
            modifications[-1] = v
        residue.modifications = modifications

        for p, link in list(residue.substituent_links.items()):
            link.break_link(refund=True)
            link.parent_position = -1
            link.apply()


def drop_configuration(residue):
    if _resolve_special_base_type(residue) is None:
        residue.configuration = (None,)


water_composition = {"O": 1, "H": 2}


class MonosaccharideResidue(Monosaccharide):

    @classmethod
    def from_monosaccharide(cls, monosaccharide, configuration=False, stem=True, ring=False):
        """Construct an instance of :class:`MonosaccharideResidue` from an instance
        of |Monosaccharide|. This function attempts to preserve derivatization if possible.

        This function will create a *deep copy* of `monosaccharide`.

        Parameters
        ----------
        monosaccharide : Monosaccharide
            The monosaccharide to be converted
        configuration : bool, optional
            Whether or not to preserve |Configuration|. Defaults to |False|
        stem : bool, optional
            Whether or not to preserve |Stem|. Defaults to |True|
        ring : bool, optional
            Whether or not to preserve |RingType|. Defaults to |False|

        Returns
        -------
        MonosaccharideResidue
        """
        residue = monosaccharide.clone(monosaccharide_type=cls)
        premass = residue.mass()

        deriv = has_derivatization(monosaccharide)
        strip_derivatization(residue)
        if _resolve_special_base_type(monosaccharide) is None:
            if not configuration:
                residue.configuration = (None,)
            if not stem:
                residue.stem = (None,)
        if not ring:
            residue.ring_start = residue.ring_end = None
        if deriv:
            derivatize(residue, deriv)
        if residue.mass() != premass and not deriv:
            residue.composition += water_composition
        return residue

    def __init__(self, *args, **kwargs):
        super(MonosaccharideResidue, self).__init__(*args, **kwargs)
        self.composition -= water_composition
        self.anomer = None

    def clone(self, *args, **kwargs):
        kwargs.setdefault("monosaccharide_type", MonosaccharideResidue)
        residue = super(MonosaccharideResidue, self).clone(*args, **kwargs)
        return residue

    def __repr__(self):  # pragma: no cover
        return "MonosaccharideResidue(%s)" % self.name()

    def __str__(self):  # pragma: no cover
        return to_iupac_lite(self)

    def __hash__(self):  # pragma: no cover
        """Obtain a hash value from `self` based on :meth:`MonosaccharideResidue.name`.

        Returns
        -------
        int
        """
        return hash(self.name())

    def open_attachment_sites(self, max_occupancy=0):
        sites, unknowns = super(
            MonosaccharideResidue, self).open_attachment_sites(max_occupancy)
        return sites[:-3], unknowns

    def __eq__(self, other):
        '''
        Test for equality between :class:`MonosaccharideResidue` instances by comparing
        the result of :meth:`MonosaccharideResidue.name` calls between `self` and `other`.

        :meth:`MonosaccharideResidue.name` is an alias of :func:`to_iupac_lite` called on `self`
        '''
        if (other is None):
            return False
        if not isinstance(other, MonosaccharideResidue):
            return False
        return self.name() == other.name()

    from_iupac_lite = staticmethod(from_iupac_lite)
    name = to_iupac_lite

    drop_stem = drop_stem
    drop_positions = drop_positions
    drop_configuration = drop_configuration

water_mass = Composition("H2O").mass


class GlycanComposition(dict, SaccharideCollection):

    @classmethod
    def from_glycan(cls, glycan):
        inst = cls()
        glycan = tree(glycan)
        inst.extend(glycan)
        inst.reducing_end = glycan.reducing_end
        deriv = has_derivatization(glycan.root)
        if deriv:
            inst._composition_offset += (
                deriv.total_composition() - deriv.attachment_composition_loss()) * 2
        return inst

    def __init__(self, *args, **kwargs):
        self._reducing_end = kwargs.pop("reducing_end", None)
        dict.__init__(self)
        self._mass = None
        self._composition_offset = Composition("H2O")
        self.update(*args, **kwargs)

    def __setitem__(self, key, value):
        if isinstance(key, basestring):
            key = from_iupac_lite(key)
        if key.reducing_end is not None:
            self.reducing_end = key.reducing_end
            key.reducing_end = None
        dict.__setitem__(self, key, value)
        self._mass = None

    def __getitem__(self, key):
        if isinstance(key, basestring):
            key = from_iupac_lite(key)
        return dict.__getitem__(self, key)

    def mass(self, average=False, charge=0, mass_data=None):
        if self._mass is not None:
            return self._mass
        mass = self._composition_offset.mass
        for residue_type, count in list(self.items()):
            mass += residue_type.mass(average=average, charge=charge, mass_data=mass_data) * count
        if self._reducing_end is not None:
            mass += self._reducing_end.mass(average=average, charge=charge, mass_data=mass_data)
        self._mass = mass
        return mass

    def update(self, *args, **kwargs):
        if len(args) == 1:
            if isinstance(args[0], dict):
                args = list(args)
                for name, count in args[0].items():
                    self[name] = count
            else:
                for name, count in args:
                    self[name] = count
        for name, count in kwargs.items():
            self[name] = count
        self._mass = None

    def extend(self, *args):
        if not isinstance(args[0], MonosaccharideResidue):
            if isinstance(args[0], (Monosaccharide)):
                args = map(MonosaccharideResidue.from_monosaccharide, args)
            elif isinstance(args[0], Glycan):
                args = map(
                    MonosaccharideResidue.from_monosaccharide,
                    [node for node in args[0] if node.node_type is MonosaccharideResidue.node_type])
            else:
                raise TypeError(
                    "Can't convert {} to MonosaccharideResidue".format(
                        type(args[0])))
        for residue in args:
            self[residue] += 1

    def __iadd__(self, other):
        for elem, cnt in (other.items()):
            self[elem] += cnt
        return self

    def __add__(self, other):
        result = self.clone()
        for elem, cnt in other.items():
            result[elem] += cnt
        return result

    def __radd__(self, other):
        return self + other

    def __isub__(self, other):
        for elem, cnt in other.items():
            self[elem] -= cnt
        return self

    def __sub__(self, other):
        result = self.clone()
        for elem, cnt in other.items():
            result[elem] -= cnt
        return result

    def __rsub__(self, other):
        return (self - other) * (-1)

    def __mul__(self, other):
        if not isinstance(other, int):
            raise TypeError(
                'Cannot multiply Composition by non-integer',
                other)
        prod = {}
        for k, v in self.items():
            prod[k] = v * other

        return GlycanComposition(prod)

    def __rmul__(self, other):
        return self * other

    def __eq__(self, other):
        if not isinstance(other, dict):
            return False
        self_items = set([i for i in self.items() if i[1]])
        other_items = set([i for i in other.items() if i[1]])
        return self_items == other_items

    def __neg__(self):
        return -1 * self

    def __missing__(self, key):
        return 0

    def drop_stems(self):
        for t in self:
            drop_stem(t)
        return self

    def drop_positions(self):
        for t in self:
            drop_positions(t)
        return self

    def drop_configurations(self):
        for t in self:
            drop_configuration(t)

    def total_composition(self):
        comp = self._composition_offset.clone()
        for residue, count in self.items():
            comp += residue.total_composition() * count
        if self._reducing_end is not None:
            comp += self._reducing_end.total_composition()
        return comp

    def collapse(self):
        '''
        Merge redundant keys.

        After performing a structure-detail removing operation like
        :meth:`drop_positions`, :meth:`drop_configurations`, or :meth:`drop_stems`,
        monosaccharide keys may be redundant.

        `collapse` will merge keys which refer to the same type of molecule.
        '''
        items = list(self.items())
        self.clear()
        for k, v in items:
            self[k] += v

    @property
    def reducing_end(self):
        return self._reducing_end

    @reducing_end.setter
    def reducing_end(self, value):
        self._mass = None
        self._reducing_end = value

    def set_reducing_end(self, value):
        self._mass = None
        self._reducing_end = value

    @property
    def composition_offset(self):
        return self._composition_offset

    @composition_offset.setter
    def composition_offset(self, value):
        self._mass = None
        self._composition_offset = value

    def clone(self):
        return self.__class__(self)

    def serialize(self):
        return "{%s}" % '; '.join("{}:{}".format(str(k), v) for k, v in sorted(
            self.items(), key=lambda x: x[0].mass()))

    __str__ = serialize

    @classmethod
    def parse(cls, string):
        inst = cls()
        tokens = string[1:-1].split('; ')
        for token in tokens:
            residue, count = token.split(":")
            inst[from_iupac_lite(residue)] = int(count)
        return inst

    def _derivatized(self, substituent, id_base):
        self._composition_offset += (
            substituent.total_composition() -
            substituent.attachment_composition_loss() * 2) * 2
        if self._reducing_end is not None:
            _derivatize_reducing_end(self._reducing_end, substituent, id_base)
        self._mass = None

    def _strip_derivatization(self):
        self._composition_offset = Composition("H2O")
        if self._reducing_end is not None:
            _strip_derivatization_reducing_end(self._reducing_end)
        self._mass = None

from_glycan = GlycanComposition.from_glycan
parse = GlycanComposition.parse
