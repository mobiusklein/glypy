'''
:class:`GlycanComposition`, :class:`MonosaccharideResidue`, and :class:`SubstituentResidue` are
useful for working with bag-of-residues where topology and connections are not relevant, but
the aggregate composition is known. These types work with a subset of the IUPAC three letter code
for specifying compositions.


>>> g = GlycanComposition(Hex=3, HexNAc=2)
>>> g["Hex"]
3
>>> r = MonosaccharideResidue.from_iupac_lite("Hex")
>>> r
MonosaccharideResidue(Hex)
>>> g[r]
3
>>> import glypy
>>> abs(g.mass() - glypy.motifs["N-Glycan core basic 1"].mass()) < 1e-5
True
>>> g2 = GlycanComposition(Hex=5)
>>> g["@n-acetyl"] = -2 # Remove two n-acetyl groups from the composition
>>> abs(g.mass() - g2.mass()) < 1e-5
True

'''
import re
from glypy import Composition, Monosaccharide, Glycan, Substituent
from glypy.utils import tree
from glypy.utils.multimap import OrderedMultiMap

from glypy.structure.base import SaccharideCollection, MoleculeBase

from glypy.io.iupac import (
    parse_modifications, named_structures, Modification,
    substituent_from_iupac, Stem, extract_modifications,
    resolve_substituent, aminate_substituent, monosaccharide_reference as _monosaccharide_reference,
    resolve_special_base_type as _resolve_special_base_type)


from glypy.composition.composition_transform import (
    derivatize, has_derivatization, strip_derivatization,
    _derivatize_reducing_end, _strip_derivatization_reducing_end)

monosaccharide_parser_lite = re.compile(r'''(?P<modification>[a-z0-9_\-,]*)
                                       (?P<base_type>[A-Z][a-z]+)
                                       (?P<substituent>[^-]*?)$''', re.VERBOSE)


monosaccharide_residue_reference = {}


def from_iupac_lite(monosaccharide_str, residue_class=None):
    """
    Parse a string in a limited subset of IUPAC three letter code into
    an instance of :class:`MonosaccharideResidue` or :class:`SubstituentResidue`.

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
        try:
            result = SubstituentResidue.from_iupac_lite(monosaccharide_str)
            return result
        except:
            try:  # pragma: no cover
                result = MolecularComposition.from_iupac_lite(monosaccharide_str)
                return result
            except:
                raise ValueError("Cannot find pattern in {}".format(monosaccharide_str))
    if residue_class is None:
        residue_class = MonosaccharideResidue
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

        # Often, acidic monosaccharides will be suffixed "A" instead of prefixed "a".
        # Handle this here.
        if substituent == "A":
            residue.add_modification(Modification.a, position)
            continue

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
                        except ValueError:  # pragma: no cover
                            pass
                except:  # pragma: no cover
                    # The site contains a modification which can be present alongside the substituent
                    occupancy = 1
                try:
                    residue.add_substituent(
                        substituent, position, occupancy,
                        parent_loss=substituent.attachment_composition_loss(),
                        child_loss='H')
                except:  # pragma: no cover
                    residue.add_substituent(
                        substituent, -1, occupancy,
                        parent_loss=substituent.attachment_composition_loss(),
                        child_loss='H')
            else:  # pragma: no cover
                residue.add_substituent(
                    substituent, -1, occupancy,
                    parent_loss=substituent.attachment_composition_loss(),
                    child_loss='H')

    return residue_class.from_monosaccharide(residue)


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
    if isinstance(residue, (SubstituentResidue, MolecularComposition)):
        return residue.to_iupac_lite()

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
    substituent = resolve_substituent(residue, monosaccharide_residue_reference).replace("-1", "")
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
        return sites[:-2], unknowns

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
        return str(self) == str(other)

    from_iupac_lite = staticmethod(from_iupac_lite)
    name = to_iupac_lite

    drop_stem = drop_stem
    drop_positions = drop_positions
    drop_configuration = drop_configuration


monosaccharide_residue_reference.update({
    k: MonosaccharideResidue.from_monosaccharide(v) for k, v in _monosaccharide_reference.items()
    })


class FrozenMonosaccharideResidue(MonosaccharideResidue):
    '''
    A subclass of |MonosaccharideResidue| which caches the result of :func:`to_iupac_lite` and instances returned
    by :meth:`FrozenMonosaccharideResidue.clone` and :meth:`FrozenMonosaccharideResidue.from_iupac_lite`.
    Also treated as immutable after initialization through :meth:`FrozenMonosaccharideResidue.from_monosaccharide`.

    Note that directly calling :meth:`FrozenMonosaccharideResidue.from_monosaccharide` will not retrieve instances
    from the cache directly, and direct initialization using normal instance creation will neither touch the cache
    nor freeze the instance.

    This type is intended for use with :class:`FrozenGlycanComposition` to minimize the number of times
    :func:`from_iupac_lite` is called.
    '''

    _frozen = False
    _cache = {}

    @classmethod
    def from_monosaccharide(cls, monosaccharide, *args, **kwargs):
        if has_derivatization(monosaccharide):
            raise FrozenError("Cannot create Frozen derivatize type")
        inst = super(FrozenMonosaccharideResidue, cls).from_monosaccharide(monosaccharide, *args, **kwargs)
        if str(inst) not in inst._cache:
            inst._cache[str(inst)] = inst
            inst._frozen = True
        else:
            inst = inst._cache[str(inst)]
        return inst

    def __init__(self, *args, **kwargs):
        super(FrozenMonosaccharideResidue, self).__init__(*args, **kwargs)
        self._frozen = kwargs.get("_frozen", False)

    def __setattr__(self, key, value):
        if self._frozen:
            self._cache.pop(self._name, None)
            raise FrozenError("Cannot change a frozen object")
        else:
            object.__setattr__(self, key, value)

    def __repr__(self):  # pragma: no cover
        return "FrozenMonosaccharideResidue(%s)" % self.name()

    def __hash__(self):  # pragma: no cover
        """Obtain a hash value from `self` based on :meth:`MonosaccharideResidue.name`.

        Returns
        -------
        int
        """
        return hash(self._name)

    def _save_to_cache(self):
        self._cache[str(self)] = self

    def __str__(self):
        try:
            return self._name
        except:
            name = to_iupac_lite(self)
            self._name = name
            return name

    def clone(self, *args, **kwargs):
        if self._frozen:
            return self
        else:
            return super(FrozenMonosaccharideResidue, self).clone(*args, **kwargs)

    @classmethod
    def from_iupac_lite(cls, string):
        try:
            return cls._cache[string]
        except KeyError:
            return from_iupac_lite(string, cls)


class SubstituentResidue(Substituent):
    r'''
    Represent substituent molecules unassociated with a specific
    monosaccharide residue.


    .. note::
        :class:`SubstituentResidue`'s composition value includes the losses for forming a bond between
        a monosaccharide residue and the substituent.

    Attributes
    ----------
    name: str
        As in |Substituent|, but with :attr:`SubstituentResidue.sigil` prepended.
    composition: |Composition|
    links: |OrderedMultiMap|
    _order: |int|
    '''
    #: All substituent string identifiers are prefixed with this character
    #: for the :func:`from_iupac_lite` parser
    sigil = "@"

    def __init__(self, name, composition=None, id=None,
                 can_nh_derivatize=None, is_nh_derivatizable=None, derivatize=False,
                 attachment_composition=None):
        if name.startswith(SubstituentResidue.sigil):
            name = name[1:]
        elif name.startswith(MolecularComposition.sigil):
            raise TypeError("Invalid Sigil. SubstituentResidue instances must be given names with either"
                            " no sigil prefix or with '@'")
        super(SubstituentResidue, self).__init__(
            name=name, composition=composition, links=None, id=id,
            can_nh_derivatize=can_nh_derivatize, is_nh_derivatizable=is_nh_derivatizable,
            derivatize=derivatize, attachment_composition=attachment_composition)

        self.composition -= self.attachment_composition
        self.composition -= {"H": 1}
        self._name = SubstituentResidue.sigil + self._name

    def __hash__(self):
        return hash(self.name)

    def to_iupac_lite(self):
        return self.name

    __str__ = to_iupac_lite

    def __repr__(self):  # pragma: no cover
        return "SubstituentResidue(%s)" % self.name

    @classmethod
    def from_iupac_lite(cls, name):
        return cls(name)

    def __eq__(self, other):
        if (other is None):
            return False
        if not isinstance(other, SubstituentResidue):
            return False
        return self.name == other.name

    def __ne__(self, other):  # pragma: no cover
        return not self == other


class MolecularComposition(MoleculeBase):  # pragma: no cover
    sigil = "#"

    def __init__(self, name, composition):
        self.name = name
        self.composition = composition

    def mass(self, average=False, charge=0, mass_data=None):
        return self.composition.calc_mass(average=average, charge=charge, mass_data=mass_data)

    def __repr__(self):
        return "%s%s%s%s" % (
            self.sigil, self.name, self.sigil,
            ''.join("%s%d" % kv for kv in self.composition.items()))

    to_iupac_lite = __repr__

    def open_attachment_sites(self, *args, **kwargs):
        return 0

    def clone(self):
        return self.__class__(self.name, Composition(self.composition))

    def total_composition(self):
        return self.composition.clone()

    @classmethod
    def from_iupac_lite(cls, string):
        if not string.startswith(cls.sigil):
            raise TypeError("%s does not start with header %s" % (string, cls.sigil))
        _, header, composition = string.split("#")
        name = header
        return cls(name, Composition(composition))

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        try:
            return self.name == other or self.name == other.name
        except:
            return self.name == str(other)

    def __ne__(self, other):
        return not (self == other)


water_mass = Composition("H2O").mass


class GlycanComposition(dict, SaccharideCollection):
    """
    Describe a glycan  as a collection of :class:`MonosaccharideResidue` counts without
    explicit linkage information relating how each monosaccharide is connected to its neighbors.

    This class subclasses |dict|, and assumes that keys will either be :class:`MonosaccharideResidue`
    instances, :class:`SubstituentResidue` instances, or strings in `iupac_lite` format which will be parsed
    into one of these types. While other types may be used, this is not recommended. All standard |dict| methods
    are supported.

    |GlycanComposition| objects may be derivatized just as |Glycan| objects are, with
    :func:`glypy.composition.composition_transform.derivatize` and
    :func:`glypy.composition.composition_transform.strip_derivatization`.

    GlycanComposition objects also support composition arithmetic, and can be added or subtracted from each other
    or multiplied by an integer.

    As GlycanComposition is not a complete structure, they cannot be translated into text formats as
    full |Glycan| objects are. They may instead be converted to and from a short-form text notation using
    :meth:`GlycanComposition.serialize` and reconstructed from this format using :meth:`GlycanComposition.parse`.

    Attributes
    ----------
    reducing_end : |ReducingEnd|
        Describe the reducing end of the aggregate without binding it to a specific monosaccharide.
        This will contribute to composition and mass calculations.
    _composition_offset: |Composition|
        Account for the one water molecule's worth of composition left over from applying the "residue"
        transformation to each monosaccharide in the aggregate.
    """
    @classmethod
    def from_glycan(cls, glycan):
        """
        Convert a |Glycan| into a |GlycanComposition|.

        Parameters
        ----------
        glycan : Glycan
            The instance to be converted

        Returns
        -------
        GlycanComposition
        """
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
        self._charge = None
        self._composition_offset = Composition("H2O")
        self.update(*args, **kwargs)

    def __setitem__(self, key, value):
        """
        Set the quantity of `key` to `value`

        If `key` is a string, it will be passed through :func:`from_iupac_lite`

        If `key` has a reducing end value, that reducing end will be set on `self`

        Parameters
        ----------
        key : str, MonosaccharideResidue, SubstituentResidue, or MolecularComposition
            The entity to store
        value : int
            The value to store
        """
        if isinstance(key, basestring):
            key = from_iupac_lite(key)
        if key.node_type is Monosaccharide.node_type and key.reducing_end is not None:
            self.reducing_end = key.reducing_end
            key.reducing_end = None
        dict.__setitem__(self, key, int(value))
        self._mass = None

    def __getitem__(self, key):
        """
        Get the quantity of `key`

        If `key` is a string, it will be passed through :func:`from_iupac_lite`

        If `key` has a reducing end value, that reducing end will be set on `self`

        Parameters
        ----------
        key : str, MonosaccharideResidue, SubstituentResidue, or MolecularComposition
            The entity to store

        Returns
        -------
        int
        """
        if isinstance(key, basestring):
            key = from_iupac_lite(key)
        try:
            return dict.__getitem__(self, key)
        except KeyError:
            return 0

    def __delitem__(self, key):
        if isinstance(key, basestring):
            key = from_iupac_lite(key)
        dict.__delitem__(self, key)
        self._mass = None

    def mass(self, average=False, charge=0, mass_data=None):
        if self._mass is not None and charge == self._charge:
            return self._mass
        if charge == 0:
            mass = self._composition_offset.mass
            for residue_type, count in list(self.items()):
                mass += residue_type.mass(average=average, charge=0, mass_data=mass_data) * count
            if self._reducing_end is not None:
                mass += self._reducing_end.mass(average=average, charge=0, mass_data=mass_data)
            self._mass = mass
            self._charge = 0
        else:
            mass = self.total_composition().calc_mass(average=average, charge=charge, mass_data=mass_data)
            self._mass = mass
            self._charge = charge
        return mass

    def update(self, *args, **kwargs):
        if len(args) == 1:
            if isinstance(args[0], dict):
                args = list(args)
                for name, count in args[0].items():
                    if count != 0:
                        self[name] = count
            else:
                for name, count in args:
                    if count != 0:
                        self[name] = count
        for name, count in kwargs.items():
            if count != 0:
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

    def __contains__(self, key):
        if isinstance(key, basestring):
            key = from_iupac_lite(key)
        return dict.__contains__(self, key)

    def drop_stems(self):
        for t in self:
            drop_stem(t)
        self.collapse()

    def drop_positions(self):
        for t in self:
            drop_positions(t)
        self.collapse()

    def drop_configurations(self):
        for t in self:
            drop_configuration(t)
        self.collapse()

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
            self.items(), key=lambda x: x[0].mass()) if v > 0)

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
        n = 2
        for k, v in self.items():
            if k.node_type is Substituent.node_type:
                n -= v
        self._composition_offset += (
            substituent.total_composition() -
            substituent.attachment_composition_loss() * 2) * n
        if self._reducing_end is not None:
            _derivatize_reducing_end(self._reducing_end, substituent, id_base)
        self._mass = None

    def _strip_derivatization(self):
        self._composition_offset = Composition("H2O")
        if self._reducing_end is not None:
            _strip_derivatization_reducing_end(self._reducing_end)
        self._mass = None

    def _invalidate(self):
        self._mass = None
        self._charge = None

from_glycan = GlycanComposition.from_glycan
parse = GlycanComposition.parse


class FrozenGlycanComposition(GlycanComposition):
    '''
    A subclass of |GlycanComposition| which uses :class:`FrozenMonosaccharideResidue` instead
    of |MonosaccharideResidue| which reduces the number of times :func:`from_iupac_lite` is called.

    Only use this type if residue names are pre-validated, residue types will not be transformed,
    and when creating many, many instances. :func:`from_iupac_lite` invokes expensive introspection
    algorithms which can be costly when repeatedly manipulating the same residue types.
    '''

    _str = None

    def __setitem__(self, key, value):
        key = FrozenMonosaccharideResidue.from_iupac_lite(str(key))
        dict.__setitem__(self, key, value)
        self._mass = None

    def __getitem__(self, key):
        key = FrozenMonosaccharideResidue.from_iupac_lite(str(key))
        return dict.__getitem__(self, key)

    def __delitem__(self, key):
        key = FrozenMonosaccharideResidue.from_iupac_lite(str(key))
        dict.__delitem__(self, key)
        self._mass = None

    @classmethod
    def parse(cls, string):
        inst = cls()
        tokens = string[1:-1].split('; ')
        for token in tokens:
            residue, count = token.split(":")
            inst[FrozenMonosaccharideResidue.from_iupac_lite(residue)] = int(count)
        return inst

    def serialize(self):
        if self._mass is None or self._str is None:
            self._str = "{%s}" % '; '.join("{}:{}".format(str(k), v) for k, v in sorted(
                self.items(), key=lambda x: x[0].mass()) if v > 0)
        return self._str

    __str__ = serialize

    def __contains__(self, key):
        if isinstance(key, basestring):
            key = FrozenMonosaccharideResidue.from_iupac_lite(key)
        return dict.__contains__(self, key)

    def _derivatized(self, *args, **kwargs):
        raise FrozenError("Cannot derivatize a Frozen type")

    def _strip_derivatization(self, *args, **kwargs):
        raise FrozenError("Cannot derivatize a Frozen type")

    def extend(self, *args):
        if not isinstance(args[0], FrozenMonosaccharideResidue):
            if isinstance(args[0], (Monosaccharide)):
                args = map(FrozenMonosaccharideResidue.from_monosaccharide, args)
            elif isinstance(args[0], Glycan):
                args = map(
                    FrozenMonosaccharideResidue.from_monosaccharide,
                    [node for node in args[0] if node.node_type is FrozenMonosaccharideResidue.node_type])
            else:
                raise TypeError(
                    "Can't convert {} to FrozenMonosaccharideResidue".format(
                        type(args[0])))
        for residue in args:
            self[residue] += 1


class FrozenError(ValueError):
    pass
