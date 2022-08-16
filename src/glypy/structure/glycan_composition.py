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

import warnings

try:
    from collections.abc import Mapping
except ImportError:
    from collections import Mapping

from glypy.utils import tree, uid
from glypy.utils.multimap import OrderedMultiMap

from glypy.composition import Composition
from glypy.structure.base import SaccharideCollection, MoleculeBase
from glypy.structure.glycan import Glycan
from glypy.structure.monosaccharide import Monosaccharide, ReducedEnd
from glypy.structure.substituent import Substituent
from glypy.structure.constants import (Anomer, Stem, Configuration, UnknownPosition)

from glypy.io import iupac
from glypy.io.iupac import (
    monosaccharide_reference as _monosaccharide_reference,
    resolve_special_base_type as _resolve_special_base_type,
    IUPACError)

from glypy.composition.base import formula
from glypy.composition.composition_transform import (
    derivatize, has_derivatization, strip_derivatization,
    _derivatize_reducing_end, _strip_derivatization_reducing_end,
    make_counter)

from six import string_types as basestring


monosaccharide_residue_reference = {}


class IUPACLiteMonosaccharideDeserializer(iupac.SimpleMonosaccharideDeserializer):

    def monosaccharide_from_iupac(self, monosaccharide_str, residue_class=None):
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
        if residue_class is None:
            residue_class = MonosaccharideResidue
        try:
            match_dict = self.extract_pattern(monosaccharide_str)
            residue = self.build_residue(match_dict)
        except IUPACError:
            if monosaccharide_str.startswith(MolecularComposition.sigil):
                result = MolecularComposition.from_iupac_lite(monosaccharide_str)
                return result
            if monosaccharide_str == "Pent":
                warnings.warn("Pent is not proper IUPAC. Converting to Pen")
                return self.monosaccharide_from_iupac("Pen", residue_class)
            try:
                result = SubstituentResidue.from_iupac_lite(monosaccharide_str)
                return result
            except Exception:
                try:  # pragma: no cover
                    result = MolecularComposition.from_iupac_lite(monosaccharide_str)
                    return result
                except Exception:
                    raise IUPACError("Cannot find pattern in {}".format(monosaccharide_str))
        except TypeError:
            raise TypeError("Expected string, received {} ({})".format(monosaccharide_str, type(monosaccharide_str)))

        deriv = match_dict.get("derivatization", '')
        if deriv is not None and deriv != "":
            self.apply_derivatization(residue, deriv)
        return residue_class.from_monosaccharide(residue)

    def build_residue(self, match_dict):
        residue, linkage = super(IUPACLiteMonosaccharideDeserializer, self).build_residue(match_dict)
        return residue

    def __call__(self, string, residue_class=None):
        return self.monosaccharide_from_iupac(string, residue_class=residue_class)


class IUPACLiteMonosaccharideSerializer(iupac.SimpleMonosaccharideSerializer):
    def monosaccharide_to_iupac(self, residue):
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
        try:
            string = super(IUPACLiteMonosaccharideSerializer, self).monosaccharide_to_iupac(residue)
        except (AttributeError, TypeError, ValueError):
            # if the residue passed was *really* a monosaccharide then this error is valid and
            # should propagate
            if isinstance(residue, Monosaccharide):
                raise
            else:
                string = str(residue)
        return string


from_iupac_lite = IUPACLiteMonosaccharideDeserializer()


to_iupac_lite = IUPACLiteMonosaccharideSerializer(
    iupac.monosaccharide_reference,
    iupac.SubstituentSerializer(monosaccharide_residue_reference))


Monosaccharide.register_serializer('iupac_lite', to_iupac_lite)


def drop_stem(residue, force=False):
    """Drops the stem, or the carbon ring stereochemical
    classification from this monosaccharide.

    Unless ``force`` is |True|, if :func:`~glypy.io.iupac.resolve_special_base_type`
    returns a truthy value, this function will do nothing.

    Parameters
    ----------
    residue : :class:`~.Monosaccharide`
        The monosaccharide to change
    force : bool, optional
        Whether or not to override known special case named monosaccharides

    Returns
    -------
    :class:`~.Monosaccharide`
        The mutated monosaccharide
    """
    if _resolve_special_base_type(residue) is None or force:
        residue.stem = (None,)
    return residue


def drop_positions(residue, force=False):
    """Drops the position classifiers from all links and modifications
    attached to this monosaccharide.

    Unless ``force`` is |True|, if :func:`~.iupac.resolve_special_base_type`
    returns a truthy value, this function will do nothing.

    Parameters
    ----------
    residue : :class:`~.Monosaccharide`
        The monosaccharide to change
    force : bool, optional
        Whether or not to override known special case named monosaccharides

    Returns
    -------
    :class:`~.Monosaccharide`
        The mutated monosaccharide
    """
    if _resolve_special_base_type(residue) is None or force:
        modifications = OrderedMultiMap()
        for _k, v in residue.modifications.items():
            modifications[UnknownPosition] = v
        residue.modifications = modifications

        for _p, link in list(residue.substituent_links.items()):
            link.break_link(refund=True)
            link.parent_position = UnknownPosition
            link.apply()
    return residue


def drop_configuration(residue, force=False):
    """Drops the absolute stereochemical configuration of this
    monosaccharide.

    Unless ``force`` is |True|, if :func:`~.iupac.resolve_special_base_type`
    returns a truthy value, this function will do nothing.

    Parameters
    ----------
    residue : :class:`~.Monosaccharide`
        The monosaccharide to change
    force : bool, optional
        Whether or not to override known special case named monosaccharides

    Returns
    -------
    :class:`~.Monosaccharide`
        The mutated monosaccharide
    """
    if _resolve_special_base_type(residue) is None or force:
        residue.configuration = (None,)
    return residue


water_composition = Composition({"O": 1, "H": 2})


class ResidueBase(object):
    __slots__ = ()

    def drop_stem(self, force=False):
        """Drops the stem, or the carbon ring stereochemical
        classification from this monosaccharide.

        Unless ``force`` is |True|, if :func:`~.iupac.resolve_special_base_type`
        returns a truthy value, this function will do nothing.

        Parameters
        ----------
        residue : :class:`~.Monosaccharide`
            The monosaccharide to change
        force : bool, optional
            Whether or not to override known special case named monosaccharides

        Returns
        -------
        :class:`~.Monosaccharide`
            The mutated monosaccharide
        """
        return self

    def drop_positions(self, force=False):
        """Drops the position classifiers from all links and modifications
        attached to this monosaccharide.

        Unless ``force`` is |True|, if :func:`~.iupac.resolve_special_base_type`
        returns a truthy value, this function will do nothing.

        Parameters
        ----------
        residue : :class:`~.Monosaccharide`
            The monosaccharide to change
        force : bool, optional
            Whether or not to override known special case named monosaccharides

        Returns
        -------
        :class:`~.Monosaccharide`
            The mutated monosaccharide
        """
        return self

    def drop_configuration(self, force=False):
        """Drops the absolute stereochemical configuration of this
        monosaccharide.

        Unless ``force`` is |True|, if :func:`~.iupac.resolve_special_base_type`
        returns a truthy value, this function will do nothing.

        Parameters
        ----------
        residue : :class:`~.Monosaccharide`
            The monosaccharide to change
        force : bool, optional
            Whether or not to override known special case named monosaccharides

        Returns
        -------
        :class:`~.Monosaccharide`
            The mutated monosaccharide
        """
        return self

    def to_iupac_lite(self):
        """Encode this residue using `iupac_lite` notation.

        Returns
        -------
        str
        """
        return to_iupac_lite(self)

    @classmethod
    def from_iupac_lite(cls, string):
        """Parse a string of `iupac_lite` notation to produce a residue object

        Parameters
        ----------
        string : :class:`str`
            The string to parse

        Returns
        -------
        ResidueBase
        """
        return from_iupac_lite(string, residue_class=cls)


class MonosaccharideResidue(Monosaccharide, ResidueBase):
    """Represents a :class:`Monosaccharide`-like object, save that it does
    not connect to other :class:`~.Monosaccharide` objects and does not have
    properties related to topology, specifically, :attr:`anomer`.

    A single :class:`MonosaccharideResidue` has lost a water molecule from its
    composition, reflecting its residual nature. This is accounted for when dealing
    with aggreates of residues. They also have altered carbon backbone occupancies.

    :class:`MonosaccharideResidue` objects are hashable and comparable on their
    `iupac_lite` representation, which is given by :meth:`__str__` or :meth:`name`.
    """
    __slots__ = ()

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
                residue.configuration = (Configuration.x,)
            if not stem:
                residue.stem = (Stem.x,)
        if not ring:
            residue.ring_start = residue.ring_end = UnknownPosition
        if deriv:
            derivatize(residue, deriv)
        if residue.mass() != premass and not deriv:
            residue.composition += water_composition
        return residue

    def __init__(self, *args, **kwargs):
        super(MonosaccharideResidue, self).__init__(*args, **kwargs)
        self.composition -= water_composition
        self.anomer = Anomer.x

    def clone(self, *args, **kwargs):
        '''
        Copies just this |Monosaccharide| and its |Substituent| objects, creating a separate instance
        with the same data. All mutable data structures are duplicated and distinct from the original.

        Does not copy any :attr:`links` as this would cause recursive duplication of the entire |Glycan|
        graph.

        Parameters
        ----------
        prop_id: :class:`bool`
            Whether to copy :attr:`id` from ``self`` to the new instance
        fast: :class:`bool`
            Whether to use the fast-path initialization process in :meth:`MonosaccharideResidue.__init__`
        monosaccharide_type: :class:`type`
            A subclass of :class:`MonosaccharideResidue` to use

        Returns
        -------
        :class:`MonosaccharideResidue`

        '''
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
        '''
        When attaching :class:`~.Monosaccharide` instances to other objects,
        bonds are formed between the carbohydrate backbone and the other object.
        If a site is already bound, the occupying object fills that space on the
        backbone and prevents other objects from binding there.

        Currently only cares about the availability of the hydroxyl group. As there
        is not a hydroxyl attached to the ring-ending carbon, that should not be
        considered an open site.

        If any existing attached units have unknown positions, we can't provide any
        known positions, in which case the list of open positions will be a :class:`list`
        of ``-1`` s of the length of open sites.

        A :class:`MonosaccharideResidue` has two fewer open attachment sites than
        the equivalent :class:`~.Monosaccharide`

        Parameters
        ----------
        max_occupancy: int
            The number of objects that may already be bound at a site before it
            is considered unavailable for attachment.

        Returns
        -------
        :class:`list`:
            The positions open for binding
        :class:`int`:
            The number of bound but unknown locations on the backbone.
        '''
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
        if not isinstance(other, (MonosaccharideResidue, str)):
            return False
        return str(self) == str(other)

    def name(self):
        """Name this object according to `iupac_lite`.

        Returns
        -------
        str

        See Also
        --------
        :meth:`to_iupac_lite`
        """
        return to_iupac_lite(self)

    def residue_name(self):
        """Name this object according to `iupac_lite`, omitting any derivatization

        Returns
        -------
        str

        See Also
        --------
        :meth:`to_iupac_lite`
        :meth:`name`
        """
        name = self.name()
        return name.split("^")[0]

    drop_stem = drop_stem
    drop_positions = drop_positions
    drop_configuration = drop_configuration

    def copy_underivatized(self):
        """Create a copy of this residue without derivatization.

        Returns
        -------
        :class:`MonosaccharideResidue`
        """
        return from_iupac_lite.strip_derivatization(str(self), residue_class=self.__class__)

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
    __slots__ = ("_frozen", "_total_composition", "_hash", "_name", "_mass")

    _attribute_caching_slots = (
        '_total_composition', '_hash',
        '_mass'
    )

    # _frozen = False
    # _total_composition = None
    __cache = {}

    @classmethod
    def from_monosaccharide(cls, monosaccharide, *args, **kwargs):
        inst = super(FrozenMonosaccharideResidue, cls).from_monosaccharide(monosaccharide, *args, **kwargs)
        if str(inst) not in inst.get_cache():
            inst.get_cache()[str(inst)] = inst
            inst._frozen = True
        else:
            inst = inst.get_cache()[str(inst)]
        return inst

    def __init__(self, *args, **kwargs):
        self._total_composition = None
        self._mass = None
        # _name is left undefined to use a fast-path in __str__ by not testing for
        # presence first.
        # self._name = None
        # self._hash = None
        super(FrozenMonosaccharideResidue, self).__init__(*args, **kwargs)
        self._frozen = kwargs.get("_frozen", False)

    def __setattr__(self, key, value):
        try:
            is_frozen = self._frozen
        except AttributeError:
            is_frozen = False
        if is_frozen and key not in FrozenMonosaccharideResidue._attribute_caching_slots:
            self.get_cache().pop(str(self), None)
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
        try:
            return self._hash
        except AttributeError:
            self._hash = hash(str(self))
            return self._hash

    def _update_hash(self):
        self._hash = hash(str(self))
        return self._hash

    def __eq__(self, other):
        '''
        Test for equality between :class:`MonosaccharideResidue` instances by comparing
        the result of :meth:`MonosaccharideResidue.name` calls between `self` and `other`.

        :meth:`MonosaccharideResidue.name` is an alias of :func:`to_iupac_lite` called on `self`
        '''
        if isinstance(other, MonosaccharideResidue):
            try:
                return self._name == other._name
            except AttributeError:
                return str(self) == str(other)
        elif isinstance(other, str):
            return str(self) == other
        if (other is None):
            return False
        if not isinstance(other, (MonosaccharideResidue, str)):
            return False

    def _save_to_cache(self):
        self.get_cache()[str(self)] = self

    def __str__(self):
        try:
            return self._name
        except AttributeError:
            name = to_iupac_lite(self)
            self._name = name
            return name

    def clone(self, *args, **kwargs):
        '''
        Copies just this |Monosaccharide| and its |Substituent|s, creating a separate instance
        with the same data. All mutable data structures are duplicated and distinct from the original.

        Does not copy any :attr:`links` as this would cause recursive duplication of the entire |Glycan|
        graph.

        Parameters
        ----------
        prop_id: :class:`bool`
            Whether to copy :attr:`id` from ``self`` to the new instance
        fast: :class:`bool`
            Whether to use the fast-path initialization process in :meth:`Monosaccharide.__init__`
        monosaccharide_type: :class:`type`
            A subclass of :class:`Monosaccharide` to use

        Returns
        -------
        :class:`Monosaccharide`

        '''
        if self._frozen and kwargs.get(
                "monosaccharide_type",
                FrozenMonosaccharideResidue) is FrozenMonosaccharideResidue:
            return self
        else:
            return super(FrozenMonosaccharideResidue, self).clone(*args, **kwargs)

    def __getstate__(self):
        state = super(FrozenMonosaccharideResidue, self).__getstate__()
        state['_name'] = str(self)
        state['_total_composition'] = self.total_composition()
        return state

    def __setstate__(self, state):
        self._frozen = False
        self._total_composition = state.get('_total_composition')
        self._name = state.get('_name')
        self._hash = hash(str(self))
        super(FrozenMonosaccharideResidue, self).__setstate__(state)

    @classmethod
    def get_cache(cls):
        return cls.__cache

    @classmethod
    def from_iupac_lite(cls, string):
        cache = cls.get_cache()
        try:
            return cache[string]
        except KeyError:
            result = from_iupac_lite(string, residue_class=cls)
            if string not in cache:
                for k, v in cache.items():
                    if v == result:
                        cache[string] = v
                        break
                else:
                    cache[string] = result
            return result

    def total_composition(self):
        if self._frozen:
            if self._total_composition is None:
                self._total_composition = super(FrozenMonosaccharideResidue, self).total_composition()
                self._mass = None
            return self._total_composition
        else:
            return super(FrozenMonosaccharideResidue, self).total_composition()

    def mass(self, average=False, charge=0, mass_data=None, substituents=True):
        '''
        Calculates the total mass of ``self``.

        Parameters
        ----------
        average: bool, optional, defaults to False
            Whether or not to use the average isotopic composition when calculating masses.
            When ``average == False``, masses are calculated using monoisotopic mass.
        charge: int, optional, defaults to 0
            If charge is non-zero, m/z is calculated, where m is the theoretical mass, and z is ``charge``
        mass_data: dict, optional
            If mass_data is None, standard NIST mass and isotopic abundance data are used. Otherwise the
            contents of mass_data are assumed to contain elemental mass and isotopic abundance information.
            Defaults to :const:`None`.
        substituents: bool, optional, defaults to True
            Whether or not to include substituents' masses.
        Returns
        -------
        :class:`float`

        See also
        --------
        :func:`glypy.composition.composition.calculate_mass`
        '''
        if not average and charge == 0 and mass_data is None and self._frozen:
            if self._mass is None:
                self._mass = self.total_composition().calc_mass()
            return self._mass
        return self.total_composition().calc_mass(average=average, charge=charge, mass_data=mass_data)


class SubstituentResidue(Substituent, ResidueBase):
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

    def __init__(self, name, composition=None, id=None, links=None,
                 can_nh_derivatize=None, is_nh_derivatizable=None, derivatize=False,
                 attachment_composition=None):
        if name.startswith(SubstituentResidue.sigil):
            name = name[1:]
        elif name.startswith(MolecularComposition.sigil):
            raise TypeError("Invalid Sigil. SubstituentResidue instances must be given names with either"
                            " no sigil prefix or with '@'")
        super(SubstituentResidue, self).__init__(
            name=name, composition=composition, links=links, id=id,
            can_nh_derivatize=can_nh_derivatize, is_nh_derivatizable=is_nh_derivatizable,
            derivatize=derivatize, attachment_composition=attachment_composition)

        self._residue_name = SubstituentResidue.sigil + self._name
        self.composition -= self.attachment_composition
        self.composition -= {"H": 1}
        self._hash = None

    def __hash__(self):  # pragma: no cover
        """Obtain a hash value from `self` based on :attr:`name`.

        Returns
        -------
        int
        """
        try:
            if self._hash is None:
                self._hash = hash(self._residue_name)
            return self._hash
        except AttributeError:
            return hash(self._residue_name)

    def __getstate__(self):
        state = super(SubstituentResidue, self).__getstate__()
        state['_residue_name'] = self._residue_name
        return state

    def __setstate__(self, state):
        super(SubstituentResidue, self).__setstate__(state)
        self._residue_name = state.get("_residue_name")

    def to_iupac_lite(self):
        return self._residue_name

    __str__ = to_iupac_lite

    def __repr__(self):  # pragma: no cover
        return "SubstituentResidue(%s)" % self._residue_name

    @classmethod
    def from_iupac_lite(cls, name):
        return cls(name)

    def __eq__(self, other):
        if (other is None):
            return False
        if isinstance(other, str):
            return other == self._residue_name
        if not isinstance(other, SubstituentResidue):
            return False
        return self.name == other.name

    def __ne__(self, other):  # pragma: no cover
        return not self == other

    def _backsolve_original_composition(self):
        comp = super(SubstituentResidue, self)._backsolve_original_composition()
        comp += {"H": 1}
        return comp

    def copy_underivatized(self):
        inst = self.clone()
        strip_derivatization(inst)
        return inst


class MolecularComposition(MoleculeBase, ResidueBase):  # pragma: no cover
    sigil = "#"

    __slots__ = ('name', 'composition', '_hash')

    def __init__(self, name, composition):
        self.name = name
        self.composition = composition
        self._hash = None

    def mass(self, average=False, charge=0, mass_data=None):
        return self.composition.calc_mass(average=average, charge=charge, mass_data=mass_data)

    def __repr__(self):
        return "%s%s%s%s" % (
            self.sigil, self.name, self.sigil,
            formula(self.composition))

    to_iupac_lite = __repr__

    def __reduce__(self):
        return self.__class__, (self.name, self.composition)

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

    def __hash__(self):  # pragma: no cover
        """Obtain a hash value from `self` based on :attr:`name`.

        Returns
        -------
        int
        """
        try:
            if self._hash is None:
                self._hash = hash(self.name)
            return self._hash
        except AttributeError:
            return hash(self.name)

    def __eq__(self, other):
        try:
            return self.name == other or self.name == other.name
        except AttributeError:
            return self.name == str(other)

    def __ne__(self, other):
        return not (self == other)


class _CompositionBase(dict):

    def _setitem_fast(self, key, value):
        dict.__setitem__(self, key, value)

    def _getitem_fast(self, key):
        try:
            return dict.__getitem__(self, key)
        except KeyError:
            return 0

    def __reduce__(self):
        return self.__class__, (), self.__getstate__()

    def __getstate__(self):
        d = {
            'mapping': dict(self),
            'reducing_end': self._reducing_end,
            'composition_offset': self._composition_offset
        }
        return d

    def __setstate__(self, state):
        self.update(state['mapping'])
        self._reducing_end = state['reducing_end']
        self._composition_offset = state['composition_offset']

    @classmethod
    def _empty(cls):
        inst = cls.__new__(cls)
        inst._composition_offset = water_composition.clone()
        inst._reducing_end = None
        inst._mass = None
        return inst

    def _update_from_typed_map(self, template, copy_nodes=False):
        if copy_nodes:
            for name, count in template.items():
                self._setitem_fast(name.clone(), count)
        else:
            for name, count in template.items():
                self._setitem_fast(name, count)
        reduced = template.reducing_end
        if reduced is not None:
            self.reducing_end = reduced.clone()
        self._mass = None

    def serialize(self):
        """Convert a glycan composition into a curly brace-enclosed string specifying
        pairs of `iupac_lite` and a integer count.

        If the glycan is reduced, it will be appended to the closing brace following a
        `$` character.

        Returns
        -------
        str
        """
        form = "{%s}" % '; '.join("{}:{}".format(str(k), v) for k, v in sorted(
            self.items(), key=lambda x: (x[0].mass(), str(x[0]))) if v != 0)
        reduced = self._reducing_end
        if reduced is not None:
            form = "%s$%s" % (form, formula(reduced.total_composition()))
        return form

try:
    from glypy._c.structure.glycan_composition import _CompositionBase
except ImportError:
    pass


class GlycanComposition(_CompositionBase, SaccharideCollection):
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
    reducing_end : ReducedEnd
        Describe the reducing end of the aggregate without binding it to a specific monosaccharide.
        This will contribute to composition and mass calculations.
    _composition_offset: CComposition
        Account for the one water molecule's worth of composition left over from applying the "residue"
        transformation to each monosaccharide in the aggregate.
    """
    _monosaccharide_type = MonosaccharideResidue

    _key_parser = staticmethod(from_iupac_lite)

    @classmethod
    def from_glycan(cls, glycan):
        """
        Convert a |Glycan| into a |GlycanComposition|.

        Parameters
        ----------
        glycan : :class:`~.Glycan`
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
        """Initialize a :class:`GlycanComposition` using the provided objects or keyword
        arguments, imitating the :class:`dict` initialization signature.

        If a :class:`Mapping` is provided as a positional argument, it will be used as a
        template. If arbitrary keyword arguments are provided, they will be interpreted
        using :meth:`update`. As a special case, if another :class:`GlycanComposition` is
        provided, its :attr:`reducing_end` attribute will also be copied.

        Parameters
        ----------
        *args:
            Arbitrary positional arguments
        **kwargs:
            Arbitrary keyword arguments
        """
        # dict.__init__ just calls C update method. Expensive parameter parsing
        # _CompositionBase.__init__(self)
        self._reducing_end = None
        self._mass = None
        if args or kwargs:
            self.update(*args, **kwargs)
        if args:
            template = args[0]
            if isinstance(template, GlycanComposition):
                reduced = template.reducing_end
                if reduced is not None:
                    self.reducing_end = reduced.clone()
                self._composition_offset = template._composition_offset.clone()
            else:
                self._composition_offset = water_composition.clone()
        else:
            self._composition_offset = water_composition.clone()

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
            key = self._key_parser(key)
        if key.node_type is Monosaccharide.node_type and key.reducing_end is not None:
            self.reducing_end = key.reducing_end
            key = key.clone()
            key.reducing_end = None
        _CompositionBase.__setitem__(self, key, int(value))
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
            key = self._key_parser(key)
        try:
            return _CompositionBase.__getitem__(self, key)
        except KeyError:
            return 0

    def __delitem__(self, key):
        if isinstance(key, basestring):
            key = self._key_parser(key)
        _CompositionBase.__delitem__(self, key)
        self._mass = None

    def mass(self, average=False, charge=0, mass_data=None):
        '''
        Calculates the total mass of ``self``.

        Parameters
        ----------
        average: bool, optional, defaults to False
            Whether or not to use the average isotopic composition when calculating masses.
            When ``average == False``, masses are calculated using monoisotopic mass.
        charge: int, optional, defaults to 0
            If charge is non-zero, m/z is calculated, where m is the theoretical mass, and z is ``charge``
        mass_data: dict, optional
            If mass_data is None, standard NIST mass and isotopic abundance data are used. Otherwise the
            contents of mass_data are assumed to contain elemental mass and isotopic abundance information.
            Defaults to :const:`None`.

        Returns
        -------
        :class:`float`

        See also
        --------
        :func:`glypy.composition.composition.calculate_mass`
        '''
        if self._mass is not None and charge == 0:
            return self._mass
        if charge == 0:
            mass = self._composition_offset.mass
            for residue_type, count in self.items():
                mass += residue_type.mass(average=average, charge=0, mass_data=mass_data) * count
            if self._reducing_end is not None:
                mass += self._reducing_end.mass(average=average, charge=0, mass_data=mass_data)
            self._mass = mass
        else:
            mass = self.total_composition().calc_mass(average=average, charge=charge, mass_data=mass_data)
        return mass

    def update(self, *args, **kwargs):
        if len(args) == 1:
            if isinstance(args[0], Mapping):
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

        return self.__class__(prod)

    def __rmul__(self, other):
        return self * other

    def __eq__(self, other):
        if isinstance(other, basestring):
            return str(self) == other
        if not isinstance(other, Mapping):
            return False
        self_items = set([i for i in self.items() if i[1]])
        other_items = set([i for i in other.items() if i[1]])
        return self_items == other_items

    def __ne__(self, other):
        return not (self == other)

    def __neg__(self):
        return -1 * self

    def __missing__(self, key):
        return 0

    def __contains__(self, key):
        if isinstance(key, basestring):
            key = self._key_parser(key)
        return _CompositionBase.__contains__(self, key)

    def drop_stems(self):
        for t in self:
            drop_stem(t)
        self.collapse()
        return self

    def drop_positions(self):
        for t in self:
            drop_positions(t)
        self.collapse()
        return self

    def drop_configurations(self):
        for t in self:
            drop_configuration(t)
        self.collapse()
        return self

    def total_composition(self):
        '''
        Computes the sum of the composition of all |Monosaccharide| objects in ``self``

        Returns
        -------
        :class:`~glypy.composition.Composition`
        '''

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
        return self

    def query(self, query, exact=True, **kwargs):
        """Return the total count of all residues in `self` which
        match `query` using :func:`glypy.io.nomenclature.identity.is_a`

        Parameters
        ----------
        query : :class:`~.MonosaccharideResidue` or :class:`str`
            A monosaccharide residue or a string which will be converted into one by
            :func:`from_iupac_lite` to test for an `is-a` relationship with.
        exact : bool, optional
            Passed to :func:`~.is_a`. Explicitly |True| by default
        **kwargs
            Passed to :func:`~.is_a`

        Returns
        -------
        int
            The total count of all residues which satisfy the `is-a` relationship

        See Also
        --------
        :func:`glypy.io.nomenclature.identity.is_a`

        """
        from glypy.io.nomenclature.identity import is_a
        if isinstance(query, basestring):
            query = self._key_parser(query)
        count = 0
        for key, value in self.items():
            if is_a(key, query, exact=exact, **kwargs):
                count += value
        return count

    def reinterpret(self, references, exact=True, **kwargs):
        """Aggregate the counts of all residues in `self` for each
        monosaccharide in `references` satisfying an `is-a` relationship,
        collapsing multiple residues to a single key. Any residue not
        aggregated will be preserved as-is.

        .. note::
            The order of ``references`` matters as any residue matched by
            a reference will not be considered for later references.

        Parameters
        ----------
        references : :class:`Iterable` of :class:`~.MonosaccharideResidue`
            The monosaccharides with which to test for an `is-a` relationship
        exact : bool, optional
            Passed to :func:`~.is_a`. Explicitly |True| by default
        **kwargs
            Passed to :func:`~.is_a`

        Returns
        -------
        :class:`~.GlycanComposition`
            self after key collection and collapse
        """
        from glypy.io.nomenclature.identity import is_a
        new_counts = []
        pairs = list(self.items())
        remaining_pairs = []
        for ref in references:
            count = 0
            for key, value in pairs:
                if is_a(key, ref, exact=exact, **kwargs):
                    count += value
                else:
                    remaining_pairs.append((key, value))
            if count > 0:
                new_counts.append((ref, count))
            pairs = remaining_pairs
            remaining_pairs = []
        self.clear()
        for key, value in new_counts:
            self[key] = value
        for key, value in pairs:
            self[key] = value
        return self

    @property
    def reducing_end(self):
        return self._reducing_end

    @reducing_end.setter
    def reducing_end(self, value):
        self._invalidate()
        self._reducing_end = value

    def set_reducing_end(self, value):
        self._invalidate()
        self._reducing_end = value

    def _invalidate(self):
        self._mass = None

    @property
    def composition_offset(self):
        return self._composition_offset

    @composition_offset.setter
    def composition_offset(self, value):
        self._invalidate()
        self._composition_offset = value


    def clone(self, propogate_composition_offset=True, copy_nodes=True):
        dup = self._empty()
        dup._update_from_typed_map(self, copy_nodes=copy_nodes)
        if not propogate_composition_offset:
            dup._composition_offset = Composition('H2O')
        else:
            dup._composition_offset = self._composition_offset.clone()
        return dup

    # inheriting from dict overwrites MoleculeBase.copy
    def copy(self, *args, **kwargs):
        return self.clone(*args, **kwargs)

    def __str__(self):
        return self.serialize()

    @classmethod
    def _get_parse_tokens(cls, string):
        string = str(string)
        parts = string.split('$')
        if len(parts) == 1:
            tokens = parts[0]
            reduced = None
        elif len(parts) == 2:
            tokens, reduced = parts
        else:
            raise ValueError("Could not interpret %r" % string)
        tokens = tokens[1:-1].split('; ')
        return tokens, reduced

    def _handle_reduction_and_derivatization(self, reduced, deriv):
        if reduced:
            reduced = ReducedEnd(Composition(reduced))
            self.reducing_end = reduced
        if deriv:
            self._derivatized(deriv.clone(), make_counter(uid()), include_reducing_end=False)

    @classmethod
    def parse(cls, string):
        """Parse a :class:`str` into a :class:`GlycanComposition`.

        This will parse the format produced by :meth:`serialize`

        Parameters
        ----------
        string : :class:`str`
            The string to parse

        Returns
        -------
        :class:`GlycanComposition`
        """
        tokens, reduced = cls._get_parse_tokens(string)
        inst = cls._empty()
        deriv = None
        for token in tokens:
            try:
                residue, count = _parse_name_count(token)
            except ValueError:
                if string == "{}":
                    return inst
                else:
                    raise ValueError("Malformed Token, %s" % (token,))
            key = cls._key_parser(residue)
            if "^" in residue:
                _deriv = has_derivatization(key)
                if _deriv:
                    deriv = _deriv
            inst._setitem_fast(key, count)
        inst._handle_reduction_and_derivatization(reduced, deriv)
        return inst

    def _derivatized(self, substituent, id_base, include_reducing_end=True):
        n = 2
        items = list(self.items())
        self.clear()
        for k, v in items:
            self._setitem_fast(k, v)
            if k.node_type is Substituent.node_type:
                n -= v
        self._composition_offset += (
            substituent.total_composition() -
            substituent.attachment_composition_loss() * 2) * n
        if self._reducing_end is not None and include_reducing_end:
            _derivatize_reducing_end(self._reducing_end, substituent, id_base)
        self.collapse()
        self._invalidate()

    def _strip_derivatization(self):
        self._composition_offset = Composition("H2O")
        if self._reducing_end is not None:
            _strip_derivatization_reducing_end(self._reducing_end)
        self.collapse()
        self._invalidate()

    def _invalidate(self):
        self._mass = None


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

    _monosaccharide_type = FrozenMonosaccharideResidue

    _key_parser = staticmethod(FrozenMonosaccharideResidue.from_iupac_lite)

    def __setitem__(self, key, value):
        key = self._key_parser(str(key))
        _CompositionBase.__setitem__(self, key, value)
        self._invalidate()

    def __getitem__(self, key):
        if not isinstance(key, FrozenMonosaccharideResidue):
            key = self._key_parser(str(key))
        return _CompositionBase.__getitem__(self, key)

    def __delitem__(self, key):
        key = self._key_parser(str(key))
        _CompositionBase.__delitem__(self, key)
        self._invalidate()

    @classmethod
    def parse(cls, string):
        tokens, reduced = cls._get_parse_tokens(string)
        inst = cls._empty()
        deriv = None
        key_parser = cls._key_parser
        for token in tokens:
            try:
                residue, count = _parse_name_count(token)
            except ValueError:
                if string == "{}":
                    return inst
                else:
                    raise ValueError("Malformed Token, %s" % (token,))
            key = key_parser(residue)
            if "^" in residue:
                _deriv = has_derivatization(key)
                if _deriv:
                    deriv = _deriv
            inst._setitem_fast(key, count)
        inst._handle_reduction_and_derivatization(reduced, deriv)
        return inst

    def serialize(self):
        if self._str is None:
            self._str = super(FrozenGlycanComposition, self).serialize()
        return self._str

    __str__ = serialize

    def __contains__(self, key):
        if isinstance(key, basestring):
            key = self._key_parser(key)
        return _CompositionBase.__contains__(self, key)

    def thaw(self):
        """Convert this :class:`FrozenGlycanComposition` into a :class:`GlycanComposition`
        that is not frozen.

        Returns
        -------
        :class:`GlycanComposition`
        """
        return GlycanComposition.parse(self)

    def extend(self, *args):
        if not isinstance(args[0], FrozenMonosaccharideResidue):
            if isinstance(args[0], (Monosaccharide)):
                args = map(FrozenMonosaccharideResidue.from_monosaccharide, args)
            elif isinstance(args[0], Glycan):
                args = map(
                    FrozenMonosaccharideResidue.from_monosaccharide,
                    [node for node in args[0]
                     if node.node_type is FrozenMonosaccharideResidue.node_type])
            else:
                raise TypeError(
                    "Can't convert {} to FrozenMonosaccharideResidue".format(
                        type(args[0])))
        for residue in args:
            self[residue] += 1

    def _validate(self):
        '''Populate the caching fields used for common behaviors, e.g.
        mass and string representation.
        '''
        if self._mass is None:
            self.mass()
        if self._str is None:
            self.serialize()

    def _invalidate(self):
        '''Clear the caching fields, forcing them to all be recalculated
        when next requested.
        '''
        self._mass = None
        self._str = None
        self._total_composition = None

    def clone(self, propogate_composition_offset=True, copy_nodes=False):
        dup = self._empty()
        dup._update_from_typed_map(self, copy_nodes=copy_nodes)
        if not propogate_composition_offset:
            dup._composition_offset = Composition('H2O')
        else:
            dup._composition_offset = self._composition_offset.clone()
        return dup


class FrozenError(ValueError):
    pass


class HashableGlycanComposition(FrozenGlycanComposition):
    def __str__(self):
        self._validate()
        # Directly use internal cache variable to save time calling
        # the super method chain
        return self._str

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        if isinstance(other, HashableGlycanComposition):
            if self._str is not None:
                if other._str is not None:
                    return self._str == other._str
                return self._str == str(other)
            return str(other) == str(self)
        else:
            return super(HashableGlycanComposition, self).__eq__(other)


def _parse_name_count(string):
    name, count = string.split(":")
    count = int(count)
    return name, count


try:
    _has_c = True
    from glypy._c.utils import get_parse_tokens, parse_name_count as _parse_name_count
    GlycanComposition._get_parse_tokens = get_parse_tokens
except ImportError:
    _has_c = False
