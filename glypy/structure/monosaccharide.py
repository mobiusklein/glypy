import logging
try:
    from itertools import chain, izip_longest
except ImportError:
    from itertools import chain, zip_longest as izip_longest
from collections import deque

from six import string_types as basestring

from glypy.io.format_constants_map import anomer_map, superclass_map
from glypy.utils import (
    invert_dict,
    identity as ident_op,
    cyclewarning,
    uid,
    make_struct)
from glypy.utils.multimap import OrderedMultiMap
from glypy.utils.enum import EnumValue
from glypy.composition import Composition, calculate_mass
from glypy.composition.structure_composition import monosaccharide_composition
from glypy.composition.structure_composition import modification_compositions

from .constants import (
    Anomer, Configuration, Stem, SuperClass, Modification, RingType, UnknownPosition,
    NoPosition)
from .substituent import Substituent
from .link import Link
from .base import SaccharideBase
from .stereochemistry import stereocode


anomer_map = invert_dict(anomer_map)
superclass_map = invert_dict(superclass_map)

logger = logging.getLogger(__name__)
logger.setLevel("DEBUG")

debug = True


def _get_standard_composition(monosaccharide):
    '''Used to get initial composition for a given monosaccharide
    |Superclass| and modifications.

    Used during initialization of a |Monosaccharide|.

    Parameters
    ----------
    monosaccharide: :class:`Monosaccharide`
        The |Monosaccharide| object to read attributes from

    Returns
    -------
    :class:`~glypy.composition.composition.Composition`:
        The baseline composition from `monosaccharide.superclass` + `monosaccharide.modifications`
    '''
    base = monosaccharide_composition[monosaccharide.superclass.name]
    modifications = list(monosaccharide.modifications.items())
    double_bond_count = 0
    keto_group = False
    for mod_pos, mod_val in modifications:
        # Don't set the reducing end here
        if isinstance(mod_val, ReducedEnd):
            monosaccharide.reducing_end = mod_val
            continue
        # In case the modification is not properly transformed elsewhere
        elif mod_val is Modification.aldi:  # pragma: no cover
            monosaccharide.reducing_end = True
            continue
        elif mod_val is Modification.keto and mod_pos != 0:
            keto_group = True
        if mod_val == Modification.en:
            double_bond_count += 1
            if double_bond_count % 2 == 0:
                continue
        # Using global constant modification
        try:
            base += modification_compositions[mod_val](mod_pos)
        # Using object-oriented modification
        except Exception:  # pragma: no cover
            base += mod_val.composition

    # The presence of a keto group implicitly converts the C1 carbonyl oxygen into
    # a hydroxyl and balances C1 with another hydrogen.
    if keto_group:
        base += Composition("H2")
    return base


def traverse(monosaccharide, visited=None, apply_fn=ident_op):
    '''
    A low-level depth-first traversal method for unwrapped residue
    graphs

    Parameters
    ----------
    monosaccharide: :class:`Monosaccharide`
        Residue to start traversing from
    visited: set or None
        The collection of node ids to ignore, having already visited them. If |None|, it defaults
        to the empty set.
    apply_fn: function
        Function to apply to each residue before yielding them

    Yields
    -------
    :class:`Monosaccharide`
    '''
    if visited is None:
        visited = set()
    yield apply_fn(monosaccharide)
    visited.add(monosaccharide.id)
    outnodes = sorted(list(monosaccharide.links.items()), key=lambda x: x[1][monosaccharide].degree(), reverse=True)
    for link_pos, link in outnodes:
        child = link[monosaccharide]
        if child.id in visited:
            continue
        for grandchild in traverse(child, visited=visited, apply_fn=apply_fn):
            yield grandchild


def _traverse_debug(monosaccharide, visited=None, apply_fn=ident_op):  # pragma: no cover
    '''
    A low-level depth-first traversal method for unwrapped residue
    graphs when the :attr:`id` attribute may be masking duplicate residues

    Parameters
    ----------
    monosaccharide: :class:`Monosaccharide`
        Residue to start traversing from
    visited: set or None
        The collection of node ids to ignore, having already visited them. If |None|, it defaults
        to the empty set.
    apply_fn: function
        Function to apply to each residue before yielding them

    Yields
    -------
    :class:`Monosaccharide`
    '''
    if visited is None:
        visited = set()
    yield apply_fn(monosaccharide)
    visited.add(id(monosaccharide))
    outnodes = sorted(list(monosaccharide.links.items()), key=lambda x: x[1][monosaccharide].degree(), reverse=True)
    for link_pos, link in outnodes:
        child = link[monosaccharide]
        if id(child) in visited:
            continue
        for grandchild in _traverse_debug(child, visited=visited, apply_fn=apply_fn):
            yield grandchild


def graph_clone(monosaccharide, visited=None):
    '''
    Low-level depth-first duplication method for unwrapped residue graphs
    which contain monosaccharides with substituent parent nodes.

    Parameters
    ----------
    residue: :class:`Monosaccharide`
        The root of the graph to clone
    visited: set or None
        The collection of node ids to ignore, having already visited them. If |None|, it defaults
        to the empty set.

    Returns
    -------
    :class:`Monosaccharide`:
        The root of a newly duplicated and identical residue graph
    '''
    llen = len
    index = {}
    visited = set() if visited is None else visited
    clone_root = monosaccharide.clone(prop_id=True)
    index[clone_root.id] = clone_root
    node_stack = deque([(clone_root, monosaccharide)])
    node_stack_append = node_stack.append
    node_stack_pop = node_stack.pop
    while llen(node_stack) > 0:
        clone, ref = node_stack_pop()
        if ref.id in visited:
            continue
        visited.add(ref.id)
        links = [link for pos, link in ref.links.items()]
        for link in links:
            terminal = link[ref]
            if terminal.id in visited or terminal.node_type != Monosaccharide.node_type:
                continue
            # Handle cycles where the same node is linked many times
            if terminal.id in index:
                clone_terminal = index[terminal.id]
                # clone_terminal.maybe_cyclic = True
                cyclewarning()
            else:
                index[terminal.id] = clone_terminal = terminal.clone(prop_id=True)
            if link.is_child(terminal):
                link.clone(clone, clone_terminal)
            else:
                link.clone(clone_terminal, clone)
            node_stack_append((clone_terminal, terminal))
            i = 0
            for p, subst in clone_terminal.substituents():
                j = 0
                for sp, child in subst.children():
                    if child.node_type == Monosaccharide.node_type:
                        orig = terminal.substituent_links[p][i].child.children()[j][1]
                        node_stack_append((child, orig))
                    j += 1
                i += 1
    return clone_root


def release(monosaccharide):
    '''Break all monosaccharide-monosaccharide links on `monosaccharide`, returning
    them as a list. Breaking is done with `refund=True`

    Parameters
    ----------
    monosaccharide : Monosaccharide
        |Monosaccharide| to break all links on

    Returns
    -------
    list of tuple(link, (link.parent, link.child))
    '''
    return [(link, link.try_break_link(refund=True)) for link in list(monosaccharide.links.values())]


def toggle(monosaccharide):
    '''A simple generator for declaratively masking and masking a residue's links. The
    first iteration masks all links. The second unmasks them. Calls :func:`release`

    Parameters
    ----------
    monosaccharide : Monosaccharide
        |Monosaccharide| to mask links on
    '''
    links = release(monosaccharide)
    yield links
    for link, termini in links:
        link.try_apply()
    yield False


def depth(monosaccharide, visited=None):
    '''
    Calculate  the distance from `monosaccharide` to its furthest grand-child node.
    '''
    if visited is None:
        visited = set()
    if monosaccharide.id in visited:  # pragma: no cover
        return 0
    visited.add(monosaccharide.id)
    depth_count = 1
    if(len(monosaccharide.links) > 1):
        depth_count += max(depth(ch, visited) for p, ch in monosaccharide.children())
    return depth_count


class Monosaccharide(SaccharideBase):
    '''
    Represents a single monosaccharide molecule, and its relationships with other
    molcules through |Link| objects. |Link| objects stored in :attr:`links` for connections to other
    |Monosaccharide| instances, building a |Glycan|
    structure as a graph of |Monosaccharide| objects. |Link|
    objects connecting the |Monosaccharide| instance to |Substituent|
    objects are stored in :attr:`substituent_links`.

    Both :attr:`links` and :attr:`substituent_links` are instances of
    |OrderedMultiMap| objects where the key is the index of the
    carbon atom in the carbohydrate backbone that hosts the bond.
    An index of `x` or `-1` represents an unknown location.


    .. warning::

        While |Monosaccharide| objects expose their :attr:`.modifications`, :attr:`.links`, and
        :attr:`.substituent_links` attributes as mutable, you should treat them as **read-only**.
        The methods for altering their contents, :meth:`add_substituent`, :meth:`add_monosaccharide`,
        :meth:`add_modification`, :meth:`drop_substituent`, :meth:`drop_monosaccharide`, and
        :meth:`drop_modification` are all responsible for handling these mutations for you. |Link| methods
        like :meth:`Link.apply` and :meth:`Link.break_link` are used internally.


    Attributes
    ----------
    anomer: :class:`Anomer`
        An entry of :class:`~.Anomer` that corresponds to the linkage type
        of the carbohydrate backbone. Is an entry of a class based on :class:`Enum`
    superclass: :class:`~.SuperClass`
        An entry of :class:`~.SuperClass` that corresponds to the number of
        carbons in the carbohydrate backbone of the monosaccharide. Controls the base composition of the
        instance and the number of positions open to be linked to or modified. Is an entry of a class
        based on :class:`~.Enum`
    configuration: :class:`~.Configuration` or {'d', 'l', 'x', 'missing', None}
        An entry of |Configuration| which corresponds to the optical
        stereomer state of the instance. Is an entry of a class based on :class:`Enum`. May possess
        more than one value.
    stem: :class:`~.Stem`
        Corresponds to the bond conformation of the carbohydrate backbone. Is an entry of a class based
        on :class:`~.Enum`. May possess more than one value.
    ring_start: |int|
        The index of the carbon of the carbohydrate backbone that starts a ring. A value of ``-1``, ``'x'``, or
        |None| corresponds to an unknown start. A value of ``0`` refers to a linear chain.
    ring_end:  |int|
        The index of the carbon of the carbohydrate backbone that ends a ring. A value of ``-1``, ``'x'``, or
        |None| corresponds to an unknown ends. A value of ``0`` refers to a linear chain.
    stereocode: :class:`~.Stereocode`
        The stereochemistry of all carbons of the monosaccharide's backbone ring/chain.
    reducing_end: :class:`ReducedEnd`
        The reducing end terminal group of the monosaccharide if the monosaccharide is uncyclized
    modifications: |OrderedMultiMap|
        The mapping of sites to |Modification| entries. Directly modifies the instance's :attr:`composition`
    links: |OrderedMultiMap|
        The mapping of sites to |Link| entries that refer to other |Monosaccharide| instances
    substituent_links: |OrderedMultiMap|
        The mapping of sites to |Link| entries that refer to
        |Substituent| instances.
    composition: |Composition|
        An instance of |Composition| corresponding to the elemental composition
        of ``self`` and its immediate modifications. If not provided, this will be inferred
        from field values.
    reduced: :class:`ReducedEnd`
        An instance of ReducedEnd, or the value |True|, represents a reduced sugar. May be inferred
        from `modifications` if "aldi" is present

    '''
    _serializers = {}

    __slots__ = (
        "id", "_anomer", "_configuration", "_stem", "_superclass",
        "ring_start", "ring_end", "links", "substituent_links",
        "modifications", "composition",
        "_reducing_end", "_degree",
        "_checked_for_reduction"
    )

    def __init__(self, anomer=None, configuration=None, stem=None,
                 superclass=None, ring_start=UnknownPosition, ring_end=UnknownPosition,
                 modifications=None, links=None, substituent_links=None,
                 composition=None, reduced=None, id=None, fast=False):
        """Create a new :class:`Monosaccharide` using the provided attributes, leaving
        any unspecified attributes in an explicit unknown state.

        If ``fast`` is :const:`True`, all parameter values are assumed to be :class:`~.EnumValue`
        instances from the correct classes, and reduction checking is done quickly. This option
        is :const:`False` by default on explicit calls but set to :const:`True` for certain internal
        cases where this precondition can be ensured.

        Parameters
        ----------
        anomer : :class:`str` or :class:`~.EnumValue`, optional
            Specify the anomeric carbon state. Translates through :class:`~.Anomer`
        configuration : :class:`str` or :class:`~.EnumValue` or :class:`tuple` optional
            Specify the stereoisomeric configuration state(s) of the carbon rings. Will
            be converted to a :class:`tuple` if not provided as one, and elements will be
            translated through :class:`~.consants.Configuration`
        stem : :class:`str` or :class:`~.EnumValue` or :class:`tuple`, optional
            Specify the base carbon ring stereochemical template. Will
            be converted to a :class:`tuple` if not provided as one, and elements will be
            translated through :class:`~.consants.Stem`
        superclass : :class:`str` or :class:`~.EnumValue`, optional
            Specify the carbon chain length. Translates through :class:`~.SuperClass`.
        ring_start : int, optional
            Specify the carbon ring start position. If missing, defaults to :const:`UnknownPosition`
        ring_end : int, optional
            Specify the carbon ring closing position. If missing, defaults to :const:`UnknownPosition`
        modifications : OrderedMultiMap, optional
            Specify the collection of :class:`~.constant.Modification` values located at different carbon
            positions.
        links : OrderedMultiMap, optional
            Specify the collection of :class:`~.Link` instances connecting this monosaccharide to other
            monosaccharides
        substituent_links : OrderedMultiMap, optional
            Specify the collection of :class:`~.Link` instances connecting this monosaccharide to its
            substituents
        composition : :class:`~.Composition`, optional
            The chemical composition of the carbon ring to use prior to applying linkage changes and
            modifications. Automatically derived from `superclass` if omitted.
        reduced : :class:`ReducedEnd` or :class:`bool`, optional
            Whether the monosaccharide is reduced, or if a :class:`ReducedEnd` is given, specifies the
            reduction type.
        id : :class:`int`, optional
            The value to set for the node's :attr:`id`. If not provided, one will be generated using
            :func:`~.glypy.utils.uid`.
        fast : bool, optional
            Whether or not to use the fast initialization branch or not. Should only be used if all
            parameters have been properly typed and valid.
        """

        if modifications is None:  # pragma: no cover
            modifications = OrderedMultiMap()
        if links is None:
            links = OrderedMultiMap()

        if id is None:
            id = uid()

        self.modifications = modifications
        self._reducing_end = None
        self._checked_for_reduction = False

        if fast:
            self._anomer = anomer
            self._configuration = tuple(configuration)
            self._stem = tuple(stem)
            self._superclass = superclass
            self._fast_reduce(reduced)
        else:
            self.anomer = anomer
            self.configuration = configuration
            self.stem = stem
            self.superclass = superclass
            self.reducing_end = reduced

        self.ring_start = ring_start
        self.ring_end = ring_end
        self.links = links
        self.substituent_links = OrderedMultiMap() if substituent_links\
            is None else substituent_links
        self.id = id
        if composition is None:
            composition = _get_standard_composition(self)
        self.composition = composition
        self._degree = len(self.links) + len(self.substituent_links)

    @property
    def anomer(self):
        return self._anomer

    @anomer.setter
    def anomer(self, value):
        self._anomer = Anomer[value]

    @property
    def configuration(self):
        return self._configuration

    @configuration.setter
    def configuration(self, value):
        if isinstance(value, (tuple, list)):
            self._configuration = tuple(Configuration[v] for v in value)
        else:
            self._configuration = (Configuration[value],)

    @property
    def stem(self):
        return self._stem

    @stem.setter
    def stem(self, value):
        if isinstance(value, (tuple, list)):
            self._stem = tuple(Stem[v] for v in value)
        else:
            self._stem = (Stem[value],)

    @property
    def superclass(self):
        return self._superclass

    @superclass.setter
    def superclass(self, value):
        self._superclass = SuperClass[value]

    @property
    def stereocode(self):
        return stereocode(self)

    def __root__(self):
        return self

    def clone(self, prop_id=False, fast=True, monosaccharide_type=None):
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
        if monosaccharide_type is None:
            monosaccharide_type = Monosaccharide
        modifications = OrderedMultiMap()
        for site, modification_list in self.modifications.lists():
            for mod in modification_list:
                if isinstance(mod, ReducedEnd):
                    continue
                try:
                    modifications[site] = mod
                except Exception:  # pragma: no cover
                    modifications[site] = mod.clone()

        monosaccharide = monosaccharide_type(
            superclass=self._superclass,
            stem=self._stem,
            configuration=self._configuration,
            ring_start=self.ring_start,
            ring_end=self.ring_end,
            modifications=modifications,
            anomer=self._anomer,
            reduced=self._reducing_end.clone(prop_id=prop_id) if self._reducing_end is not None else None,
            id=self.id if prop_id else None,
            fast=fast)
        for pos, link in self.substituent_links.items():
            sub = link.to(self)
            dup = sub.clone(prop_id=prop_id)
            link.clone(monosaccharide, dup)
        return monosaccharide

    @property
    def ring_type(self):
        """The size of the ring-shape of the carbohydrate, as computed
        by :attr:`ring_end` - :attr:`ring_start`.

        Returns
        -------
        :class:`~.EnumValue`:
            The appropriate value of :class:`~.RingType`
        """
        null_positions = (UnknownPosition, NoPosition)
        if self.ring_start in null_positions or self.ring_end in null_positions:
            return RingType.x
        try:
            diff = self.ring_end - self.ring_start
            if diff == 4:
                return RingType.pyranose
            elif diff == 3:
                return RingType.furanose
            elif diff == 0:
                return RingType.open
        except TypeError:
            pass
        return RingType.x

    @property
    def reducing_end(self):
        """Return the reducing end type of
        `self` or |None| if `self` is not reduced. The reducing end value
        can also be found in :attr:`modifications`.

        If `Modification.aldi` is present, it will be converted into
        an instance of :class:`.ReducedEnd` with default arguments.

        TODO: Remove Redundancy between `aldi` check and reduction.

        Returns
        -------
        ReducedEnd or None
        """
        if self._reducing_end is None:
            if not self._checked_for_reduction:
                for pos, mod in list(self.modifications.items()):
                    if isinstance(mod, ReducedEnd):
                        self._reducing_end = mod
                        break
                    elif mod is Modification.aldi:  # pragma: no cover
                        self.modifications.popv(Modification.aldi)
                        self.reducing_end = ReducedEnd()
                        break
            else:
                return self._reducing_end
        if self._reducing_end is None:
            self._checked_for_reduction = True
        return self._reducing_end

    @reducing_end.setter
    def reducing_end(self, value):
        """Sets the reducing end type of `self`.

        If `value` is |None|, then this residue is not reduced.
        Else, if `value` ``is True``, then this residue's reducing end
        is set to an instance of :class:`.ReducedEnd` with default parameters.
        Else, this residue's reducing_end is set to `value`.

        This process will scan :attr:`modifications` for a pre-existing `reducing_end`
        value to replace.

        Parameters
        ----------
        value: True, None, or ReducedEnd

        """
        red_end = self.reducing_end
        if red_end is not None:
            self.modifications.pop(1, red_end)
        # Setting to None will un-reduce the sugar, but should not
        # put a None in :attr:`modifications`
        self._checked_for_reduction = False
        if value:
            if value is True:
                value = ReducedEnd()
            self.modifications[1] = value
            self._reducing_end = value
        else:
            self._reducing_end = value

    def _fast_reduce(self, value):
        """Expedite adding a reducing end to this monosaccharide. Assumes
        that `value` is a :class:`ReducedEnd` and that this monosaccharide is not
        already reduced.

        Parameters
        ----------
        value : ReducedEnd
            Description
        """
        if value is not None:
            self.modifications[1] = value
            self._reducing_end = value

    def __getitem__(self, position):
        '''
        Gets the collection of alterations made to the carbohydrate
        backbone at `position`. This queries :attr:`modifications`, :attr:`links`, and
        :attr:`substituent_links`.

        Returns
        -------
        :class:`dict`
        '''
        return {
            'modifications': self.modifications[position],
            'substituent_links': self.substituent_links[position],
            'links': self.links[position]
        }

    def _is_full(self, max_occupancy=0):
        return self._remaining_capacity() >= 0

    def _remaining_capacity(self):
        bonds = self.degree(deep=True)
        max_size = self.superclass.value - 1
        if self.ring_type == RingType.open:
            max_size += 1
        return max_size - bonds

    def open_attachment_sites(self, max_occupancy=0):
        '''
        When attaching :class:`Monosaccharide` instances to other objects,
        bonds are formed between the carbohydrate backbone and the other object.
        If a site is already bound, the occupying object fills that space on the
        backbone and prevents other objects from binding there.

        Currently only cares about the availability of the hydroxyl group. As there
        is not a hydroxyl attached to the ring-ending carbon, that should not be
        considered an open site.

        If any existing attached units have unknown positions, we can't provide any
        known positions, in which case the list of open positions will be a :class:`list`
        of ``-1`` s of the length of open sites.

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
        slots, unknowns = self._backbone_occupancy_list()

        null_positions = [
            UnknownPosition, 'x', NoPosition
        ]
        open_slots = []
        can_determine_positions = unknowns > 0 or (self.ring_end in null_positions)

        for i in range(len(slots)):
            if slots[i] <= max_occupancy and (i + 1) != self.ring_end:
                open_slots.append(
                    (i + 1) if not can_determine_positions else UnknownPosition)

        if self.ring_end in null_positions:
            open_slots.pop()

        return open_slots, unknowns

    def _backbone_occupancy_list(self):
        slots = [0] * self.superclass.value
        unknowns = 0
        null_positions = (UnknownPosition, 'x')
        for pos, obj in chain(self.modifications.items(),
                              self.links.items(),
                              self.substituent_links.items()):
            if obj == Modification.keto:
                continue
            if pos in null_positions:  # pragma: no cover
                unknowns += 1
            else:
                slots[pos - 1] += 1
        return slots, unknowns

    def occupied_attachment_sites(self):
        slots, unknowns = self._backbone_occupancy_list()
        occupied_sites = []
        for i, occupied in enumerate(slots, start=1):
            if occupied:
                occupied_sites.append(i)
        if self.ring_end != UnknownPosition:
            occupied_sites.append(self.ring_end)
        for i in range(unknowns):
            occupied_sites.append(UnknownPosition)
        return occupied_sites

    def total_attachement_sites(self):
        slots = self.superclass.value
        if self.ring_type != RingType.open:
            slots -= 1
        return slots

    def is_occupied(self, position):
        '''
        Checks to see if a particular backbone position is occupied by a :class:`~.constants.Modification`,
        :class:`~.substituent.Substituent`, or :class:`~.link.Link` to another :class:`Monosaccharide`.

        Parameters
        ----------
        position: int
            The position to check for occupancy. Passing -1 checks for undetermined attachments.

        Returns
        -------
        int:
            The number of occupants at ``position``

        Raises
        ------
        IndexError:
            When the position is less than 1 or exceeds the limits of the carbohydrate backbone's size.

        '''

        null_positions = (UnknownPosition, NoPosition, 'x')

        if (position > self.superclass.value) or (position < 1 and position not in null_positions):
            raise IndexError("Index out of bounds")
        # The unknown position is always available
        if position in (UnknownPosition, 'x'):  # pragma: no cover
            return 0
        n_occupants = len(self.links[position]) +\
            len(self.modifications[position]) +\
            len(self.substituent_links[position])
        if Modification.keto in self.modifications[position]:
            n_occupants -= 1
        return n_occupants

    def add_modification(self, modification, position, max_occupancy=0):
        '''
        Adds a modification instance to :attr:`modifications` at the site given
        by ``position``. This directly modifies :attr:`composition`, consequently
        changing :meth:`mass`

        Parameters
        ----------
        modification: str or Modification
            The modification to add. If passed a |str|, it will be
            translated into an instance of :class:`~glypy.structure.constants.Modification`
        position: int or 'x'
            The location to add the :class:`~.constants.Modification` to.
        max_occupancy: int, optional
            The maximum number of items acceptable at `position`. defaults to :const:`1`

        Raises
        ------
        IndexError
            ``position`` exceeds the bounds set by :attr:`superclass`.
        ValueError
            ``position`` is occupied by more than ``max_occupancy`` elements

        Returns
        -------
        :class:`Monosaccharide`:
            `self`, for chain calls
        '''
        if self.is_occupied(position) > max_occupancy:
            raise ValueError("Site is already occupied")
        self._checked_for_reduction = False
        is_keto = modification == Modification.keto
        if modification is Modification.aldi:  # pragma: no cover
            self.reducing_end = ReducedEnd()
        else:
            try:
                self.composition += modification_compositions[modification](position)
                self.modifications[position] = Modification[modification]
            # OO Modification
            except Exception:  # pragma: no cover
                self.composition += modification.composition
                self.modifications[position] = modification
        if is_keto:
            n_ketos = 0
            for position, modification in self.modifications.items():
                if modification == Modification.keto:
                    n_ketos += 1
            if n_ketos == 1:
                self.composition += Composition("H2")
        return self

    def drop_modification(self, position, modification):
        '''
        Remove the `modification` at `position`

        Parameters
        ----------
        position: int
            The position to drop the modification from
        modification: Modification
            The Modification to remove.

        Raises
        ------
        IndexError:
            If `position` is not a valid carbohydrate backbone position

        ValueError:
            If `modification` is not found at `position`

        Returns
        -------
        :class:`Monosaccharide`:
            `self`, for chain calls
        '''
        if position > self.superclass.value or position < 0 and position not in (UnknownPosition, 'x'):
            raise IndexError("Index out of bounds")
        try:
            self.modifications.pop(position, modification)
        except IndexError:
            raise ValueError("Modification {} not found at {}".format(modification, position))
        is_keto = modification == Modification.keto
        self._checked_for_reduction = False
        if modification is Modification.aldi:  # pragma: no cover
            self.reducing_end = None
        else:
            try:
                self.composition = self.composition - \
                    modification_compositions[modification](position)
            # OO Modification
            except Exception:  # pragma: no cover
                self.composition = self.composition - modification.composition
        if is_keto:
            n_ketos = 0
            for position, modification in self.modifications.items():
                if modification == Modification.keto:
                    n_ketos += 1
            if n_ketos == 0:
                self.composition -= Composition("H2")
        return self

    def add_substituent(self, substituent, position=-1, max_occupancy=0,
                        child_position=1, parent_loss=None, child_loss=None):
        '''
        Adds a :class:`~glypy.structure.substituent.Substituent` and associated :class:`~glypy.structure.link.Link`
        to :attr:`substituent_links` at the site given by ``position``. This new substituent is included when
        calculating mass with substituents included.

        >>> from glypy import monosaccharides
        >>> hex = monosaccharides.Hex
        >>> hexnac = monosaccharides.HexNAc
        >>> hex.add_substituent("n-acetyl", 2, parent_loss="OH")
        RES 1b:x-xx-HEX-1:5 2s:n-acetyl LIN 1:1d(2+1)2n
        >>> hexnac == hex
        True

        Parameters
        ----------
        substituent: str or Substituent
            The substituent to add. If passed a |str| it will be
            translated into an instance of |Substituent|.
        position: int or 'x'
            The location to add the |Substituent| link to :attr:`substituent_links`. Defaults to -1
        child_position: int
            The location to add the link to in `substituent` :attr:`links`. Defaults to -1. Substituent
            indices are currently not checked.
        max_occupancy: int, optional
            The maximum number of items acceptable at ``position``. Defaults to :const:`1`
        parent_loss: Composition or str
            The elemental composition removed from ``self``
        child_loss: Composition or str
            The elemental composition removed from ``substituent``

        Raises
        ------
        IndexError
            ``position`` exceeds the bounds set by :attr:`superclass`.
        ValueError
            ``position`` is occupied by more than ``max_occupancy`` elements


        Returns
        -------
        :class:`Monosaccharide`:
            `self`, for chain calls
        '''
        if self.is_occupied(position) > max_occupancy:
            raise ValueError("Site is already occupied")
        if isinstance(substituent, basestring):
            substituent = Substituent(substituent)
        if child_loss is None:
            child_loss = Composition("H")
        if parent_loss is None:
            parent_loss = substituent.attachment_composition_loss()
        Link(parent=self, child=substituent,
             parent_position=position, child_position=child_position,
             parent_loss=parent_loss, child_loss=child_loss)
        return self

    def drop_substituent(self, position, substituent=None, refund=True):
        '''
        Remove the `substituent` at `position`.

        If `substituent` is |None|, then the first substituent found at `position` is
        removed.

        >>> from glypy import monosaccharides
        >>> hex = monosaccharides.Hex
        >>> hexnac = monosaccharides.HexNAc
        >>> hexnac.drop_substituent(2)
        RES 1b:x-xx-HEX-1:5
        >>> hexnac == hex
        True

        Parameters
        ----------
        position: int
            The position to drop the modification from
        substituent: Substituent
            The |Substituent| to remove. If |None|, the first substituent found at
            `position` will be removed
        refund: bool
            Passed to :meth:`~.Link.break_link`

        Raises
        ------
        IndexError:
            If `position` is not a valid carbohydrate backbone position

        ValueError:
            If `substituent` is not found at `position`

        Returns
        -------
        :class:`Monosaccharide`:
            `self`, for chain calls
        '''
        if position > self.superclass.value or position < 0 and position not in (UnknownPosition, 'x'):
            raise IndexError("Index out of bounds")
        if isinstance(substituent, basestring):
            substituent = Substituent(substituent)
        link_obj = None
        for substituent_link in self.substituent_links[position]:
            if substituent is None or substituent_link.child.name == substituent.name:
                link_obj = substituent_link
                break
        if link_obj is None:
            if substituent is not None:
                msg = "No matching substituent found at {position}.".format(
                    position=position)
            else:
                msg = "No substituents found at {position}.".format(
                    position=position)
            raise IndexError(msg)

        link_obj.break_link(refund=refund)
        return self

    def add_monosaccharide(self, monosaccharide, position=-1, max_occupancy=0,
                           child_position=-1, parent_loss=None, child_loss=None):
        '''
        Adds a |Monosaccharide| and associated |Link| to :attr:`links` at the site given by
        ``position``.

        >>> from glypy import monosaccharides
        >>> hexnac = monosaccharides.HexNAc
        >>> hex = monosaccharides.Hex
        >>> hexnac.add_monosaccharide(hex, 1)
        RES 1b:x-xx-HEX-1:5 2s:n-acetyl LIN 1:1d(2+1)2n
        >>> hexnac.links[1][0].child
        RES 1b:x-xx-HEX-1:5

        Parameters
        ----------
        monosaccharide: Monosaccharide
            The monosaccharide to add.
        position: int or 'x'
            The location to add the |Monosaccharide| link to :attr:`links`. Defaults to -1
        child_position: int
            The location to add the link to in `monosaccharide`'s :attr:`links`. Defaults to -1.
        max_occupancy: int, optional
            The maximum number of items acceptable at ``position``. Defaults to :const:`1`
        parent_loss: Composition or str
            The elemental composition removed from ``self``
        child_loss: Composition or str
            The elemental composition removed from ``monosaccharide``

        Raises
        ------
        IndexError
            ``position`` exceeds the bounds set by :attr:`superclass`.
        ValueError
            ``position`` is occupied by more than ``max_occupancy`` elements


        Returns
        -------
        :class:`Monosaccharide`:
            `self`, for chain calls
        '''
        if parent_loss is None:
            parent_loss = Composition("H")
        if child_loss is None:
            child_loss = Composition("OH")
        if self.is_occupied(position) > max_occupancy:
            raise ValueError("Parent Site is already occupied")  # pragma: no cover
        if monosaccharide.is_occupied(child_position) > max_occupancy:
            raise ValueError("Child Site is already occupied")  # pragma: no cover
        Link(parent=self, child=monosaccharide,
             parent_position=position, child_position=child_position,
             parent_loss=parent_loss, child_loss=child_loss)
        return self

    def drop_monosaccharide(self, position, refund=True):
        '''
        Remove the glycosidic bond at `position`, detatching a connected |Monosaccharide|

        If there is more than one glycosidic bond at `position`, an error will be raised.

        >>> from glypy import glycans
        >>> n_linked_core = glycans["N-Linked Core"]
        >>> n_linked_core.root.drop_monosaccharide(4)
        RES 1b:b-dglc-HEX-1:5 2s:n-acetyl LIN 1:1d(2+1)2n
        >>> n_linked_core.mass()
        221.08993720321

        Parameters
        ----------
        position: int
            The position to drop the modification from
        refund: bool
            Passed to :meth:`~.Link.break_link`

        Raises
        ------
        IndexError:
            If `position` is not a valid carbohydrate backbone position

        ValueError:
            If no |Link| or more than one |Link| is found at `position`

        Returns
        -------
        :class:`Monosaccharide`:
            `self`, for chain calls
        '''
        if position > self.superclass.value or position < 0 and position not in (UnknownPosition, 'x'):
            raise IndexError("Index out of bounds")
        link_obj = None
        if len(self.links[position]) > 1:
            raise ValueError("Too many monosaccharides found")
        for link in self.links[position]:
            link_obj = link
            break
        if link_obj is None:
            raise ValueError(
                "No matching monosaccharide found at {position}".format(position=position))

        link_obj.break_link(refund=refund)
        return self

    def serialize(self, name='glycoct'):
        """Convert this object into text using the requested textual encoding

        Parameters
        ----------
        name: :class:`str`, optional
            The name of the textual encoding, e.g. "glycoct" or "iupac"

        Returns
        -------
        :class:`str`
        """
        return self._serializers[name.lower()](self)

    @classmethod
    def register_serializer(cls, name, method):
        """Add `method` as `name` to the set of serializers to pick from
        in :meth:`serialize`

        Parameters
        ----------
        name : str
            The name of the serializer
        method : Callable
            A callable object that when called with a :class:`Monosaccharide` returns a :class:`str`
        """
        cls._serializers[name.lower()] = method

    @classmethod
    def available_serializers(cls):
        """Get the list of available serialization formats

        Returns
        -------
        :class:`list` of :class:`str`
        """
        return list(cls._serializers.keys())

    def _flat_equality(self, other, lengths=True):
        '''
        Test for equality of all scalar-ish features that do not
        require recursively comparing links which in turn compare their
        connected units.
        '''
        llen = len
        flat = (self.anomer == other.anomer) and\
            (self.ring_start == other.ring_start) and\
            (self.ring_end == other.ring_end) and\
            (self.superclass == other.superclass) and\
            (self.modifications) == (other.modifications) and\
            (self.configuration) == (other.configuration) and\
            (self.stem) == (other.stem)
        if lengths:
            flat = flat and\
                llen(self.links) == llen(other.links) and\
                llen(self.substituent_links) == llen(other.substituent_links) and\
                self.total_composition() == other.total_composition()
        return flat

    def exact_ordering_equality(self, other, substituents=True, visited=None):
        '''
        Performs equality testing between two monosaccharides where
        the exact position (and ordering by sort) of links must to match between
        the input |Monosaccharide| objects

        Returns
        -------
        |bool|
        '''
        if visited is None:
            visited = set()
        if (self.id, other.id) in visited:  # pragma: no cover
            return True

        visited.add((self.id, other.id))

        if self._flat_equality(other):
            if substituents:
                for a_sub, b_sub in zip(self.substituents(), other.substituents()):
                    if a_sub != b_sub:  # pragma: no cover
                        return False
            for a_mod, b_mod in zip(self.modifications.items(), other.modifications.items()):
                if a_mod != b_mod:  # pragma: no cover
                    return False
            for a_child, b_child in zip(self.children(True), other.children(True)):
                a_parent_pos, a_link = a_child
                b_parent_pos, b_link = b_child
                if (a_parent_pos != b_parent_pos) or (a_link.child_position != b_link.child_position):
                    return False
                if not a_link.child.exact_ordering_equality(b_link.child,
                                                            substituents=substituents,
                                                            visited=visited):  # pragma: no cover
                    return False
            return True
        return False

    def topological_equality(self, other, substituents=True, visited=None):
        '''
        Performs equality testing between two monosaccharides where
        the exact ordering of child links does not have to match between
        the input |Monosaccharide|s, so long as an exact match of the
        subtrees is found

        Returns
        -------
        |bool|
        '''
        if visited is None:
            visited = set()
        if (self.id, other.id) in visited:  # pragma: no cover
            return True
        if self._flat_equality(other) and (not substituents or self._match_substituents(other)):
            taken_b = set()
            b_children = list(other.children(links=True))
            a_children = list(self.children(links=True))
            for a_pos, a_link in a_children:
                a_child = a_link.child
                matched = False
                for b_pos, b_link in b_children:
                    b_child = b_link.child
                    if (b_pos, b_child.id) in taken_b:
                        continue
                    if a_child.topological_equality(b_child,
                                                    substituents=substituents,
                                                    visited=visited):
                        matched = True
                        taken_b.add((b_pos, b_child.id))
                        break
                if not matched and len(a_children) > 0:
                    return False
            if len(taken_b) != len(b_children):  # pragma: no cover
                return False
            return True
        return False

    def _match_substituents(self, other):
        '''
        Helper method for matching substituents in an order-independent
        fashion. Used by :meth:`topological_equality`
        '''
        taken_b = set()
        b_num_unknown = len(other.substituent_links[-1])
        unknown_cntr = 0
        b_substituents = list(other.substituents())
        cntr = 0
        for a_pos, a_substituent in self.substituents():
            matched = False
            cntr += 1
            for b_pos, b_substituent in b_substituents:
                if b_pos in taken_b:
                    if b_pos != UnknownPosition:  # pragma: no cover
                        continue
                    else:  # pragma: no cover
                        if unknown_cntr < b_num_unknown:
                            unknown_cntr += 1
                        else:
                            continue
                if b_substituent == a_substituent:
                    matched = True
                    taken_b.add(b_pos)
                    break
            if not matched and cntr > 0:
                return False
        if len(taken_b) + unknown_cntr != len(b_substituents):
            return False
        return True

    def __eq__(self, other):
        '''
        Test for equality between :class:`Monosaccharide` instances.
        First try scalar equality of fields, and then compare descendants.
        '''
        if (other is None):
            return False
        if not isinstance(other, Monosaccharide):
            return False
        return self.exact_ordering_equality(other)

    def __hash__(self):
        return hash(self.id)

    def __ne__(self, other):
        return not (self == other)

    __repr__ = serialize

    def __getstate__(self):
        state = dict()
        state['_anomer'] = self._anomer
        state['_superclass'] = self._superclass
        state['_stem'] = self._stem
        state['_configuration'] = self._configuration
        state['_degree'] = self._degree
        state['ring_start'] = self.ring_start
        state['ring_end'] = self.ring_end
        state['id'] = self.id
        state['modifications'] = self.modifications
        state['links'] = self.links
        state['substituent_links'] = self.substituent_links
        state['composition'] = self.composition
        state['_reducing_end'] = self._reducing_end
        return state

    def __setstate__(self, state):
        '''
        Does some testing to upgrade outdated, but equivalent
        modification models.
        '''
        self._checked_for_reduction = False
        self._anomer = state['_anomer']
        self._superclass = state['_superclass']
        self._stem = state['_stem']
        self._configuration = state['_configuration']
        self._degree = state.get("_degree")
        self.ring_start = state['ring_start']
        self.ring_end = state['ring_end']
        self.id = state['id']

        self.modifications = state['modifications']
        self.links = state['links']
        self.substituent_links = state['substituent_links']
        self.composition = state["composition"]
        reduced = state.get('_reducing_end', None)
        if self._degree is None:
            self._degree = len(self.links) + len(self.substituent_links)
        # Make sure that if "aldi" is present, to replace it with
        # the default ReducedEnd
        if reduced is None:
            if self.modifications.popv(Modification.aldi) is not None:  # pragma: no cover
                # Deduct the modification mass from the main composition
                self.composition -= {"H": 2}
                reduced = True
        self._reducing_end = None
        self.reducing_end = reduced

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
        if charge == 0:
            mass = calculate_mass(
                self.composition, average=average, charge=0, mass_data=mass_data)
            if substituents:
                subst_visited = set()
                for substituent_link in self.substituent_links.values():
                    subst = substituent_link[self]
                    if subst.id in subst_visited:
                        continue
                    subst_visited.add(subst.id)
                    mass += subst.mass(
                        average=average, charge=0, mass_data=mass_data)
            if self._reducing_end is not None:
                mass += self._reducing_end.mass(
                    average=average, charge=0, mass_data=mass_data)
        else:
            mass = self.total_composition().calc_mass(average=average, charge=charge, mass_data=mass_data)
        return mass

    def total_composition(self):
        '''
        Computes the sum of the composition of ``self`` and each of its linked
        :class:`~glypy.structure.substituent.Substituent`\ s

        Returns
        -------
        :class:`~glypy.composition.Composition`
        '''
        comp = self.composition.clone()
        for pos, sub in self.substituents():
            comp += sub.total_composition()
        red_end = self.reducing_end
        if red_end is not None:
            comp += red_end.total_composition()
        return comp

    def children(self, links=False):
        '''
        Returns an iterator over the :class:`Monosaccharide` instancess which are considered
        the descendants of `self`

        >>> from glypy import glycans
        >>> n_linked_core = glycans["N-Linked Core"]
        >>> ch = n_linked_core.root.children()
        >>> ch[0]
        (4, RES 1b:b-dglc-HEX-1:5 2s:n-acetyl LIN 1:1d(2+1)2n)
        >>>

        Parameters
        ----------
        links: bool
            Whether to return the Link objects, or their children.
            Defaults to False

        Returns
        -------
        :class:`list` of
        position: int
            Location of the bond to the child |Monosaccharide|
        child: Monosaccharide
            |Monosaccharide| at `position`
        '''
        if links:
            result = [(pos, link) for pos, link in self.links.items() if not link.is_child(self)]
        else:
            result = [(pos, link.child) for pos, link in self.links.items() if not link.is_child(self)]
        return result

    def parents(self, links=False):
        '''
        Returns an iterator over the :class:`Monosaccharide` instances which are considered
        the ancestors of `self`.

        links: bool
            Whether to return the Link objects, or their parents.
            Defaults to False

        Returns
        -------
        :class:`list` of
        position: int
            Location of the bond to the parent |Monosaccharide|
        parent: Monosaccharide
            |Monosaccharide| at `position`
        '''
        if links:
            result = [(pos, link) for pos, link in self.links.items() if not link.is_parent(self)]
        else:
            result = [(pos, link.parent) for pos, link in self.links.items() if not link.is_parent(self)]
        return result

    def substituents(self):
        '''
        Returns an iterator over all substituents attached
        to :obj:`self` by a :class:`~.link.Link` object stored in :attr:`~.substituent_links`

        Returns
        -------
        :class:`list` of
        position: int
            Location of the bond to the substituent
        substituent: Substituent
            |Substituent| at `position`
        '''
        subst_visited = set()
        results = []
        for pos, link in self.substituent_links.items():
            subst = link.to(self)
            if subst.id in subst_visited:
                continue
            subst_visited.add(subst.id)
            results.append((pos, subst))
        return results

    def degree(self, deep=False):
        '''
        Return the "graph theory" degree of this molecule

        Returns
        -------
        int
        '''
        if deep:
            res = len(self.links) + len(self.substituent_links)
            if hasattr(self, "_degree"):
                assert res == self._degree
        return self._degree

    def __iter__(self):
        return iter(self.children())

    def has_undefined_linkages(self):
        for link in self.links.values():
            if link.parent_position == UnknownPosition or link.child_position == UnknownPosition:
                return True
        for link in self.substituent_links.values():
            if link.parent_position == UnknownPosition or link.child_position == UnknownPosition:
                return True
        return False


class ReducedEnd(object):
    """Represents the composition shift and conformation change created
    by reducing a |Monosaccharide|.

    Attributes
    ----------
    composition : |Composition|
        The elemental composition of the reducing end reduction modification.
    links: |OrderedMultiMap|
        The attached substituents
    valence: |int|
        Number of substituents this node can host
    id: |int|
        Unique identifier

    There is also a class attribute, :attr:`name` for comparison with :attr:`~.Modification.aldi`
    """

    node_type = object()
    name = 'aldi'

    def __init__(self, composition=None, substituents=None, valence=1, id=None):
        if composition is None:
            composition = Composition("H2")
        else:
            composition = Composition(composition)
        self.composition = composition
        self.base_composition = self.composition.clone()
        self.links = substituents or OrderedMultiMap()
        self.valence = valence
        self.id = id or uid()
        self._degree = len(self.links)

    def is_occupied(self, position):
        '''
        Checks to see if a particular backbone position is occupied by a or :class:`.Substituent`.

        Parameters
        ----------
        position: int
            The position to check for occupancy. Passing -1 checks for undetermined attachments.

        Returns
        -------
        numeric:
            The number of occupants at ``position``, or `float('inf')` if `position` exceeds :attr:`valence`
        '''
        if position > self.valence:  # pragma: no cover
            return float('inf')
        else:
            return len(self.links[position])

    def add_substituent(self, substituent, position=-1, max_occupancy=0,
                        child_position=1, parent_loss=None, child_loss=None):
        '''
        Adds a |Substituent| and associated |Link|
        to :attr:`substituent_links` at the site given by ``position``. This new substituent is included when
        calculating mass with substituents included

        Parameters
        ----------
        substituent: str or Substituent
            The substituent to add. If passed a |str|, it will be
            translated into an instance of |Substituent|
        position: int or 'x'
            The location to add the |Substituent| link to :attr:`substituent_links`. Defaults to -1
        child_position: int
            The location to add the link to in `substituent`'s :attr:`links`. Defaults to -1. Substituent
            indices are currently not checked.
        max_occupancy: int, optional
            The maximum number of items acceptable at ``position``. Defaults to :const:`1`
        parent_loss: Composition or str
            The elemental composition removed from ``self``
        child_loss: Composition or str
            The elemental composition removed from ``substituent``

        Raises
        ------
        IndexError
            ``position`` exceeds the bounds set by :attr:`superclass`.
        ValueError
            ``position`` is occupied by more than ``max_occupancy`` elements
        '''
        if self.is_occupied(position) > max_occupancy:  # pragma: no cover
            raise ValueError("Site is already occupied")
        if child_loss is None:
            child_loss = Composition("H")
        if parent_loss is None:
            parent_loss = Composition("H")
        if isinstance(substituent, basestring):
            substituent = Substituent(substituent)
        Link(parent=self, child=substituent,
             parent_position=position, child_position=child_position,
             parent_loss=parent_loss, child_loss=child_loss)
        return self

    def drop_substituent(self, position, substituent=None, refund=True):
        '''
        Remove the `substituent` at `position`.

        If `substituent` is |None|, then the first substituent found at `position` is
        removed.

        Parameters
        ----------
        position: int
            The position to drop the modification from
        substituent: Substituent
            The |Substituent| to remove. If |None|, the first substituent found at
            `position` will be removed
        refund: bool
            Passed to :meth:`~.Link.break_link`

        Raises
        ------
        IndexError:
            If `position` exceeds :attr:`valence`

        ValueError:
            If `substituent` is not found at `position`

        Returns
        -------
        ReducedEnd:
            `self` for chaining calls
        '''
        if position > self.valence:  # pragma: no cover
            raise IndexError("Index out of bounds")
        if isinstance(substituent, basestring):
            substituent = Substituent(substituent)
        link_obj = None
        for substituent_link in self.links[position]:
            if substituent_link.child.name == substituent.name or substituent is None:
                link_obj = substituent_link
                break
        if link_obj is None:
            if substituent is not None:
                msg = "No matching substituent found at {position}.".format(
                    position=position)
            else:
                msg = "No substituents found at {position}.".format(
                    position=position)
            raise IndexError(msg)

        link_obj.break_link(refund=refund)
        return self

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
        mass = calculate_mass(
            self.composition, average=average, charge=charge, mass_data=mass_data)

        for pos, link in self.links.items():
            mass += link[self].mass(average=average, charge=charge, mass_data=mass_data)
        return mass

    def clone(self, prop_id=True):
        """Make a deep copy of `self`.

        Parameters
        ----------
        prop_id: bool
            Whether to copy over :attr:`id`.

        Returns
        -------
        ReducedEnd
        """
        result = ReducedEnd(
            Composition(self.base_composition), substituents=None,
            valence=self.valence, id=self.id if prop_id else None)
        for position, link in self.links.items():
            result.add_substituent(link.child.clone(), position)
        return result

    def children(self):
        '''
        Returns an iterator over the nodes which are considered
        the descendants of ``self``.
        '''
        for pos, link in self.links.items():
            if link.is_child(self):
                continue
            yield (pos, link.child)

    def total_composition(self):
        '''
        Computes the sum of the composition of `self` and each of its linked
        :class:`~.substituent.Substituent`s

        Returns
        -------
        :class:`~glypy.composition.Composition`
        '''
        comp = self.composition
        for pos, sub in self.children():
            comp = comp + sub.total_composition()
        return comp

    def __repr__(self):  # pragma: no cover
        rep = "<ReducedEnd {}>".format(self.total_composition())
        return rep

    def __eq__(self, other):
        '''
        Test for equality with `other`, with special handling for `EnumValue` comparisons
        for dealing with legacy format objects. If compared to the string "aldi", returns |True|,
        otherwise does deep comparison of composition and substituents.
        '''
        if other is None:
            return False
        elif isinstance(other, (int, EnumValue)):
            return False
        elif isinstance(other, str):
            return other == self.name
        composition_equal = self.composition == other.composition
        if not composition_equal:
            return False
        for spair, opair in izip_longest(self.links.items(), other.links.items()):
            if spair[0] != opair[0]:  # pragma: no cover
                return False
            if not spair[1]._flat_equality(opair[1]):  # pragma: no cover
                return False
        return True

    def __hash__(self):
        return hash(self.name)

    def __ne__(self, other):
        return not self == other

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._degree = state.get("_degree", len(self.links))


MonosaccharideOccupancy = make_struct(
    "MonosaccharideOccupancy", ["id", "open_sites", "unknown_sites", "occupied_sites"])


def build_monosaccharide_occupancy(cls, monosaccharide):
    open_sites, unknown_sites = monosaccharide.open_attachment_sites()
    occupied_sites = monosaccharide.occupied_attachment_sites()
    return cls(monosaccharide.id, open_sites, unknown_sites, occupied_sites)


MonosaccharideOccupancy.build = classmethod(
    build_monosaccharide_occupancy)


class AnnotatedMonosaccharide(Monosaccharide):
    __slots__ = ('annotations', )

    def __init__(self, anomer=None, configuration=None, stem=None,
                 superclass=None, ring_start=UnknownPosition, ring_end=UnknownPosition,
                 modifications=None, links=None, substituent_links=None,
                 composition=None, reduced=None, id=None, fast=False, annotations=None):
        super(AnnotatedMonosaccharide, self).__init__(
            anomer, configuration, stem, superclass, ring_start, ring_end,
            modifications, links, substituent_links, composition, reduced,
            id, fast)
        if annotations is None:
            annotations = {}
        self.annotations = annotations

    def __getstate__(self):
        state = super(AnnotatedMonosaccharide, self).__getstate__()
        state['annotations'] = self.annotations
        return state

    def __setstate__(self, state):
        super(AnnotatedMonosaccharide, self).__setstate__(state)
        self.annotations = state.get("annotations", {})

    def clone(self, prop_id=False, fast=True, monosaccharide_type=None):
        '''
        Copies just this |Monosaccharide| and its |Substituent|s, creating a separate instance
        with the same data. All mutable data structures are duplicated and distinct from the original.

        Does not copy any :attr:`links` as this would cause recursive duplication of the entire |Glycan|
        graph.

        Returns
        -------
        Monosaccharide

        '''
        if monosaccharide_type is None:
            monosaccharide_type = self.__class__
        monosaccharide = super(AnnotatedMonosaccharide, self).clone(prop_id, fast, monosaccharide_type)
        monosaccharide.annotations = self.annotations.copy()
        return monosaccharide
