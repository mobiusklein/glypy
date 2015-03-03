from uuid import uuid4
from itertools import chain

from .constants import Anomer, Configuration, Stem, SuperClass, Modification, ReducingEnd, RingType
from .substituent import Substituent
from .link import Link
from .base import SaccharideBase

from ..io.format_constants_map import anomer_map, superclass_map
from ..utils import invert_dict, make_counter, StringIO, identity as ident_op
from ..utils.multimap import OrderedMultiMap
from ..composition import Composition, calculate_mass
from ..composition.structure_composition import monosaccharide_composition
from ..composition.structure_composition import modification_compositions

anomer_map = invert_dict(anomer_map)
superclass_map = invert_dict(superclass_map)


def get_standard_composition(monosaccharide):
    '''Used to get initial composition for a given monosaccharide
    ``Superclass`` and modifications

    Parameters
    ----------
    monosaccharide: Monosaccharide

    Returns
    -------
    Composition
    '''
    base = monosaccharide_composition[monosaccharide.superclass]
    for mod_pos, mod_val in monosaccharide.modifications.items():
        base = base + modification_compositions[mod_val](mod_pos)
    return base


def traverse(monosaccharide, visited=None, apply_fn=ident_op):
    if visited is None:
        visited = set()
    yield apply_fn(monosaccharide)
    visited.add(monosaccharide.id)
    for p, link in monosaccharide.links.items():
        child = link[monosaccharide]
        if child.id in visited:
            continue
        for grandchild in traverse(child, visited=visited, apply_fn=apply_fn):
            yield grandchild


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

    Attributes
    ----------
    anomer: :class:`Anomer` or {'alpha', 'beta', 'uncyclized', 'x', 'missing', None}
        An entry of :class:`~pygly2.structure.constants.Anomer` that corresponds to the linkage type
        of the carbohydrate backbone. Is an entry of a class based on :class:`EnumMeta`
    superclass: :class:`SuperClass` or {'tri', 'tet', 'pen', 'hex', 'hep', 'oct', 'non' 'missing', 'x', None}
        An entry of :class:`~pygly2.structure.constants.SuperClass` that corresponds to the number of
        carbons in the carbohydrate backbone of the monosaccharide. Controls the base composition of the
        instance and the number of positions open to be linked to or modified. Is an entry of a class
        based on :class:`EnumMeta`
    configuration: :class:`Configuration` or {'d', 'l', 'x', 'missing', None}
        An entry of :class:`~pygly2.structure.constants.Configuration` which corresponds to the optical
        stereomer state of the instance. Is an entry of a class based on :class:`EnumMeta`. May possess
        more than one value.
    stem: :class:`Stem` or {"gro", "ery", "rib", "ara", "all", "alt", "glc", "man", "tre", "xyl", "lyx", "gul", "ido", "gal", "tal", "thr", 'x', 'missing', None}
        Corresponds to the bond conformation of the carbohydrate backbone. Is an entry of a class based
        on :class:`EnumMeta`. May possess more than one value.
    ring_start: int or {0, -1, 'x', None}
        The index of the carbon of the carbohydrate backbone that starts a ring. A value of :const:`-1`, :const:`'x'`, or
        :const:`None` corresponds to an unknown start. A value of :const:`0` refers to a linear chain.
    ring_end:  int or {0, -1, 'x', None}
        The index of the carbon of the carbohydrate backbone that ends a ring. A value of :const:`-1`, :const:`'x'`, or
        :const:`None` corresponds to an unknown ends. A value of :const:`0` refers to a linear chain.
    reducing_end: :class:`int`
        The index of the carbon which hosts the reducing end.
    modifications: |OrderedMultiMap|
        The mapping of sites to :class:`~pygly2.structure.constants.Modification` entries. Directly modifies
        the instance's :attr:`Monosaccharide.composition`
    links: |OrderedMultiMap|
        The mapping of sites to :class:`~pygly2.structure.link.Link` entries that refer to other :class:`Monosaccharide` instances
    substituent_links: |OrderedMultiMap|
            The mapping of sites to :class:`~pygly2.structure.link.Link` entries that refer to
            :class:`~pygly2.structure.substituent.Substituent` instances.
    composition: |Composition|
        An instance of :class:`~pygly2.composition.composition.Composition` corresponding to the elemental composition
        of ``self`` and its immediate modifications.

    '''

    def __init__(self, anomer=None, configuration=None, stem=None,
                 superclass=None, ring_start=None, ring_end=None,
                 modifications=None, links=None, substituent_links=None,
                 composition=None, id=None):
        self.anomer = anomer
        self.configuration = configuration
        self.stem = stem
        self.superclass = superclass
        self.ring_start = ring_start
        self.ring_end = ring_end
        self.modifications = OrderedMultiMap(
        ) if modifications is None else modifications
        self.links = OrderedMultiMap() if links is None else links
        self.substituent_links = OrderedMultiMap() if substituent_links\
            is None else substituent_links
        self.id = id or uuid4().int
        if composition is None:
            composition = get_standard_composition(self)
        self.composition = composition
        self._reducing_end = None

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
        if isinstance(value, (list, tuple)):
            self._configuration = tuple(Configuration[v] for v in value)
        else:
            self._configuration = (Configuration[value],)

    @property
    def stem(self):
        return self._stem

    @stem.setter
    def stem(self, value):
        if isinstance(value, (list, tuple)):
            self._stem = tuple(Stem[v] for v in value)
        else:
            if value not in Stem and value is not None:
                raise Exception(value)
            self._stem = (Stem[value],)

    @property
    def superclass(self):
        return self._superclass

    @superclass.setter
    def superclass(self, value):
        if value not in SuperClass and value is not None:
            raise Exception(value)
        self._superclass = SuperClass[value]

    def clone(self, prop_id=False):
        '''
        Copies just this |Monosaccharide| and its |Substituent|s, creating a separate instance
        with the same data. All mutable data structures are duplicated and distinct from the original.

        Does not copy any :attr:`links` as this would cause recursive duplication of the entire |Glycan|
        graph.

        Returns
        -------

        :class:`Monosaccharide`

        '''
        modifications = OrderedMultiMap()
        for k, v in self.modifications.items():
            modifications[k] = Modification[v]
        monosaccharide = Monosaccharide(
            superclass=self.superclass,
            stem=self.stem,
            configuration=self.configuration,
            ring_start=self.ring_start,
            ring_end=self.ring_end,
            modifications=modifications,
            anomer=self.anomer,
            id=self.id if prop_id else None)
        for pos, link in self.substituent_links.items():
            sub = link.to(self)
            dup = sub.clone()
            Link(monosaccharide, dup, link.parent_position, link.child_position,
                 link.parent_loss, link.child_loss)

        return monosaccharide

    # @property
    # def ring_start(self):
    #     return self._ring_start

    # @ring_start.setter
    # def ring_start(self, value):
    #     self._ring_start = value

    # @property
    # def ring_end(self):
    #     return self._ring_end

    # @ring_end.setter
    # def ring_end(self, value):
    #     self._ring_end = value

    @property
    def ring_type(self):
        try:
            diff = self.ring_end - self.ring_start
            if diff == 4:
                return RingType.pyranose
            elif diff == 3:
                return RingType.furanose
            elif diff == 0:
                return RingType.open
        except TypeError:
            return RingType.x


    @property
    def reducing_end(self):
        if self._reducing_end is None:
            for pos, mod in self.modifications.items():
                if mod == ReducingEnd or mod == Modification.aldi:
                    self._reducing_end = pos
                    break
        return self._reducing_end

    @reducing_end.setter
    def reducing_end(self, value):
        if value in {-1, 'x'}:
            raise ValueError("Cannot set the reducing end to be indeterminate")
        if value > self.superclass.value:
            raise IndexError("Index out of bounds")
        pos = self.reducing_end
        if pos is not None:
            self.drop_modification(pos, ReducingEnd)
        self.add_modification(ReducingEnd, value)

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
        max_occupancy: :class:`int`
            The number of objects that may already be bound at a site before it
            is considered unavailable for attachment.

        Returns
        -------
        list : The positions open for binding
        :class:`int` : The number of bound but unknown locations on the backbone.
        '''
        slots = [0] * self.superclass.value
        unknowns = 0
        for pos, obj in chain(self.modifications.items(),
                              self.links.items(),
                              self.substituent_links.items()):
            if pos in {-1, 'x'}:
                unknowns += 1
            else:
                slots[pos - 1] += 1

        reducing_end = self.reducing_end
        if reducing_end is not None:
            slots[reducing_end - 1] -= 1

        open_slots = []
        can_determine_positions = unknowns > 0 or (self.ring_end in {-1, 'x'})

        for i in range(len(slots)):
            if slots[i] <= max_occupancy and (i + 1) != self.ring_end:
                open_slots.append(
                    (i + 1) if not can_determine_positions else -1)

        if self.ring_end in {-1, 'x'}:
            open_slots.pop()

        return open_slots, unknowns

    def is_occupied(self, position):
        '''
        Checks to see if a particular backbone position is occupied by a :class:`~.constants.Modification`,
        :class:`~.substituent.Substituent`, or :class:`~.link.Link` to another :class:`Monosaccharide`.

        Parameters
        ----------
        position: :class:`int`
            The position to check for occupancy. Passing -1 checks for undetermined attachments.

        Returns
        -------
        :class:`int`: The number of occupants at ``position``

        Raises
        ------
        IndexError: When the position is less than 1 or exceeds the limits of the carbohydrate backbone's size.

        '''

        if (position > self.superclass.value) or (position < 1 and position not in {'x', -1, None}):
            raise IndexError("Index out of bounds")
        # The unknown position is always available
        if position in {-1, 'x'}:
            return 0
        n_occupants = len(self.links[position]) +\
            len(self.modifications[position]) +\
            len(self.substituent_links[position])
        if position == self.reducing_end:
            n_occupants -= 1
        return n_occupants

    def add_modification(self, modification, position, max_occupancy=0):
        '''
        Adds a modification instance to :attr:`modifications` at the site given
        by ``position``. This directly modifies :attr:`composition`, consequently
        changing :meth:`mass`

        Parameters
        ----------
        position: |int| or 'x'
            The location to add the :class:`~.constants.Modification` to.
        modification: |str| or Modification
            The modification to add. If passed a |str|, it will be 
            translated into an instance of :class:`~pygly2.structure.constants.Modification`
        max_occupancy: |int|, optional
            The maximum number of items acceptable at `position`. defaults to :const:`1`

        Raises
        ------
        IndexError
            ``position`` exceeds the bounds set by :attr:`superfamily`.
        ValueError
            ``position`` is occupied by more than ``max_occupancy`` elements

        '''
        if self.is_occupied(position) > max_occupancy:
            raise ValueError("Site is already occupied")
        self.composition = self.composition + \
            modification_compositions[modification](position)
        self.modifications[position] = Modification[modification]
        return self

    def drop_modification(self, position, modification):
        if position > self.superclass.value:
            raise IndexError("Index out of bounds")
        self.modifications.pop(position, modification)
        self.composition = self.composition - \
            modification_compositions[modification](position)
        return self

    def add_substituent(self, substituent, position=-1, max_occupancy=0,
                        child_position=1, parent_loss=None, child_loss=None):
        '''
        Adds a :class:`~pygly2.structure.substituent.Substituent` and associated :class:`~pygly2.structure.link.Link`
        to :attr:`substituent_links` at the site given by ``position``. This new substituent is included when calculating mass
        with substituents included

        Parameters
        ----------
        substituent: |str| or |Substituent|
            The substituent to add. If passed a |str|, it will be 
            translated into an instance of |Substituent|
        position: |int| or 'x'
            The location to add the |Substituent| link to :attr:`substituent_links`. Defaults to -1
        child_position: |int|
            The location to add the link to in `substituent`'s :attr:`links`. Defaults to -1. Substituent
            indices are currently not checked.
        max_occupancy: |int|, optional
            The maximum number of items acceptable at ``position``. Defaults to :const:`1`
        parent_loss: |Composition|
            The elemental composition removed from ``self``
        child_loss: |Composition|
            The elemental composition removed from ``substituent``

        Raises
        ------
        IndexError
            ``position`` exceeds the bounds set by :attr:`superfamily`.
        ValueError
            ``position`` is occupied by more than ``max_occupancy`` elements
        '''
        if self.is_occupied(position) > max_occupancy:
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
        if position > self.superclass.value:
            raise IndexError("Index out of bounds")
        if isinstance(substituent, basestring):
            substituent = Substituent(substituent)
        link_obj = None
        for substituent_link in self.substituent_links[position]:
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

    def add_monosaccharide(self, monosaccharide, position=-1, max_occupancy=0,
                           child_position=-1,
                           parent_loss=None, child_loss=None):
        '''
        Adds a |Monosaccharide| and associated |Link| to :attr:`links` at the site given by
        ``position``.

        Parameters
        ----------
        monosaccharide: |Monosaccharide|
            The monosaccharide to add.
        position: |int| or 'x'
            The location to add the |Monosaccharide| link to :attr:`links`. Defaults to -1
        child_position: |int|
            The location to add the link to in `monosaccharide`'s :attr:`links`. Defaults to -1.
        max_occupancy: |int|, optional
            The maximum number of items acceptable at ``position``. Defaults to :const:`1`
        parent_loss: |Composition|
            The elemental composition removed from ``self``
        child_loss: |Composition|
            The elemental composition removed from ``monosaccharide``

        Raises
        ------
        IndexError
            ``position`` exceeds the bounds set by :attr:`superfamily`.
        ValueError
            ``position`` is occupied by more than ``max_occupancy`` elements
        '''
        if parent_loss is None:
            parent_loss = Composition(H=1)
        if child_loss is None:
            child_loss = Composition(H=1, O=1)
        if self.is_occupied(position) > max_occupancy:
            raise ValueError("Parent Site is already occupied")
        if monosaccharide.is_occupied(child_position) > max_occupancy:
            raise ValueError("Child Site is already occupied")
        Link(parent=self, child=monosaccharide,
             parent_position=position, child_position=child_position,
             parent_loss=parent_loss, child_loss=child_loss)
        return self

    def drop_monosaccharide(self, position, refund=True):
        if position > self.superclass.value:
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

    def to_glycoct(self, res_index=None, lin_index=None, complete=True):
        '''
        If ``complete`` is True, returns the Monosaccharide instance formatted as condensed GlycoCT text.
        Otherwise, returns each line of the representation in a partitioned list:
        1. ``list`` of each line in the *RES* section, the monosaccharide residue itself,
        and each of its substituents
        2. ``list`` of each line in the *LIN* section connecting the monosaccharide to each substituent.
        Does *not* include links to other objects from `self.links` to avoid recursively building an
        entire glycan graph.
        3. ``int`` corresponding to the numerical index of the monosaccharide residue to be used to refer
        to it when building monosaccharide to monosaccharide *LIN* entries.

        ``complete = False`` is used by :meth:`~.glycan.Glycan.to_glycoct` to build the full graph in parts.

        Parameters
        ----------
        res_index : function, optional
            A function to yield index numbers for the *RES* section. If not provided, defaults to :func:`~pygly2.utils.make_counter`
        lin_index : function, optional
            A function to yield index numbers for the *LIN* section. If not provided, defaults to :func:`~pygly2.utils.make_counter`
        complete : |bool|, optional
            Format the object completely, returning a self-contained condensed GlycoCT string. Defaults to :const:`True`

        Returns
        -------
        :class:`str`
        '''
        residue_template = "{ix}b:{anomer}{conf_stem}{superclass}-{ring_start}:{ring_end}{modifications}"
        if res_index is None:
            res_index = make_counter()
        if lin_index is None:
            lin_index = make_counter()

        # This index is reused many times
        monosaccharide_index = res_index()

        # Format individual fields
        anomer = anomer_map[self.anomer]
        conf_stem = ''.join("-{0}{1}".format(c.name, s.name)
                            for c, s in zip(self.configuration, self.stem))
        superclass = "-" + superclass_map[self.superclass]
        modifications = '|'.join(
            "{0}:{1}".format(k, v.name) for k, v in self.modifications.items())
        modifications = "|" + modifications if modifications != "" else ""

        # The complete monosaccharide residue line
        residue_str = residue_template.format(ix=monosaccharide_index, anomer=anomer, conf_stem=conf_stem,
                                              superclass=superclass, modifications=modifications,
                                              ring_start=self.ring_start or 'x', ring_end=self.ring_end or 'x')
        res = [residue_str]
        lin = []
        # Construct the substituent lines
        # and their links
        for lin_pos, link_obj in self.substituent_links.items():
            sub = link_obj.to(self)
            sub_index = res_index()
            subst_str = str(sub_index) + sub.to_glycoct()
            res.append(subst_str)
            lin.append(
                link_obj.to_glycoct(lin_index(), monosaccharide_index, sub_index))

        # Completely render the data if `complete`
        if complete:
            buff = StringIO()
            buff.write("RES\n")
            buff.write('\n'.join(res))
            if len(lin) > 0:
                buff.write("\nLIN\n")
                buff.write("\n".join(lin))
            return buff.getvalue()
        else:
            return [res, lin, monosaccharide_index]

    def _flat_equality(self, other, lengths=True):
        '''
        Test for equality of all scalar-ish features that do not
        require recursively comparing links which in turn compare their
        connected units.
        '''
        flat = (self.anomer == other.anomer) and\
            (self.ring_start == other.ring_start) and\
            (self.ring_end == other.ring_end) and\
            (self.superclass == other.superclass) and\
            (self.modifications) == (other.modifications) and\
            (self.configuration) == (other.configuration) and\
            (self.stem) == (other.stem)
        if lengths:
            flat = flat and\
                len(self.links) == len(other.links) and\
                len(self.substituent_links) == len(other.substituent_links) and\
                self.total_composition() == other.total_composition()
        return flat

    def exact_ordering_equality(self, other, substituents=True):
        if self._flat_equality(other):
            if substituents:
                for a_sub, b_sub in zip(self.substituents(), other.substituents()):
                    if a_sub != b_sub:
                        return False
            for a_mod, b_mod in zip(self.modifications.items(), other.modifications.items()):
                if a_mod != b_mod:
                    return False
            for a_child, b_child in zip(self.children(), other.children()):
                if a_child[0] != b_child[0]:
                    return False
                if not a_child[1].exact_ordering_equality(b_child[1], substituents=substituents):
                    return False
            return True
        return False

    def topological_equality(self, other, substituents=True):
        '''
        Performs equality testing between two monosaccharides where
        the exact ordering of child links does not have to match between
        the input |Monosaccharide|s, so long as an exact match of the
        subtrees is found

        Returns
        -------
        |bool|
        '''
        if self._flat_equality(other) and (not substituents or self._match_substituents(other)):
            taken_b = set()
            b_children = list(other.children())
            a_children = list(self.children())
            for a_pos, a_child in a_children:
                matched = False
                for b_pos, b_child in b_children:
                    if (b_pos, b_child.id) in taken_b:
                        continue
                    if a_child.topological_equality(b_child, substituents=substituents):
                        matched = True
                        taken_b.add((b_pos, b_child.id))
                        break
                if not matched and len(a_children) > 0:
                    return False
            if len(taken_b) != len(b_children):
                return False
            return True
        return False

    def _match_substituents(self, other):
        '''
        Helper method for matching substituents in an order-independent
        fashion. Used by :meth:`topological_equality`
        '''
        taken_b = set()
        b_substituents = list(other.substituents())
        cntr = 0
        for a_pos, a_substituent in self.substituents():
            matched = False
            cntr += 1
            for b_pos, b_substituent in b_substituents:
                if b_pos in taken_b:
                    continue
                if b_substituent == a_substituent:
                    matched = True
                    taken_b.add(b_pos)
                    break
            if not matched and cntr > 0:
                return False
        if len(taken_b) != len(b_substituents):
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

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):  # pragma: no cover
        return self.to_glycoct().replace("\n", ' ')

    def mass(self, substituents=True, average=False, charge=0, mass_data=None):
        '''
        Calculates the total mass of ``self``. If ``substituents=True`` it will include
        the masses of its substituents.

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
        :func:`pygly2.composition.composition.calculate_mass`
        '''
        mass = calculate_mass(
            self.composition, average=average, charge=charge, mass_data=mass_data)
        if substituents:
            for link_pos, substituent_link, in self.substituent_links.items():
                mass += substituent_link[self].mass(
                    average=average, charge=charge, mass_data=mass_data)
        return mass

    def total_composition(self):
        '''
        Computes the sum of the composition of ``self`` and each of its linked
        :class:`~pygly2.structure.substituent.Substituent`s

        Returns
        -------
        :class:`~pygly2.composition.Composition`
        '''
        comp = self.composition
        for p, sub in self.substituents():
            comp = comp + sub.total_composition()
        return comp

    def children(self):
        '''
        Returns an iterator over the :class:`Monosaccharide` instancess which are considered
        the descendants of `self`.  Alias for `__iter__`
        '''
        for pos, link in self.links.items():
            if link.is_child(self):
                continue
            yield (pos, link.child)

    def parents(self):
        '''
        Returns an iterator over the :class:`Monosaccharide` instances which are considered
        the ancestors of `self`.
        '''
        for pos, link in self.links.items():
            if link.is_parent(self):
                continue
            yield (pos, link.parent)

    def substituents(self):
        '''
        Returns an iterator over all substituents attached1
        to :obj:`self` by a :class:`~.link.Link` object stored in :attr:`~.substituent_links`
        '''
        for pos, link in self.substituent_links.items():
            yield pos, link.to(self)

    def order(self):
        return len(list(self.children())) + len(self.substituent_links)

    def __iter__(self):
        return self.children()
