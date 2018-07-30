import itertools
import collections

from . import monosaccharide, constants
from .link import LinkMaskContext
from glypy.composition import Composition, structure_composition
from glypy.utils.multimap import OrderedMultiMap

RingType = constants.RingType
SuperClass = constants.SuperClass
Configuration = constants.Configuration
Stem = constants.Stem
Anomer = constants.Anomer

UnknownPosition = constants.UnknownPosition

Monosaccharide = monosaccharide.Monosaccharide
ReducedEnd = monosaccharide.ReducedEnd
Modification = monosaccharide.Modification
Link = monosaccharide.Link
graph_clone = monosaccharide.graph_clone
traverse = monosaccharide.traverse
toggle = monosaccharide.toggle
SuperClass = constants.SuperClass
modification_compositions = structure_composition.modification_compositions

chain = itertools.chain
Counter = collections.Counter


class CrossRingPair(object):
    @classmethod
    def can_crossring_fragment(cls, link):
        """Test if the given |Link| object's :attr:`child` supports
        cross ring fragments.

        Parameters
        ----------
        link : Link
            |Link| object whose child will be cross ring cleaved

        Returns
        -------
        bool
        """
        child = link.child
        ring_type = child.ring_type
        return ring_type is not RingType.x and ring_type is not RingType.open

    @classmethod
    def from_link(cls, link):
        """Generate an instance of :class:`CrossRingPair` for each method
        in which the :attr:`child` of `link` can be cross ring cleaved.

        Parameters
        ----------
        link : Link
            |Link| object whose child will be cross ring cleaved

        Yields
        ------
        :class:`CrossRingPair` :
            An instance of :class:`CrossRingPair` for the next viable
            set of ring coordinates
        """
        if not cls.can_crossring_fragment(link):
            return []
        return [cls(link.child, c1, c2) for c1, c2 in
                enumerate_cleavage_pairs(link.child)]

    def __init__(self, residue, c1, c2):
        self.residue = residue
        self.id = residue.id
        self.toggler = None
        self.parent = None
        self.child = None
        self.active = False
        self.cleave_1 = c1
        self.cleave_2 = c2

    def break_link(self, **kwargs):
        """Emulate |Link| interface for passing state with a call to
        :func:`crossring_fragments`.

        Call :func:`crossring_fragments` with :attr:`cleave_1` and :attr:`cleave_2` on
        :attr:`residue` in-place.

        Sets :attr:`parent` to the resulting `X` fragment and :attr:`child` to the resulting
        `A` fragment.

        Sets :attr:`toggler` to an instance of :class:`glypy.structure.link.LinkMaskContext`
        and calls its :meth:`mask` method to hide the links attached to :attr:`residue` after
        they have been duplicated on the resulting fragments.

        Returns
        -------
        :class:`CrossRingFragment` parent: The X fragment from :func:`crossring_fragments`
        :class:`CrossRingFragment` child: The A fragment from :func:`crossring_fragments`
        """
        a, x = crossring_fragments(self.residue, self.cleave_1, self.cleave_2, copy=False)
        # Bind the toggler late to pick up any changes to the residue's links after instantiation
        self.toggler = LinkMaskContext(self.residue)
        self.toggler.mask()
        self.parent = x
        self.child = a
        self.parent._pair = self
        self.child._pair = self
        self.active = True
        return self.parent, self.child

    def apply(self):
        """Emulate |Link| interface for passing state with the removal of
        links to the associated :class:`CrossRingFragment` objects.

        Calls :meth:`CrossRingFragment.release` on :attr:`parent` and :attr:`child`.

        Calls :meth:`LinkMaskContext.unmask`.

        .. warning::
            This method cannot find |Link| objects that are already hidden under another
            :class:`LinkMaskContext`. If using multiple masks on adjacent nodes, make sure
            to rerun :meth:`CrossRingPair.release` again after all links are unmasked.

        """
        self.toggler.unmask()
        self.parent.release()
        self.child.release()
        self.active = False
        self.toggler = None

    def release(self):
        """Convenience function to call :meth:`CrossRingFragment.release` on
        :attr:`parent` and :attr:`child`
        """
        self.parent.release()
        self.child.release()
        self.active = False

    def __repr__(self):  # pragma: no cover
        rep = "<CrossRingPair id={} {},{}>".format(self.id, self.cleave_1, self.cleave_2)
        return rep

    def is_attached(self, deep=False):
        return len(self.child.links) + len(self.parent.links) > 0


class CrossRingFragment(Monosaccharide):
    '''
    Describes a fragment formed by cleaving across two bonds of the ring of a
    cyclical monosaccharide. Behaves in all respects like an instance of |Monosaccharide|,
    in addition to the following adjustments.

    Attributes
    ----------
    cleave_1, cleave_2: int
        Ring sites that were cleaved
    kind: str
        Whether the fragment is reducing (X) or non-reducing (A)
    contains: list
        The list of positions on the monosaccharide included in this fragment
    _link_cache: list
        A list of all |Link| objects connecting this fragment to other |Monosaccharide| objects
        regardless of whether they are attached. Convenient for ease of attaching and
        detatching to a |Glycan|
    '''
    def __init__(self, composition, cleave_1, cleave_2, contains, kind,
                 modifications=None, anomer=Anomer.x,
                 stem=(Stem.x,), configuration=(Configuration.x,), id=None,
                 link_cache=None, source=None):
        _base_composition = composition.clone()
        super(CrossRingFragment, self).__init__(
            modifications=modifications, ring_start=None, ring_end=None,
            substituent_links=None, links=None, superclass=SuperClass.x, anomer=anomer,
            stem=stem, configuration=configuration, id=id, composition=composition, fast=True)
        self._source = source
        self._base_composition = _base_composition
        self.kind = kind
        self.cleave_1 = cleave_1
        self.cleave_2 = cleave_2
        self.contains = contains
        self._link_cache = link_cache or []

    def attach(self):
        '''
        Calls :meth:`glypy.structure.link.Link.apply` for each |Link| in :attr:`_link_cache`
        '''
        for link in self._link_cache:
            if not link.is_attached():
                link.apply()

    def release(self):
        '''
        Calls :meth:`glypy.structure.link.Link.break_link` with `refund=True` for each |Link| in :attr:`_link_cache`
        '''
        for link in self.links.values():
            if link.is_attached():
                link.break_link(refund=True)
        for link in self._link_cache:
            if link.is_attached():
                link.break_link(refund=True)

    def clone(self, prop_id=False, fast=True):
        """Clone this instance, excluding :attr:`links`

        Returns
        -------
        CrossRingFragment

        See Also
        --------
        :meth:`glypy.structure.monosaccharide.Monosaccharide.clone`
        """
        modifications = OrderedMultiMap()
        for k, v in self.modifications.items():
            if isinstance(v, ReducedEnd):
                continue
            modifications[k] = v

        residue = CrossRingFragment(
            self._base_composition.clone(), self.cleave_1, self.cleave_2,
            self.contains, self.kind,
            stem=self._stem,
            configuration=self._configuration,
            modifications=modifications,
            anomer=self._anomer,
            id=self.id if prop_id else None)
        for pos, link in self.substituent_links.items():
            sub = link.to(self)
            dup = sub.clone()
            link.clone(residue, dup)

        return residue

    def __repr__(self):  # pragma: no cover
        return "CrossRingFragment({kind}({c1}, {c2}) {contains} {mass})".format(
            kind=self.kind, c1=self.cleave_1, c2=self.cleave_2, contains=self.contains, mass=self.mass())

    def __eq__(self, other):
        """Test for equality between `self` and `other`.

        Parameters
        ----------
        other

        Returns
        -------
        bool
        """
        if other is None:
            return False
        res = super(CrossRingFragment, self).__eq__(other)

        return res

    def __ne__(self, other):
        return not self == other

    def graph_mass(self, average=False, mass_data=None):
        '''
        Calcules the total mass of all connected residues.

        See Also
        --------
        glypy.structure.glycan.Glycan.mass
        '''
        charge = 0
        return sum(
            node.mass(average=average, charge=charge, mass_data=mass_data) for node in traverse(self))

    def open_attachment_sites(self, max_occupancy=0):
        slots = Counter()
        for i in self.contains:
            slots[i] = 0
        unknowns = 0
        null_positions = (UnknownPosition, 'x')
        for pos, obj in chain(self.modifications.items(),
                              self.links.items(),
                              self.substituent_links.items()):
            if obj == Modification.keto:
                continue
            if pos in null_positions:
                unknowns += 1
            else:
                slots[pos] += 1

        open_slots = []

        for i, count in slots.items():
            if count <= max_occupancy:
                open_slots.append(i)

        return open_slots, unknowns


def crossring_fragments(monosaccharide, c1, c2, attach=True, copy=True):
    '''
    Generate cross-ring fragments from `monosaccharide` cutting bonds at `c1` and `c2`. If
    `attach` is |True|, bonds will be attached between the resulting :class:`CrossRingFragment` and
    other |Monosaccharide| objects. Otherwise they will be created but not `apply`'d. If `copy` is |True|,
    then attached subtrees will be cloned, otherwise they will be shared with the original residue.

    Parameters
    ----------
    monosaccharide: Monosaccharide
        Residue to be cleaved
    c1, c2: int
        Sites to be cleaved at
    attach: bool
        Whether or not |Link| objects between the resulting :class:`CrossRingFragment` objects should
        be made active by `apply` or just saved. Defaults to |True|
    copy: bool
        Whether or not to clone subtrees or share the originals. Defaults to |True|

    Returns
    -------
    tuple of CrossRingFragment
    '''
    ring_type = monosaccharide.ring_type
    if ring_type is RingType.x or ring_type is RingType.open:
        raise TypeError("Cannot cleave an open or unknown carbohydrate backbone")

    c1_segment, c1_include, c2_segment, c2_include = cleave_ring(
        monosaccharide, c1, c2)

    c1_fragment = pack_fragment(c1_segment, c1, c2, c1_include, monosaccharide, attach=attach, copy=copy)
    c2_fragment = pack_fragment(c2_segment, c1, c2, c2_include, monosaccharide, attach=attach, copy=copy)

    ring_start = monosaccharide.ring_start
    ring_end = monosaccharide.ring_end

    a_fragment = x_fragment = None
    if ring_start in c1_fragment.contains:
        x_fragment = c1_fragment
        x_fragment.kind = "X"
        a_fragment = c2_fragment
        a_fragment.kind = "A"
    else:
        x_fragment = c2_fragment
        x_fragment.kind = "X"
        a_fragment = c1_fragment
        a_fragment.kind = "A"

    # If the ring_end is cleaved, then the O component of
    # the ring should go with the X ion
    if ring_end - (ring_start - 1) == c2:
        c1_fragment.composition -= {"O": 1}
        c1_fragment._base_composition -= {"O": 1}
        c2_fragment.composition += {"O": 1}
        c2_fragment._base_composition += {"O": 1}

    return a_fragment, x_fragment


def enumerate_cleavage_pairs(residue):
    '''
    Enumerate all positions `residue` can be cross-ring cleaved at
    '''
    ring_size = (residue.ring_end - residue.ring_start) + 2
    for c1, c2 in itertools.combinations(range(0, ring_size), 2):
        if (c2 - c1 > 1) and ((c2 - c1) + 1 != ring_size):
            yield c1, c2


def pack_fragment(fragment_parts, c1, c2, include, residue, attach=True, copy=True):
    '''
    Given a collection of parts from an unrolled and cleaved ring, assemble them
    into a :clas:`CrossRingFragment` and apply any links or modifications observed
    in the parent residue.

    Parameters
    ----------
    fragment_parts: list
        A list of |dict| objects containing carbohydrate backbone compositions, modification and
        link attachments.
    c1, c2: int
        Sites cleaved at
    include: list
        A list of |int| values referring to the carbons included in this fragment
    residue: Monosaccharide
        The parent |Monosaccharide|. Passed to provide access to :attr:`.Monosaccharide.links`
        for constructing subtrees.

    Returns
    -------
    CrossRingFragment
    '''
    fragment_data = {
        "composition": sum((c['backbone'] for c in fragment_parts), {}),
        "modifications": OrderedMultiMap(),
        "substituent_links": OrderedMultiMap(),
        "links": OrderedMultiMap(),
        "contains": include
    }

    for i in fragment_parts:
        ilink = i['links']
        if ilink:
            fragment_data['links'].update(ilink)

        imodifications = i['modifications']
        if imodifications:
            fragment_data['modifications'].update(imodifications)

        isubstituents = i['substituent_links']
        if isubstituents:
            fragment_data['substituent_links'].update(isubstituents)

    fragment_object = CrossRingFragment(
        fragment_data['composition'], c1, c2, fragment_data['contains'],
        '?', fragment_data['modifications'], stem=residue._stem,
        configuration=residue._configuration, id=residue.id, source=residue
    )

    composition_shift = Composition()
    for pos, mod in fragment_data["modifications"].items():
        composition_shift += modification_compositions[mod](pos)
    fragment_object.composition += composition_shift

    for pos, link in fragment_data['substituent_links'].items():
        subst = link[residue]
        link.clone(fragment_object, subst.clone())

    links = []
    # The id values of all neighboring residues and the parent residue
    edges = {node[residue].id for node in residue.links.values()} | {residue.id}
    for pos, link in fragment_data['links'].items():
        if copy:
            # Copy all structures along the current link excluding all neighbors
            # but the current one to prevent cycles.
            subtree = graph_clone(link[residue], visited=edges - {link[residue].id})
        else:
            subtree = link[residue]
        # Attach the fragment to the subtree and save the link object to the link cache

        if link.is_parent(residue):
            links.append(link.clone(fragment_object, subtree, attach=attach))
        else:
            links.append(link.clone(subtree, fragment_object, attach=attach))
    fragment_object._link_cache = links

    return fragment_object


def unroll_ring(residue):
    '''
    Flattens `residue` into a list of |dict| objects containing
    information about the `i`th position on the carbohydrate backbone

    Parameters
    ----------
    residue: Monosaccharide
        The carbohydrate to unroll

    Returns
    -------
    list of dict
    '''
    ring_size = residue.superclass.value
    nodes = [None for i in range(ring_size)]
    for i in range(ring_size):
        segment = {
            "index": i + 1,
            "backbone": Composition({'H': 2, 'C': 1, 'O': 1}),
            "modifications": OrderedMultiMap(),
            "substituent_links": OrderedMultiMap(),
            "links": OrderedMultiMap()
        }
        for mod in (residue.modifications[i + 1]):
            segment['modifications'][i + 1] = mod

        for sub_link in residue.substituent_links[i + 1]:
            segment['substituent_links'][i + 1] = sub_link

        for glyco_link in residue.links[i + 1]:
            segment['links'][i + 1] = glyco_link

        nodes[i] = segment

    return nodes


def cleave_ring(residue, c1, c2):
    '''
    Given a |Monosaccharide| `residue` and two ring cleavage sites, `(c1, c2)`, calculate
    the linear structures produced by cleaving the ring of `residue`, including any non-ring
    carbons with the appropriate fragment. Calls :func:`slice_ring` and :func:`unroll_ring`

    Parameters
    ----------
    residue: Monosaccharide
        The residue to be cleaved
    c1, c2: int
        The cleavage sites

    Returns
    -------
    c1_segment: list of dicts of backbone information
    c1_include: list of ints
    c2_segment: list of dicts of backbone information
    c2_include: list of ints
    '''

    linear = unroll_ring(residue)
    ring_start = residue.ring_start
    ring_end = residue.ring_end

    # Slice the ring region of the linearized carbohydrate
    c1_segment, c1_include, c2_segment, c2_include = slice_ring(
        linear[ring_start - 1: ring_end], c1, c2)

    # If the ring's start is associated with c1, add all positions prior to ring start
    # to c1.
    if ring_start in c1_include and ring_start - 1 not in c1_include and ring_start != 1:
        c1_segment = linear[:ring_start - 1] + c1_segment
        c1_include = list(i + 1 for i in range(ring_start + 1)) + c1_include
    # Alternatively if the ring start is in c2, add the preceding residues to it
    elif ring_start in c2_include and ring_start - 1 not in c2_include and ring_start != 1:
        c2_segment = linear[:ring_start - 1] + c2_segment
        c2_include = list(i + 1 for i in range(ring_start + 1)) + c2_include

    # If the ring's end is associated with c2, add all positions proceeding ring end
    # to c2
    if ring_end in c2_include:
        c2_segment = c2_segment + linear[ring_end:]
        c2_include = c2_include + \
            list(i + 1 for i in range(ring_end, len(linear)))
    # else they must be in c1
    else:
        c1_segment = c1_segment + linear[ring_end:]
        c1_include = c1_include + \
            list(i + 1 for i in range(ring_end, len(linear)))

    return c1_segment, c1_include, c2_segment, c2_include


def slice_ring(linear, c1, c2):
    '''
    Given a linear carbohydrate and two cut sites, extract the two regions
    produced by cutting the sequence like a ring buffer.
    '''
    # The first fragment is defined as the region between the first and
    # second cut.
    c1_segment = linear[c1:c2]
    c1_include = [c['index'] for c in c1_segment]
    # Add the positions after the second cut up to the end of the ring
    # and the positions from the start of the ring to the first cut to the
    # second fragment.
    c2_segment = [linear[i]
                  for i in list(range(c2, len(linear))) + list(range(c1))
                  ]
    c2_include = [c['index'] for c in c2_segment]
    return c1_segment, c1_include, c2_segment, c2_include
