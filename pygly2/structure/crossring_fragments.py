import itertools

from . import structure_composition, monosaccharide, constants
from ..composition import Composition
from ..utils.multimap import OrderedMultiMap

RingType = constants.RingType
Monosaccharide = monosaccharide.Monosaccharide
graph_clone = monosaccharide.graph_clone
traverse = monosaccharide.traverse
SuperClass = constants.SuperClass
modification_compositions = structure_composition.modification_compositions


class CrossRingFragment(Monosaccharide):
    '''
    Describes a fragment formed by cleaving across two bonds of the ring of a
    cyclical monosaccharide. Behaves in all respects like an instance of |Monosaccharide|,
    in addition to the following adjustments.

    Attributes
    ----------
    c1, c2: int
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
                 modifications=None, anomer=None,
                 stem=None, configuration=None, id=None, link_cache=None):
        super(CrossRingFragment, self).__init__(
            modifications=modifications, ring_start=None, ring_end=None,
            substituent_links=None, links=None, superclass=None, anomer=anomer,
            stem=stem, configuration=configuration, id=id, composition=composition)
        self.kind = kind
        self.cleave_1 = cleave_1
        self.cleave_2 = cleave_2
        self.contains = contains
        self._link_cache = link_cache or []

    def attach(self):
        '''
        Calls :meth:`pygly2.structure.link.Link.apply` for each |Link| in :attr:`_link_cache`
        '''
        for link in self._link_cache:
            if not link.is_attached():
                link.apply()

    def release(self):
        '''
        Calls :meth:`pygly2.structure.link.Link.break_link` with `refund=True` for each |Link| in :attr:`_link_cache`
        '''
        for link in self._link_cache:
            if link.is_attached():
                link.break_link(refund=True)

    def __repr__(self):
        return "CrossRingFragment({kind}({c1}, {c2}) {contains} {mass})".format(
            kind=self.kind, c1=self.cleave_1, c2=self.cleave_2, contains=self.contains, mass=self.graph_mass())

    def graph_mass(self, average=False, charge=0, mass_data=None):
        '''
        Calcules the total mass of all connected residues.

        See Also
        --------
        pygly2.structure.glycan.Glycan.mass
        '''
        return sum(
            node.mass(average=average, charge=charge, mass_data=mass_data) for node in traverse(self))


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
        c1_fragment.composition = Composition(c1_fragment.composition) - {"O": 1}
        c2_fragment.composition = Composition(c2_fragment.composition) + {"O": 1}

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
        fragment_data['links'].update(i['links'])
        fragment_data['modifications'].update(i['modifications'])
        fragment_data['substituent_links'].update(i['substituent_links'])

    fragment_object = CrossRingFragment(
        fragment_data['composition'], c1, c2, fragment_data['contains'],
        '?', fragment_data['modifications'], stem=residue.stem,
        configuration=residue.configuration, id=residue.id
    )

    composition_shift = Composition()
    for pos, mod in fragment_data["modifications"].items():
        composition_shift = composition_shift + modification_compositions[mod](pos)
    fragment_object.composition = fragment_object.composition + composition_shift

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
        The cyclical carbohydrate to unroll

    Returns
    -------
    list of dict
    '''
    ring_size = residue.superclass.value
    nodes = [None for i in range(ring_size)]
    for i in range(ring_size):
        segment = {
            "index": i + 1,
            "backbone": Composition("COH2"),
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
        c1_segment = linear[:ring_start-1] + c1_segment
        c1_include = list(i + 1 for i in range(ring_start + 1)) + c1_include
    # Alternatively if the ring start is in c2, add the preceding residues to it
    elif ring_start in c2_include and ring_start - 1 not in c2_include and ring_start != 1:
        c2_segment = linear[:ring_start-1] + c2_segment
        c2_include = list(i + 1 for i in range(ring_start + 1)) + c2_include

    # If the ring's end is associated withc2, add all positions proceeding ring end
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
