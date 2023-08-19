import itertools
import re
from dataclasses import dataclass, field
from typing import DefaultDict, Dict, FrozenSet, List, Set, Optional, Deque, TYPE_CHECKING, Tuple, Union

from glypy.composition import Composition

from glypy.utils import tree as treep

if TYPE_CHECKING:
    from glypy.structure import Glycan, Link, Monosaccharide
    from glypy.structure.glycan_composition import HashableGlycanComposition


_fragment_shift = {
    "B": Composition(O=1, H=2),
    "Y": Composition(),
    "C": Composition(),
    "Z": Composition(H=2, O=1),
}


def _link_ids_splitter(fragment, link_ids, kind):  # pragma: no cover
    ion_types = re.findall(r"(\d+,\d+)?(\S)", kind)
    links_broken = link_ids

    pairings = zip(ion_types, links_broken)

    fragment.link_ids = {
        link_id: ion_type for ion_type,
        link_id in pairings if ion_type[0] == ""}
    fragment.crossring_cleavages = {
        node_id: ion_type for ion_type,
        node_id in pairings if ion_type[0] != ""}


latex_symbol_map = {
    "a": r"\alpha",
    "b": r"\beta",
    "c": r"\sigma",
    "d": r"\delta",
    "e": r"\epsilon",
    "f": r"\zeta",
    "g": r"\gamma",
    "h": r"\eta",
    "i": r"\iota",
    "j": r"\lambda",
    "k": r"\kappa",
    "l": r"\lambda",
    "m": r"\mu",
    "n": r"\nu",
    "o": r"\omicron",
    "p": r"\pi",
    "q": r"\sigma",
    "r": r"\rho",
    "s": r"\upsilon",
    "t": r"\tau",
    "u": r"\upsilon",
    "v": r"\omega",
}


class GlycanFragment(object):
    '''
    A simple container for a fragment ion, produced by :meth:`Glycan.fragments`

    Attributes
    ----------
    kind: |str|
        One of A, B, C, X, Y, or Z for each link broken or ring cleaved

    link_ids: |dict| of |int| -> |tuple|
        The :attr:`id` value of each link cleaved to the corresponding cleavage type

    included_nodes: |list| of |int|
        The :attr:`id` value of each |Monosaccharide| contained in the fragment

    mass: |float|
        The mass of the fragment

    name: |str|
        The fragment name under the branching nomenclature

    crossring_cleavages: |dict| of |int| -> |tuple|
        The :attr:`id` value of each link cleaved to the corresponding cleavage type, including ring coordinates

    composition: |Composition|
        The elemental composition of the fragment

    See Also
    --------
    :meth:`Glycan.fragments`
    '''
    __slots__ = [
        "mass",
        "kind",
        "included_nodes",
        "link_ids",
        "name",
        "crossring_cleavages",
        "composition",
    ]

    mass: float
    kind: str
    included_nodes: Set[int]
    link_ids: List[int]
    name: str
    crossring_cleavages: Dict
    composition: Optional[Composition]

    def __init__(self, kind, link_ids, included_nodes, mass,
                 name=None, crossring_cleavages=None, composition=None):
        self.mass = mass
        self.kind = kind
        self.link_ids = link_ids
        self.included_nodes = included_nodes
        self.crossring_cleavages = crossring_cleavages
        self.name = name
        self.composition = composition

    def copy(self):
        """Create a deep copy of the fragment.

        Returns
        -------
        GlycanFragment
        """
        return self.__class__(
            self.kind, self.link_ids.copy(), self.included_nodes.copy(), self.mass,
            self.name, self.crossring_cleavages.copy(), self.composition.copy())

    def clone(self):
        """Create a deep copy of the fragment.

        Returns
        -------
        GlycanFragment
        """
        return self.copy()

    def is_reducing(self):
        """Is this fragment from the reducing end

        Returns
        -------
        |bool|
        """
        return bool(set(self.kind) & set("XYZ"))

    def is_non_reducing(self):
        """Is this fragment from the non-reducing end

        Returns
        -------
        |bool|
        """
        return bool(set(self.kind) & set("ABC"))

    def is_internal(self):
        """Is the fragment composed of both reducing and non-reducing bond
        cleavage events.

        Returns
        -------
        bool
        """
        return bool(self.is_reducing() and self.is_non_reducing())

    @property
    def fname(self):
        """A `LaTeX` formatted version of the fragment's name

        Returns
        -------
        str
        """
        buff = []
        for c in self.name:
            if c in latex_symbol_map:
                buff.append("$_{}$".format(latex_symbol_map[c]))
            else:
                buff.append(c)
        return ''.join(buff)

    def _make_tex_name(self):
        buff = []
        for c in self.name:
            if c in latex_symbol_map:
                buff.append("_{}".format(latex_symbol_map[c]))
            else:
                buff.append(c)
        return ''.join(buff)

    @property
    def series(self):
        """An alias for :attr:`kind`

        Returns
        -------
        str
        """
        return self.kind

    def __getstate__(self):  # pragma: no cover
        d = {}
        for a in self.__slots__:
            d[a] = getattr(self, a)
        return d

    def __setstate__(self, state):  # pragma: no cover
        for a, v in state.items():
            setattr(self, a, v)

    def __hash__(self):
        return hash((self.name, int(self.mass)))

    def __eq__(self, other):  # pragma: no cover
        for f in self.__slots__:
            if getattr(self, f) != getattr(other, f, NotImplemented):
                return False
        return True

    def __ne__(self, other):  # pragma: no cover
        return not self == other

    def __repr__(self):  # pragma: no cover
        rep = "<GlycanFragment "
        for f in self.__slots__:
            rep += " {}={}".format(f, getattr(self, f))
        rep += ">"
        return rep

    @property
    def __dict__(self):
        return self.__getstate__()

    @property
    def break_count(self):
        """The number of cleavage events to create this fragment.

        Returns
        -------
        int
        """
        return len(self.link_ids) + len(self.crossring_cleavages)

    @property
    def residues_contained(self):
        """The number of :class:`~.Monosaccharide` nodes included in the fragmented sub-tree.

        Returns
        -------
        int
        """
        return len(self.included_nodes)

    @classmethod
    def to_glycan_compositions(cls, glycan: 'Glycan', fragments: List['GlycanFragment'],
                               by_series: bool = True) -> Union[
                                   DefaultDict['HashableGlycanComposition', List['GlycanFragment']],
                                   DefaultDict[str, Dict['HashableGlycanComposition',
                                                         List['GlycanFragment']]]
                                ]:
        """
        From a list of :class:`GlycanFragment` instances, build
        :class:`~glypy.structure.glycan_composition.HashableGlycanComposition`
        instances corresponding to those fragments, and return a mapping relating
        them.

        Parameters
        ----------
        glycan : :class:`~.glypy.structure.glycan.Glycan`
            The glycan the fragments came from.
        fragments : :class:`list` of :class:`~
        """
        from glypy.structure.glycan_composition import (
            FrozenMonosaccharideResidue, HashableGlycanComposition)
        index_to_residue = {
            node.id: FrozenMonosaccharideResidue.from_monosaccharide(
                node, False, False, False, False
            )
            for node in glycan
        }

        compositions: DefaultDict[HashableGlycanComposition, List['GlycanFragment']] = DefaultDict(list)

        for frag in fragments:
            gc = HashableGlycanComposition()
            for node_id in frag.included_nodes:
                gc[index_to_residue[node_id]] += 1
            compositions[gc].append(frag)
        if not by_series:
            return compositions

        _shift_cache = {}
        results: DefaultDict[str, Dict[HashableGlycanComposition, List[GlycanFragment]]] = DefaultDict(dict)
        for gc, frags in compositions.items():
            frags.sort(key=lambda x: x.kind)
            for key, subset in itertools.groupby(frags, lambda x: x.kind):
                if key not in _shift_cache:
                    comp_shift = Composition.sum([_fragment_shift[k] for k in key])
                    _shift_cache[key] = comp_shift
                else:
                    comp_shift = _shift_cache[key]
                tmp = gc.clone()
                tmp.composition_offset = tmp.composition_offset - comp_shift
                results[key][tmp] = list(subset)
        return results


Fragment = GlycanFragment


class GlycanSubstructure(object):
    """Represent a subgraph of a :class:`~.Glycan` produced
    from another :class:`~.Glycan` by breaking linkages.

    Supports the :func:`~.root` and :func:`~.tree` protocols.

    Attributes
    ----------
    tree : :class:`~.Glycan`
        The substructure
    include_nodes : set
        The :attr:`~.Monosaccharide.id` of all included monosaccharides
    link_ids : list
        The :attr:`~.Link.id` of all linkages broken to form this substructure
    parent_breaks : dict
        Mapping from :attr:`~.Link.id` to :attr:`~.Monosaccharide.id` if that monosaccharide
        was the :attr:`~.Link.parent` of that linkage and that linkage was broken.
    child_breaks : dict
        Mapping from :attr:`~.Link.id` to :attr:`~.Monosaccharide.id` if that monosaccharide
        was the :attr:`~.Link.child` of that linkage and that linkage was broken.
    crossring_cleavages : dict
        Mapping from :attr:`~.Monosaccharide.id` to :class:`str` denoting the ring coordinate pair that
        was cleaved to produce any included crossring cleavage events.

    """

    tree: 'Glycan'
    include_nodes: Set[int]
    link_ids: List[int]
    parent_breaks: Dict[int, int]
    child_breaks: Dict[int, int]
    crossring_cleavages: Dict

    def __init__(self, tree, include_nodes, link_ids,
                 parent_breaks, child_breaks, crossring_cleavages=None):
        self.tree = tree
        self.include_nodes = include_nodes
        self.link_ids = link_ids
        self.parent_breaks = parent_breaks
        self.child_breaks = child_breaks
        self.crossring_cleavages = crossring_cleavages or {}

    def __hash__(self):
        return hash(frozenset(self.include_nodes))

    def __eq__(self, other):
        return treep(self) == treep(other)

    def __ne__(self, other):
        return not self == other

    def contains_reducing_end(self):
        """Whether or not the sub-tree includes a :class:`ReducedEnd`

        Returns
        -------
        bool
        """
        return self.tree.root.reducing_end is not None

    def to_fragments(self, kind="BY", average=False, charge=None, mass_data=None, include_composition=True,
                     traversal_method='dfs'):
        """Transform an instance of :class:`GlycanSubstructure` into every combination of
        :class:`GlycanFragment` allowed under `kind`.

        Parameters
        ----------
        kind : Iterable, optional
            The types of fragments to emit. Defaults to "BY"
        average : bool, optional
            Calculate masses with average isotopic composition
        charge : int, optional
            Calculate `m/z` instead of neutral mass, with `z = charge`
        mass_data : dict, optional
            If mass_data is None, standard NIST mass and isotopic abundance data are used. Otherwise the
            contents of mass_data are assumed to contain elemental mass and isotopic abundance information.
            Defaults to :const:`None`.
        include_composition: bool, optional
            Whether or not to populate the `composition` attribute of the fragment. Defaults to :const:`True`

        Yields
        -------
        :class:`GlycanFragment`
        """
        parent_type = set("YZ") & set(kind)
        child_type = set("BC") & set(kind)
        crossring_type = set("AX") & set(kind)
        parent_shifts = [(lid, parent_type)
                         for lid in self.parent_breaks.keys()]
        child_shifts = [(lid, child_type)
                        for lid in self.child_breaks.keys()]

        shift_masses = {k: _fragment_shift[k].calc_mass(average=average,
                                                        charge=charge,
                                                        mass_data=mass_data)
                        for k in parent_type | child_type}
        crossring_contained = {
            kind for pos,
            kind in self.crossring_cleavages.values()}
        # Do not proceed if crossring fragments are included and
        # were not requested.
        if len(crossring_contained & crossring_type) == 0 and len(
                crossring_contained) != 0:

            return

        all_shifts = parent_shifts + child_shifts
        all_link_ids = [i for i, t in all_shifts]
        frag_types = [t for i, t in all_shifts]
        base_mass = self.tree.mass(
            average=average,
            charge=charge,
            mass_data=mass_data, method=traversal_method)
        if include_composition:
            base_composition = self.tree.total_composition(method=traversal_method)
        # product of splat of empty list is a list of the empty list. So a fragment with
        # no glycosidic cleavages still enters this outer loop, letting only crossring-cleavage
        # GlycanSubstructure instances through without issue
        for shift_set in itertools.product(*frag_types):
            mass_offset = 0.0
            if include_composition:
                composition_offset = Composition()
            link_ids = {}
            # The type of fragment being produced, expressed a collection of ABCXYZs
            kind = [] + [''.join(cr_kind)
                         for cr_kind in self.crossring_cleavages.values()]
            i = 0
            shift_set = list(shift_set)
            for shift in shift_set:
                link_id = all_link_ids[i]
                shift = shift[0]
                mass_offset -= shift_masses[shift]
                if include_composition:
                    composition_offset -= _fragment_shift[shift]
                link_ids[link_id] = ("", shift)
                kind.append(shift)
                i += 1

            if include_composition:
                fragment_composition = base_composition + composition_offset
            else:
                fragment_composition = None

            yield GlycanFragment(kind=''.join(kind), link_ids=link_ids, included_nodes=self.include_nodes,
                                 mass=base_mass + mass_offset, name=None,
                                 crossring_cleavages=self.crossring_cleavages,
                                 composition=fragment_composition)

    def __repr__(self):  # pragma: no cover
        rep = ("<GlycanSubstructure include_nodes={} link_ids={} parent_breaks={}"
               " child_breaks={} crossring_cleavages={}>\n{}").format(
                   self.include_nodes, self.link_ids, self.parent_breaks,
                   self.child_breaks, self.crossring_cleavages, self.tree)
        return rep

    def __root__(self):  # pragma: no cover
        return self.tree.root

    def __tree__(self):  # pragma: no cover
        return self.tree


Subtree = GlycanSubstructure


def flatten(x: List) -> List:
    return [
        z for y in x
        for z in ([y] if not isinstance(y, list) else flatten(y))
    ]


@dataclass
class DepthFirstLinkTraversal:
    path: List["Link"] = field(default_factory=list)
    children: List["DepthFirstLinkTraversal"] = field(default_factory=list)
    visited: Set[int] = field(default_factory=set)
    size: int = 0

    @classmethod
    def build_from(
        cls,
        root: "Monosaccharide",
        starting_link: Optional["Link"] = None,
        traverse_children: bool = True,
        max_size: int = -1,
        visited: Optional[Set[int]] = None
    ):
        self = cls(visited=visited if visited else set())
        node_stack: Deque["Monosaccharide"] = Deque()
        if starting_link is not None:
            node_stack.append(starting_link[root])
            self.path.append(starting_link)
            self.visited.add(starting_link.id)
        else:
            node_stack.append(root)
        while len(node_stack) > 0:
            if max_size >= 0 and len(self.path) >= max_size:
                break
            node = node_stack.pop()
            children = [
                c[1]
                for c in (
                    node.children(
                        True) if traverse_children else node.links.items()
                )
                if c[1].id not in self.visited
            ]
            if len(children) == 1:
                child = children[0]
                self.path.append(child)
                self.visited.add(child.id)
                node_stack.appendleft(child[node])
            elif len(children) > 1:
                for child in children:
                    self.visited.add(child.id)
                    self.children.append(
                        cls.build_from(
                            node,
                            child,
                            traverse_children=traverse_children,
                            max_size=max_size - len(self.path),
                            visited=self.visited.copy()
                        )
                    )
        return self

    def __iter__(self):
        for link in self.path:
            yield [link]
        for links in itertools.product(*[[None] + list(c) for c in self.children]):
            links = [li for li in links if li is not None]
            if not links:
                continue
            yield links

    @classmethod
    def generate_y_fragments(cls, glycan: 'Glycan',
                             include_composition: bool=False,
                             traversal_method: str="index", **kwargs):
        fragments = list(
            itertools.chain.from_iterable(
                (y_fragments_from_links(
                    links,
                    include_composition=include_composition,
                    traversal_method=traversal_method, **kwargs
                ) for links in map(flatten, cls.build_from(glycan.root)))
            )
        )
        link_index = {link.id: link for _, link in glycan.iterlinks()}
        for frag in fragments:
            frag.name = glycan.name_fragment(frag, link_index=link_index)
        return fragments

    @classmethod
    def _generate_b_fragments(cls, glycan: "Glycan", include_composition: bool=False,
                             traversal_method: str="index", max_size: int = 4, **kwargs):
        """Incomplete"""
        seen = set()
        bond_set_iterator = (
            map(lambda x: flatten(x) if isinstance(x, list) else [x],
                itertools.chain.from_iterable(
                    itertools.chain.from_iterable(cls.build_from(
                        leaf,
                        traverse_children=False,
                        max_size=max_size
                    ) for leaf in glycan.leaves()),
                )
            )
        )

        def c(x):
            return b_fragments_from_links(
                x, include_composition=include_composition,
                traversal_method=traversal_method,
                max_size=max_size,
                seen=seen,
            )

        fragments = itertools.chain.from_iterable(map(c, bond_set_iterator))
        return fragments
        # return bond_set_iterator

        # return (itertools.chain.from_iterable(
        #     map(
        #         lambda x: b_fragments_from_links(
        #             x, include_composition=include_composition,
        #             traversal_method=traversal_method,
        #             max_size=max_size,
        #             seen=seen,
        #         ),
        #         map(
        #             flatten,
        #             itertools.chain.from_iterable(cls.build_from(
        #                 leaf,
        #                 traverse_children=False,
        #                 max_size=max_size
        #             ) for leaf in glycan.leaves()),
        #         ),
        #     )
        # ))


def b_fragments_from_links(links_to_break: List['Link'],
                           max_size: int = 4,
                           seen: Optional[Set[FrozenSet[int]]] = None,
                           **kwargs):
    from glypy.structure import Glycan
    subtrees = []
    if seen is None:
        seen = set()
    for link in links_to_break:
        parent, child = link.break_link(refund=True)
        child_tree = Glycan(child, index_method=None)
        child_tree._build_node_index()
        if len(child_tree) > max_size:
            continue
        subtrees.append(child_tree)

    unique_subtrees = []
    for subtree in subtrees:
        ids = frozenset([n.id for n in subtree])
        if ids in seen:
            continue
        for uids, _unique in unique_subtrees:
            if ids == uids:
                break
        else:
            unique_subtrees.append((ids, subtree))

    for ids, subtree in unique_subtrees:
        subtree = subtree.reroot(index_method=None)
        include_nodes = ids
        seen.add(ids)

        link_ids = [link.id for link in links_to_break
                    if link.parent.id in include_nodes or
                    link.child.id in include_nodes]

        parent_break_ids = {link.id: link.parent.id for link in links_to_break
                            if link.parent.id in include_nodes}

        child_break_ids = {link.id: link.child.id for link in links_to_break
                           if link.child.id in include_nodes}

        yield from Subtree(
            subtree,
            include_nodes,
            link_ids,
            parent_break_ids,
            child_break_ids).to_fragments('B', **kwargs)

    for link in links_to_break:
        link.apply()


def y_fragments_from_links(links_to_break: List['Link'], **kwargs):
    from glypy.structure import Glycan
    subtrees = []
    for link in links_to_break:
        parent, _child = link.break_link(refund=True)
        parent_tree = Glycan(parent, index_method=None)
        subtrees.append(parent_tree)

    unique_subtrees = []
    for subtree in subtrees:
        ids = {n.id for n in subtree}
        for uids, _unique in unique_subtrees:
            if ids == uids:
                break
        else:
            unique_subtrees.append((ids, subtree))

    assert len(unique_subtrees) == 1

    for ids, subtree in unique_subtrees:
        subtree = subtree.reroot(index_method=None)
        include_nodes = ids

        link_ids = [link.id for link in links_to_break
                    if link.parent.id in include_nodes or
                    link.child.id in include_nodes]

        parent_break_ids = {link.id: link.parent.id for link in links_to_break
                            if link.parent.id in include_nodes}

        child_break_ids = {link.id: link.child.id for link in links_to_break
                           if link.child.id in include_nodes}

        yield from Subtree(
            subtree,
            include_nodes,
            link_ids,
            parent_break_ids,
            child_break_ids).to_fragments('Y', **kwargs)

    for link in links_to_break:
        link.apply()


def y_fragments_to_glycan_compositions(glycan: 'Glycan',
                                       fragments: List[GlycanFragment],
                                       composition_offset: Optional[Composition] = None) -> DefaultDict[
                                           'HashableGlycanComposition',
                                           List[GlycanFragment]
                                        ]:
    from glypy.structure.glycan_composition import (
        FrozenMonosaccharideResidue, HashableGlycanComposition)
    index_to_residue = {
        node.id: FrozenMonosaccharideResidue.from_monosaccharide(
            node, False, False, False, False
        )
        for node in glycan
    }

    compositions = DefaultDict(list)

    for frag in fragments:
        gc = HashableGlycanComposition()
        for node_id in frag.included_nodes:
            gc[index_to_residue[node_id]] += 1
        if composition_offset is not None:
            gc.set_composition_offset(composition_offset)
        compositions[gc].append(frag)

    return compositions


class LabileLeafMask:
    glycan: 'Glycan'
    monosaccharides: List['Monosaccharide']
    links: List['Link']

    def __init__(self, glycan: 'Glycan', monosaccharides: List['Monosaccharide']):
        from glypy.algorithms.similarity import commutative_similarity
        self.glycan = glycan
        self.monosaccharides = monosaccharides
        self.links = []
        for leaf in self.glycan.leaves():
            _, link = next(iter(leaf.links.items()))
            for ref in self.monosaccharides:
                if commutative_similarity(ref, leaf):
                    self.links.append(link)
                    break

    def __enter__(self):
        for link in self.links:
            link.break_link(refund=True)
        return self.links

    def __exit__(self, exc_type, exc_val, exc_tb):
        for link in self.links:
            link.apply()

