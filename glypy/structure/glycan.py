import operator
import logging
import itertools
from functools import partial
from collections import deque, defaultdict, Callable
from uuid import uuid4

from .base import SaccharideBase
from .monosaccharide import Monosaccharide, graph_clone, toggle as residue_toggle, depth
from .crossring_fragments import crossring_fragments, CrossRingPair
from .fragment import Subtree
from ..utils import make_counter, identity, StringIO, chrinc
from ..composition import Composition

methodcaller = operator.methodcaller
logger = logging.getLogger("Glycan")

_fragment_direction = {
    "A": -1,
    "B": -1,
    "C": -1,
    "X": 1,
    "Y": 1,
    "Z": 1
}

MAIN_BRANCH_SYM = '-'


def fragment_to_substructure(fragment, tree):
    """Extract the substructure of `tree` which is contained in `fragment`


    >>> from glypy import glycans as glycan_factory
    >>> from glypy.structure import glycan
    >>> n_linked_core = glycan_factory["N-Linked Core"]
    >>> frag = n_linked_core.fragments().next()
    >>> frag
    <Fragment  mass=221.089937203 kind=Y included_nodes=set([1]) link_ids={1: ('', 'Y')} name=Y1 crossring_cleavages={} score=0.0>
    >>> glycan.fragment_to_substructure(frag, n_linked_core)
    RES
    1b:b-dglc-HEX-1:5
    2s:n-acetyl
    LIN
    1:1d(2+1)2n
    <BLANKLINE>
    >>>

    Parameters
    ----------
    fragment: Fragment
        The :class:`Fragment` to extract substructure for.
    tree: Glycan
        The |Glycan| to extract substructure from.

    Returns
    -------
    Glycan:
        The |Glycan| substructure defined by the nodes contained in `fragment` as
        found in `tree`
    """
    break_targets = fragment.link_ids
    crossring_targets = fragment.crossring_cleavages

    # All operations will be done on a copy of the tree of interest
    tree = tree.clone()
    crossring_targets_nodes = []
    break_targets_nodes = []
    # A point of reference known to be inside the fragment tree
    anchor = None
    for pos, link in tree.iterlinks():
        # If the current link's child was cross-ring cleaved,
        # then we must apply that cleavage and find an anchor
        if link.child.id in crossring_targets:
            ion_type = crossring_targets[link.child.id]
            c1, c2 = map(int, ion_type[0].split(","))
            target = tree.get(link.child.id)
            crossring_targets_nodes.append((target, c1, c2))
        # If this link was cleaved, break it and find an anchor
        if link.id in break_targets:
            break_targets_nodes.append(link)

    for target, c1, c2 in crossring_targets_nodes:
        a_frag, x_frag = crossring_fragments(
            target, c1, c2, attach=True, copy=False)
        residue_toggle(target).next()

        if crossring_targets[target.id][1] == "A":
            anchor = a_frag
        else:
            anchor = x_frag

    for link in break_targets_nodes:
        parent, child = link.break_link(refund=True)
        if parent.id in fragment.included_nodes:
            anchor = parent
        elif child.id in fragment.included_nodes:
            anchor = child

    # Build a new tree from the anchor
    substructure = Glycan(root=anchor, index_method=None).reroot()

    return substructure


class Glycan(SaccharideBase):

    '''
    Represents a full graph of connected |Monosaccharide| objects and their connecting bonds.

    Attributes
    ----------
    root: |Monosaccharide|
        The first monosaccharide unit of the glycan, and the reducing end if present.
    index: |list|
        A list of the |Monosaccharide| instances in `self` in the order they are encountered
        by traversal by `traversal_methods[index_method]`
    link_index: |list|
        A list of the |Link| connecting the |Monosaccharide| instances in `self` in the order they
        are encountered by traversal by `traversal_methods[index_method]`
    reducing_end: |ReducedEnd| or |None|
        The reducing end on :attr:`root`.
    branch_lengths: |dict|
        A dictionary mapping branch symbols to their lengths
    '''

    @classmethod
    def subtree_from(cls, tree, node, visited=None):
        if isinstance(node, int):
            node = tree[node]
        visited = {
            node.id for p,
            node in node.parents()} if visited is None else visited
        subtree = cls(root=node, index_method=None).clone(
            index_method=None, visited=visited)
        return subtree

    traversal_methods = {}

    def __init__(self, root=None, index_method='dfs'):
        '''
        Constructs a new Glycan from the collection of connected |Monosaccharide| objects
        rooted at `root`.

        If index_method is not |None|, the graph is indexed by the default search method
        given by `traversal_methods[index_method]`
        '''
        if root is None:
            root = Monosaccharide()
        self.root = root
        self.index = []
        self.link_index = []
        self.branch_lengths = {}
        if index_method is not None:
            self.reindex(index_method)

    def reindex(self, method='dfs'):
        '''
        Traverse the graph using the function specified by ``method``. The order of
        traversal defines the new :attr:`id` value for each |Monosaccharide|
        and |Link|.

        The order of traversal also defines the ordering of the |Monosaccharide|
        in :attr:`index` and |Link| in :attr:`link_index`.

        '''
        self.deindex()
        traversal = self._get_traversal_method(method)
        index = []
        visited = set()
        i = 1
        for node in traversal():
            addr = id(node)
            if addr in visited:
                continue
            visited.add(addr)
            index.append(node)
        for node in index:
            if node.id < 0:
                node.id = i
                i += 1

        link_index = []
        for pos, link in self.iterlinks(method=method):
            link_index.append(link)

        i = 1
        for link in link_index:
            link.id = i
            i += 1

        self.index = index
        self.link_index = link_index

        self.label_branches()

        return self

    def deindex(self):
        '''
        When combining two Glycan structures, very often their component ids will
        overlap, making it impossible to differentiate between a cycle and the new
        graph. This function mangles all of the node and link ids so that they are
        distinct from the pre-existing nodes.
        '''
        if self.index is not None and len(self.index) > 0:
            base = uuid4().int
            for node in self.index:
                node.id += base
                node.id *= -1
            for link in self.link_index:
                link.id += base
                link.id *= -1
        return self

    def reroot(self, index_method='dfs'):
        '''
        Set :attr:`root` to the node with the lowest :attr:`id`
        '''
        self.root = sorted(self, key=operator.attrgetter('id'))[0]
        if index_method is not None:
            self.reindex(index_method)
        return self

    def __getitem__(self, ix):
        '''
        Alias for :attr:`index.__getitem__`
        '''
        if self.index is not None:
            return self.index[ix]
        else:
            raise IndexError(
                "Tried to access the index of an unindexed Glycan.")

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, state):
        self.__dict__ = state
        if "_root" in state:
            self.root = self._root

    def __root__(self):
        return self.root

    def get(self, ix):
        for node in self:
            if node.id == ix:
                return node
        raise IndexError(
            "Could not find a node with the given id {}".format(ix))

    def get_link(self, ix):
        for pos, link in self.iterlinks():
            if link.id == ix or link.label == ix:
                return link
        raise IndexError(
            "Could not find a link with the given id or label {}".format(ix))

    @property
    def reducing_end(self):
        '''
        An alias for :attr:`Monosaccharide.reducing_end` for :attr:`root`
        '''
        return self.root.reducing_end

    def set_reducing_end(self, value):
        '''
        Sets the reducing end type, and configures the root residue appropriately.

        If the reducing_end is not |None|, then the following state changes are made to the root:

        .. code-block:: python

            self.root.ring_start = 0
            self.root.ring_end = 0
            self.root.anomer = "uncyclized"

        Else, the correct state is unknown:

        .. code-block:: python

            self.root.ring_start = None
            self.root.ring_end = None
            self.root.anomer = None

        '''
        self.root.reducing_end = value
        if self.reducing_end is not None:
            self.root.ring_start = 0
            self.root.ring_end = 0
            self.root.anomer = "uncyclized"
        else:
            self.root.ring_start = None
            self.root.ring_end = None
            self.root.anomer = None

    @reducing_end.setter
    def reducing_end(self, value):
        self.set_reducing_end(value)

    def depth_first_traversal(
            self, from_node=None, apply_fn=identity, visited=None):
        '''
        Make a depth-first traversal of the glycan graph. Children are explored in descending bond-order.

        This is the default traversal method for all |Glycan| objects. :meth:`~.dfs` is an alias of this method.
        Both names can be used to specify this strategy to :meth:`~._get_traversal_method`.

        Parameters
        ----------
        from_node: None or Monosaccharide
            If `from_node` is |None|, then traversal starts from the root node. Otherwise it begins
            from the given node.
        apply_fn: function
            A function applied to each node on arrival. If this function returns a non-None value,
            the result is yielded from the generator, otherwise it is ignored. Defaults to :func:`.identity`
        visited: set or None
            A :class:`set` of node ID values to ignore. If |None|, defaults to the empty `set`

        Yields
        ------
        Return Value of `apply_fn`, by default |Monosaccharide|

        See also
        --------
        Glycan.breadth_first_traversal
        '''
        sort_predicate = methodcaller("order")
        node_stack = list([self.root if from_node is None else from_node])
        visited = set() if visited is None else visited
        while len(node_stack) > 0:
            node = node_stack.pop()
            visited.add(node.id)
            res = apply_fn(node)
            if res is not None:
                yield res
            node_stack.extend(sorted((terminal for link in node.links.values()
                                      for terminal in link if terminal.id not in visited), key=sort_predicate))

    # Convenience aliases and the set up the traversal_methods entry
    dfs = depth_first_traversal
    traversal_methods['dfs'] = "dfs"
    traversal_methods['depth_first_traversal'] = "dfs"

    def breadth_first_traversal(
            self, from_node=None, apply_fn=identity, visited=None):
        '''
        Make a breadth-first traversal of the glycan graph. Children are explored in descending bond-order.

        :meth:`~.bfs` is an alias of this method.
        Both names can be used to specify this strategy to :meth:`~._get_traversal_method`.

        Parameters
        ----------
        from_node: None or Monosaccharide
            If `from_node` is |None|, then traversal starts from the root node. Otherwise it begins
            from the given node.
        apply_fn: function
            A function applied to each node on arrival. If this function returns a non-None value,
            the result is yielded from the generator, otherwise it is ignored. Defaults to :func:`.identity`
        visited: set or None
            A :class:`set` of node ID values to ignore. If |None|, defaults to the empty `set`

        Yields
        ------
        Return Value of `apply_fn`, by default |Monosaccharide|

        See also
        --------
        Glycan.depth_first_traversal
        '''
        sort_predicate = methodcaller("order")
        node_queue = deque([self.root if from_node is None else from_node])
        visited = set() if visited is None else visited
        while len(node_queue) > 0:
            node = node_queue.popleft()
            visited.add(node.id)
            res = apply_fn(node)
            if res is not None:
                yield res
            node_queue.extend(sorted((terminal for link in node.links.values()
                                      for terminal in link if terminal.id not in visited), key=sort_predicate))

    # Convenience aliases and the set up the traversal_methods entry
    bfs = breadth_first_traversal
    traversal_methods['bfs'] = "bfs"
    traversal_methods['breadth_first_traversal'] = "bfs"

    def _get_traversal_method(self, method):
        if method == 'dfs':
            return self.dfs
        elif method == 'bfs':
            return self.bfs
        elif isinstance(method, Callable):
            return partial(method, self)
        traversal = self.traversal_methods.get(method, None)
        if traversal is None:
            raise AttributeError("Unknown traversal method: {}".format(method))
        traversal = getattr(self, traversal)
        return traversal

    def __iter__(self):
        return self.iternodes()

    def iternodes(
            self, from_node=None, apply_fn=identity, method='dfs', visited=None):
        '''
        Generic iterator over nodes. :meth:`Glycan.__iter__` is an alias of this method

        Parameters
        ----------
        from_node: None or Monosaccharide
            If `from_node` is |None|, then traversal starts from the root node. Otherwise it begins
            from the given node.
        apply_fn: function
            A function applied to each node on arrival. If this function returns a non-None value,
            the result is yielded from the generator, otherwise it is ignored. Defaults to :func:`.identity`
        method: str or `function`
            Traversal method to use. See :meth:`._get_traversal_method`
        visited: set or None
            A :class:`set` of node ID values to ignore. If |None|, defaults to the empty `set`

        Yields
        ------
        Return Value of `apply_fn`, by default Monosaccharide

        See also
        --------
        depth_first_traversal
        breadth_first_traversal
        _get_traversal_method
        '''
        traversal = self._get_traversal_method(method)
        return traversal(
            from_node=from_node, apply_fn=apply_fn, visited=visited)

    def iterlinks(
            self, apply_fn=identity, substituents=False, method='dfs', visited=None):
        '''
        Iterates over all |Link| objects in |Glycan|.

        Parameters
        ----------
        substituents: bool
            If `substituents` is |True|, then include the |Link| objects in
            :attr:`substituent_links` on each |Monosaccharide|
        method: str or function
            The traversal method controlling the order of the nodes visited
        visited: None or set
            The collection of id values to ignore when traversing

        Yields
        ------
        Link
        '''
        traversal = self._get_traversal_method(method)
        links_visited = set()

        def links(obj):
            if substituents:
                for pos, link in obj.substituent_links.items():
                    res = apply_fn((pos, link))
                    if res:
                        yield res
            for pos, link in obj.links.items():
                if link.id in links_visited:
                    continue
                links_visited.add(link.id)
                res = apply_fn((pos, link))
                if res:
                    yield res

        return itertools.chain.from_iterable(
            traversal(apply_fn=links, visited=visited))

    def leaves(self, bidirectional=False, method='dfs', visited=None):
        '''
        Iterates over all |Monosaccharide| objects in |Glycan|, yielding only those
        that have no child nodes.

        Parameters
        ----------
        bidirectional: bool
            If `bidirectional` is |True|, then only |Monosaccharide| objects
            with only one entry in :attr:`links`.
        method: str or function
            The traversal method controlling the order of the nodes visited
        visited: None or set
            The collection of id values to ignore when traversing

        Yields
        ------
        |Monosaccharide|
        '''
        traversal = self._get_traversal_method(method)
        if bidirectional:
            def is_leaf(obj):
                if len(obj.links) == 1:
                    yield obj
        else:
            def is_leaf(obj):
                if len(list(obj.children())) == 0:
                    yield obj

        return itertools.chain.from_iterable(
            traversal(apply_fn=is_leaf, visited=visited))

    def label_branches(self):
        '''
        Labels each branch point with an alphabetical symbol. Also computes and stores
        each branch's length and stores it in :attr:`branch_lengths`. Sets :attr:`branch_lengths`
        of `self` and :attr:`Link.label` for each link attached to `self`.
        '''
        last_branch_label = MAIN_BRANCH_SYM
        self.branch_lengths = defaultdict(int)
        branch_parent_map = {}

        def parent_link_symbol(node):
            try:
                label = node.links[node.parents().next()[0]][0].label
                if label is None:
                    return MAIN_BRANCH_SYM
                else:
                    return label[0]
            except StopIteration:
                return MAIN_BRANCH_SYM

        for node in self:
            links = []
            for link in node.links.values():
                if link.is_child(node):
                    continue
                links.append(link)

            if len(links) == 1:
                label_key = parent_link_symbol(node)
                self.branch_lengths[label_key] += 1
                label = "{}{}".format(
                    label_key, self.branch_lengths[label_key])
                links[0].label = label
            else:
                last_label_key = label_key = parent_link_symbol(node)
                count = self.branch_lengths[last_label_key]
                for link in links:
                    # Does not handle small branches correctly.
                    # if len(list(link.child.children())) < 2:
                    #     label_key = parent_link_symbol(node)
                    #     self.branch_lengths[label_key] += 1
                    #     label = "{}{}".format(
                    #         label_key, self.branch_lengths[label_key])
                    #     link.label = label
                    # else:
                    last_branch_label = chrinc(
                        last_branch_label) if last_branch_label != MAIN_BRANCH_SYM else 'a'
                    new_label_key = last_branch_label
                    branch_parent_map[new_label_key] = last_label_key
                    self.branch_lengths[new_label_key] = count + 1
                    label = "{}{}".format(
                        new_label_key, self.branch_lengths[new_label_key])
                    link.label = label
        # Update parent branch lengths
        longest = 0
        for branch in sorted(list(self.branch_lengths.keys()), reverse=True):
            if branch == '-':
                continue
            length = self.branch_lengths[branch]
            longest = max(longest, length)
            parent = branch_parent_map[branch]
            self.branch_lengths[parent] = max(length, self.branch_lengths[parent])
        self.branch_parent_map = branch_parent_map
        self.branch_lengths["-"] = longest

    def count_branches(self):
        '''
        Count the number of branches in the Glycan tree

        Returns
        -------
        int
        '''
        count = 0
        for node in self:
            if len(node.links) > 2:
                count += 2 if count == 0 else 1
        return count

    def order(self):
        '''
        The number of nodes in the graph. :meth:`__len__` is an alias of this

        Returns
        -------
        int
        '''
        count = 0
        for node in self:
            count += 1
        return count

    __len__ = order

    def to_glycoct(self, buffer=None, close=False):
        '''
        Serialize the |Glycan| graph object into condensed GlycoCT, using
        `buffer` to store the result. If `buffer` is |None|, then the
        function will operate on a newly created :class:`~glypy.utils.StringIO` object.

        Parameters
        ----------
        buffer: file-like or None
            The stream to write the serialized structure to. If |None|, uses an instance
            of `StringIO`
        close: bool
            Whether or not to close the stream in `buffer` after writing is done

        Returns
        -------
        file-like or str if ``buffer`` is :const:`None`

        '''
        is_stringio = False
        if buffer is None:
            buffer = StringIO()
            is_stringio = True

        buffer.write("RES\n")

        res_counter = make_counter()
        lin_counter = make_counter()

        # Look-ups for mapping RES nodes to objects by section index and id,
        # respectively
        index_to_residue = {}
        residue_to_index = {}

        # Accumulator for linkage indices and mapping linkage indices to
        # dependent RES indices
        lin_accumulator = []
        dependencies = defaultdict(dict)

        # Detect cycles and avoid including the same residue twice
        visited = set()

        for node in (self):
            if node.id in visited:
                continue
            visited.add(node.id)
            try:
                res, lin, index = node.to_glycoct(
                    res_counter, lin_counter, complete=False)

                lin_accumulator.append((index, lin))
                residue_to_index[node.id] = index
                index_to_residue[index] = node

                for pos, link in node.links.items():
                    if link.is_child(node):
                        continue
                    dependencies[link.child.id][node.id] = ((lin_counter(), link))
                for line in res:
                    buffer.write(line + '\n')
            except:
                pass
        buffer.write("LIN\n")
        for res_ix, links in lin_accumulator:
            for line in links:
                buffer.write(line + '\n')
            residue = index_to_residue[res_ix]
            for pos, link in residue.links.items():
                if link.is_child(residue):
                    continue
                child_res = link.child
                ix, link = dependencies[child_res.id][residue.id]
                buffer.write(
                    link.to_glycoct(ix, res_ix, residue_to_index[child_res.id]) + "\n")

        if is_stringio:
            return buffer.getvalue()
        else:  # pragma: no cover
            if close:
                buffer.close()
            return buffer

    __repr__ = to_glycoct

    def mass(self, average=False, charge=0, mass_data=None):
        '''
        Calculates the total mass of the intact graph by querying each
        node for its mass.

        Parameters
        ----------
        average: bool
            Whether or not to use the average isotopic composition when calculating masses.
            When ``average == False``, masses are calculated using monoisotopic mass.
        charge: int
            If charge is non-zero, m/z is calculated, where m is the theoretical mass, and z is `charge`
        mass_data: dict
            If mass_data is None, standard NIST mass and isotopic abundance data are used. Otherwise the
            contents of mass_data are assumed to contain elemental mass and isotopic abundance information.


        Returns
        -------
        float

        See also
        --------
        :func:`glypy.composition.composition.calculate_mass`
        '''
        return sum(
            node.mass(average=average, charge=charge, mass_data=mass_data) for node in self)

    def total_composition(self):
        '''
        Computes the sum of the composition of all |Monosaccharide| objects in ``self``

        Returns
        -------
        :class:`~glypy.composition.Composition`
        '''

        return sum((node.total_composition() for node in self), Composition())

    def clone(self, index_method='dfs', visited=None):
        '''
        Create a copy of `self`, indexed using `index_method`, a *traversal method*  or |None|.

        Returns
        -------
        :class:`~glypy.structure.glycan.Glycan`
        '''
        clone_root = graph_clone(self.root, visited=visited)
        duplicate = Glycan(clone_root, index_method=index_method)

        return duplicate

    def __eq__(self, other):
        '''
        Two glycans are considered equal if they are identically ordered nodes.

        Parameters
        ----------
        self, other: :class:`~glypy.structure.glycan.Glycan`

        Returns
        -------
        bool

        See also
        --------
        :meth:`glypy.structure.Monosaccharide.exact_ordering_equality`
        '''
        if other is None:
            return False
        elif not isinstance(other, Glycan):
            return False
        return self.root.exact_ordering_equality(other.root)

    def topological_equality(self, other):
        '''
        Two glycans are considered equal if they are topologically equal.

        Parameters
        ----------
        self: :class:`Glycan`
        other: :class:`Glycan`

        Returns
        -------
        bool

        See also
        --------
        :meth:`glypy.structure.Monosaccharide.topological_equality`
        '''
        return self.root.topological_equality(other.root)

    def __ne__(self, other):
        return not self == other

    def substructures(self, max_cleavages=1, min_cleavages=1, inplace=False):
        '''
        Generate disjoint subtrees from this glycan by examining by removing one or
        more monosaccharide-monosaccharide bond.


        Parameters
        ----------
        max_cleavages: |int|
            The maximum number of bonds to break per substructure
        min_cleavages: |int|
            The minimum number of bonds to break per substructure
        min_size: |int|
            The minimum number of monosaccharides per substructure

        See also
        --------
        :func:`glypy.composition.composition.calculate_mass`
        '''
        structure = self
        if not inplace:
            structure = self.clone()
        for frag in structure.break_links_subtrees(
                max_cleavages):
            yield frag

    def name_fragment(self, fragment):
        '''
        Attempt to assign a full name to a fragment based on the branch and position relative to
        the reducing end along side A/B/C/X/Y/Z, according to :title-reference:`Domon and Costello`
        '''

        break_targets = fragment.link_ids
        crossring_targets = fragment.crossring_cleavages

        # Accumulator for name components
        name_parts = []
        # Collect cross-ring fragment names
        for crossring_id in crossring_targets:
            # Seek the link that holds the fragmented residue
            for link in self.link_index:
                if link.child.id == crossring_id:
                    ion_type = crossring_targets[crossring_id]
                    label = link.label
                    if _fragment_direction[ion_type[1]] > 0:
                        name = "{}{}".format(
                            ''.join(map(str, ion_type)), label.replace(MAIN_BRANCH_SYM, ""))
                        name_parts.append(name)
                    else:
                        label_key = label[0]
                        distance = int(label[1:])
                        inverted_distance = self.branch_lengths[
                            label_key] - (distance - 1)
                        name = "{}{}{}".format(
                            ''.join(map(str, ion_type)), label_key.replace(MAIN_BRANCH_SYM, ""), inverted_distance)
                        name_parts.append(name)

        # Collect glycocidic fragment names
        for break_id, ion_type in break_targets.items():
            ion_type = ion_type[1]
            if _fragment_direction[ion_type] > 0:
                link = self.link_index[break_id - 1]
                label = link.label
                name = "{}{}".format(
                    ion_type,
                    label.replace(
                        MAIN_BRANCH_SYM,
                        ""))
                name_parts.append(name)
            else:
                link = self.link_index[break_id - 1]
                label = link.label
                label_key = label[0]
                distance = int(label[1:])
                inverted_distance = self.branch_lengths[
                    label_key] - (distance - 1)
                name = "{}{}{}".format(
                    ion_type, label_key.replace(MAIN_BRANCH_SYM, ""), inverted_distance)
                name_parts.append(name)

        return '-'.join(sorted(name_parts))

    def break_links_subtrees(self, n_links):
        """Iteratively generate all subtrees from glycosidic bond cleavages, creating all
        :math:`2{L \choose n}` subtrees.

        Parameters
        ----------
        n_links : int
            Number of links to break simultaneously

        Yields
        ------
        Subtree
        """
        links = list(self.link_index)
        for breaks in itertools.combinations(links, n_links):

            subtrees = []
            for link in breaks:
                parent, child = link.break_link(refund=True)
                parent_tree = Glycan(parent, index_method=None)
                child_tree = Glycan(child, index_method=None)
                subtrees.append(parent_tree)
                subtrees.append(child_tree)

            unique_subtrees = []
            for subtree in subtrees:
                ids = {n.id for n in subtree}
                for uids, unique in unique_subtrees:
                    if ids == uids:
                        break
                else:
                    unique_subtrees.append((ids, subtree))

            for ids, subtree in unique_subtrees:
                subtree = subtree.clone(
                    index_method=None).reroot(
                    index_method=None)
                include_nodes = {n.id for n in subtree}

                link_ids = [link.id for link in breaks
                            if link.parent.id in include_nodes or
                            link.child.id in include_nodes]

                parent_break_ids = {link.id: link.parent.id for link in breaks
                                    if link.parent.id in include_nodes}

                child_break_ids = {link.id: link.child.id for link in breaks
                                   if link.child.id in include_nodes}

                yield Subtree(subtree, include_nodes, link_ids, parent_break_ids, child_break_ids)

            for link in breaks:
                link.apply()

    def crossring_subtrees(self, n_links):
        """Generate all combinations of cross ring fragments and
        glycosidic cleavages, cleaving between 1 and `n_links`
        monosaccharides paired with `n_links` - 1 to 0 glycosidic cleavages.

        Parameters
        ----------
        n_links : int
            Total number of breaks to create, between cross ring cleavages and
            complemenatary glycosidic cleavages.

        Yields
        ------
        Subtree
        """
        links = list(self.link_index)
        # Localize globals
        _str = str
        # Break at least one ring
        for i in range(1, n_links + 1):
            # Generate all combinations of i rings to break
            for link_combination in itertools.combinations(links, i):
                # Creates a list of lists of CrossRingPairs, each inner list for a
                # single Monosaccharide
                crossring_combinations = [
                    CrossRingPair.from_link(link) for link in link_combination]
                # Combinations are splatted to unwrap the outer container so that
                # the inner lists are multiplexed.
                for breaks in itertools.product(*crossring_combinations):
                    subtrees = []
                    for ring in breaks:
                        parent, child = ring.break_link()
                        parent_tree = Glycan(parent, index_method=None)
                        child_tree = Glycan(child, index_method=None)
                        subtrees.append(parent_tree)
                        subtrees.append(child_tree)

                    unique_subtrees = []
                    for tree in subtrees:
                        ids = {n.id for n in tree}
                        for uids, unique in unique_subtrees:
                            if ids == uids:
                                break
                        else:
                            unique_subtrees.append((ids, tree))

                    # If tis iteration hasn't broken all n rings, there are some
                    # breaks left over to be generated in glycosidic cleavages.
                    # Generate all possible glycosidic cleavage subtrees of the
                    # generated cross ring cleavage subtrees.
                    if n_links - i > 0:
                        for ids, subtree in unique_subtrees:
                            partitions = subtree.break_links_subtrees(n_links - i)
                            for part in partitions:
                                included_crossring = {}
                                for crossring in breaks:
                                    if crossring.id in part.include_nodes:
                                        xring_residue = part.tree.get(
                                            crossring.id)
                                        try:
                                            included_crossring[xring_residue.id] = (','.join(
                                                (_str(xring_residue.cleave_1), _str(xring_residue.cleave_2))),
                                                xring_residue.kind)
                                        except AttributeError:
                                            pass
                                            # Not a CrossRingFragment in this instance

                                part.crossring_cleavages = included_crossring
                                yield part
                    else:
                        for ids, subtree in unique_subtrees:
                            subtree = subtree.reroot(
                                index_method=None).clone(
                                index_method=None)

                            included_crossring = {}
                            include_nodes = {n.id for n in subtree}
                            for crossring in breaks:
                                if crossring.id in include_nodes:
                                    xring_residue = subtree.get(crossring.id)
                                    try:
                                        included_crossring[xring_residue.id] = (','.join(
                                            (_str(xring_residue.cleave_1), _str(xring_residue.cleave_2))),
                                            xring_residue.kind)
                                    except AttributeError:
                                        pass
                                        # Not a CrossRingFragment in this instance
                            yield Subtree(subtree, include_nodes, {}, {}, {}, crossring_cleavages=included_crossring)
                    # Re-join the ring, retracting all links from the crossring objects
                    for ring in breaks:
                        ring.apply()
                    # Clean up any lingering links trapped in the closure created by adjacent
                    # crossring pairs that were made visible by `ring.apply()`
                    for ring in breaks:
                        ring.release()

    def fragments(self, kind="BY", max_cleavages=1, average=False, charge=0, mass_data=None):
        '''
        Generate carbohydrate backbone fragments from this glycan by examining the disjoint subtrees
        created by removing one or more monosaccharide-monosaccharide bond.


        Parameters
        ----------
        kind: `sequence`
            Any `iterable` or `sequence` of characters corresponding to A/B/C/X/Y/Z
            as published by :title-reference:`Domon and Costello`
        max_cleavages: |int|
            The maximum number of bonds to break per fragment
        average: bool, optional, defaults to `False`
            Whether or not to use the average isotopic composition when calculating masses.
            When ``average == False``, masses are calculated using monoisotopic mass.
        charge: int, optional, defaults to 0
            If charge is non-zero, m/z is calculated, where m is the theoretical mass, and z is `charge`
        mass_data: dict, optional, defaults to `None`
            If mass_data is |None|, standard NIST mass and isotopic abundance data are used. Otherwise the
            contents of `mass_data` are assumed to contain elemental mass and isotopic abundance information.

        Yields
        ------
        :class:`Fragment`

        See also
        --------
        :func:`glypy.composition.composition.calculate_mass`
        :meth:`subtrees`
        :meth:`crossring_subtrees`
        :meth:`.Subtree.to_fragments`
        '''
        seen = set()
        source = self.clone()
        for i in range(1, max_cleavages + 1):
            gen = source.break_links_subtrees(i)
            if len(set("AX") & set(kind)) > 0:
                gen = itertools.chain(
                    gen,
                    source.crossring_subtrees(i))
            for subtree in gen:
                for fragment in subtree.to_fragments(kind, average=average,
                                                     charge=charge, mass_data=mass_data):
                    fragment.name = self.name_fragment(fragment)
                    if fragment.name in seen:
                        continue
                    else:
                        seen.add(fragment.name)
                    yield fragment

    def subtrees(self, max_cleavages=1, include_crossring=False):
        '''
        Generate subtrees from this tree by breaking `max_cleavages` bonds or rings.

        Parameters
        ----------
        max_cleavages: int
            The maximum number of bonds to break per fragment
        include_crossring: bool
            Whether to include cross ring cleavages

        Yields
        ------
        Subtree
        '''
        source = self.clone()
        for i in range(1, max_cleavages + 1):
            gen = source.break_links_subtrees(i)
            if include_crossring:
                gen = itertools.chain(
                    gen,
                    source.crossring_subtrees(i))
            for subtree in gen:
                yield subtree
