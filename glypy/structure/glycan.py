'''
Represent a sugar graph with pseudo-directed edges.

'''
import warnings
import operator
import logging
import itertools
from functools import partial
from collections import deque, defaultdict
try:
    from collections.abc import Callable
except ImportError:
    from collections import Callable

from glypy.utils import (
    identity,
    chrinc,
    uid)
from glypy.composition import Composition

from .base import SaccharideCollection
from .monosaccharide import (
    Monosaccharide,
    graph_clone,
    toggle as residue_toggle,
    MonosaccharideOccupancy)
from .constants import UnknownPosition, NoPosition
from .substituent import Substituent
from .crossring_fragments import crossring_fragments, CrossRingPair
from .fragment import Subtree

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

    Parameters
    ----------
    fragment: :class:`~.GlycanFragment`
        The :class:`~.GlycanFragment` to extract substructure for.
    tree: :class:`Glycan`
        The |Glycan| to extract substructure from.

    Returns
    -------
    :class:`Glycan`:
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
        next(residue_toggle(target))

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


class Glycan(SaccharideCollection):
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
    branch_parent_map: |dict|
        A dictionary mapping branch symbols to their parent branch symbols
    '''

    verbose = False

    _serializers = {}

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

    def fragment_to_substructure(self, fragment):
        """Extract the substructure of `tree` which is contained in `fragment`

        Parameters
        ----------
        fragment: :class:`~.GlycanFragment`
            The :class:`~.GlycanFragment` to extract substructure for.
        tree: :class:`Glycan`
            The |Glycan| to extract substructure from.

        Returns
        -------
        :class:`Glycan`:
            The |Glycan| substructure defined by the nodes contained in `fragment` as
            found in `tree`
        """
        return fragment_to_substructure(fragment, self)

    traversal_methods = {}

    def __init__(self, root=None, index_method='dfs', canonicalize=False):
        '''
        Constructs a new :class:`Glycan` from the collection of connected :class:`~.Monosaccharide`
        objects rooted at `root`.

        If :obj:`index_method` is not |None|, the graph is indexed by the default search method
        given by :obj:`traversal_methods[index_method]`
        '''
        if root is None:
            root = Monosaccharide()
        self.root = root
        self.index = []
        self.link_index = []
        self.branch_lengths = {}
        if index_method is not None:
            self.reindex(method=index_method)
        if canonicalize:
            self.canonicalize()

    def has_index(self):
        return bool(self.index)

    def reindex(self, method='dfs'):
        '''
        Traverse the graph using the function specified by `method`. The order of
        traversal defines the new :attr:`id` value for each |Monosaccharide|
        and |Link|.

        The order of traversal also defines the ordering of the |Monosaccharide|
        in :attr:`index` and |Link| in :attr:`link_index`.

        Prior to constructing a |Glycan| instance, component |Monosaccharide| instances
        may be labeled, converting their id field into a tuple.

        Calls :meth:`label_branches` after indexing is complete.

        Returns
        -------
        Glycan:
            self

        See Also
        --------
        deindex
        label_branches
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
            node.id = i
            i += 1
            # reindex substituents as well
            try:
                for j, subst in node.substituents():
                    subst.id = i
                    i += 1
            except AttributeError:
                if node.node_type is Substituent.node_type:
                    continue
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

    def _build_link_index(self, method='dfs'):
        link_index = []
        for pos, link in self.iterlinks(method=method):
            link_index.append(link)
        self.link_index = link_index

    def _build_node_index(self, method='dfs'):
        index = []
        for node in self.iternodes(method=method):
            index.append(node)
        self.index = index

    def deindex(self):
        '''
        When combining two Glycan structures, very often their component ids will
        overlap, making it impossible to differentiate between a cycle and the new
        graph. This function mangles all of the node and link ids so that they are
        distinct from the pre-existing nodes.

        Returns
        -------
        Glycan:
            self
        '''
        if self.index is not None and len(self.index) > 0:
            base = uid()
            for node in self.index:
                node.id += base
                node.id *= -1
                try:
                    for j, subst in node.substituents():
                        subst.id += base
                        subst.id *= -1
                except AttributeError:
                    if node.node_type is Substituent.node_type:
                        warnings.warn("A substituent has been detected as a parent "
                                      "node. This glycan may not fragment correctly")
            for link in self.link_index:
                link.id += base
                link.id *= -1
        return self

    def reroot(self, index_method='dfs'):
        """Set :attr:`root` to the node with the lowest :attr:`id`.

        Should only be used if the glycan has been indexed.

        Parameters
        ----------
        index_method : str, optional
            The name of the index method to use to reindex the glycan relative to
            the new root node. If :const:`None`, no reindexing is done. The default
            is 'dfs'

        Returns
        -------
        Glycan:
            self
        """
        self.root = sorted(iter(self), key=operator.attrgetter('id'))[0]
        if index_method is not None:
            self.reindex(method=index_method)
        return self

    def initialize_structure(self, method='dfs'):
        """Rebuild the structure index from scratch using :meth:`reindex`,
        then sets the traversal order using :meth:`canonicalize`.

        Use this method when building a :class:`Glycan` instance from a manually
        created :class:`~.Monosaccharide` graph.
        """
        self.canonicalize()
        self.reindex(method=method)
        return self

    def __getitem__(self, ix):
        '''
        Fetch a :class:`~.Monosaccharide` from :attr:`index`.

        Returns
        -------
        Monosaccharide

        Raises
        ------
        IndexError:
            If the provided `ix` exceeds the length of the index,
            or if :attr:`index` has not been populated.
        '''
        if self.index is not None:
            return self.index[ix]
        else:
            raise IndexError(
                "Tried to access the index of an unindexed Glycan.")

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, state):
        self.__dict__.update(state)

    def __root__(self):
        return self.root

    def __tree__(self):
        return self

    def get(self, ix):
        """Get a :class:`~.Monosaccharide` from this structure by its :attr:`~.Monosaccharide.id` value.

        If :attr:`index` is populated it will be iterated over, otherwise :meth:`__iter__` will be called.

        Parameters
        ----------
        ix : int or tuple
            The id value to search for.

        Returns
        -------
        Monosaccharide

        Raises
        ------
        IndexError:
            If the id value is not found
        """
        if self.index:
            iterable = self.index
        else:
            iterable = self
        for node in iterable:
            if node.id == ix:
                return node
        raise IndexError(
            "Could not find a node with the given id {}".format(ix))

    def get_link(self, ix):
        """Search for a :class:`~.Link` by :attr:`~.Link.id` value.

        This will use :meth:`iterlinks` to iterate over the linkages
        in the structure

        Parameters
        ----------
        ix : int
            The link index to search for.

        Returns
        -------
        :class:`~.Link`

        Raises
        ------
        IndexError:
            if the id value is not found
        """
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
        try:
            return self.root.reducing_end
        except AttributeError:
            if self.root.node_type is Substituent.node_type:
                return None
            else:
                raise

    def set_reducing_end(self, value):
        '''
        Sets the reducing end type, and configures the :attr:`root` |Monosaccharide| appropriately.

        If the reducing_end is not |None|, then the following state changes are made to :attr:`root`:

        .. code-block:: python

            self.root.ring_start = 0
            self.root.ring_end = 0
            self.root.anomer = "uncyclized"

        Else, the correct state is unknown:

        .. code-block:: python

            self.root.ring_start = UnknownPosition
            self.root.ring_end = UnknownPosition
            self.root.anomer = None

        .. note:: This method is called automatically when setting :attr:`reducing_end`, and does not

        need to be used explicitly.
        '''
        if self.root.node_type is Substituent.node_type:
            raise TypeError("Cannot set reducing end on this glycan composition")
        self.root.reducing_end = value
        if self.reducing_end is not None:
            self.root.ring_start = 0
            self.root.ring_end = 0
            self.root.anomer = "uncyclized"
        else:
            self.root.ring_start = UnknownPosition
            self.root.ring_end = UnknownPosition
            self.root.anomer = None

    @reducing_end.setter
    def reducing_end(self, value):
        self.set_reducing_end(value)

    def depth_first_traversal(self, from_node=None, apply_fn=identity, visited=None):
        '''
        Make a depth-first traversal of the glycan graph. Children are explored in descending bond-order.

        This is the default traversal method for all |Glycan| objects. :meth:`~.dfs` is an alias of this method.
        Both names can be used to specify this strategy to :meth:`~._get_traversal_method`.

        When selecting an iteration strategy, this strategy is specified as "dfs".

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
        node_stack = deque([self.root if from_node is None else from_node])
        visited = set() if visited is None else visited
        while len(node_stack) > 0:
            node = node_stack.pop()
            visited.add(node.id)
            if apply_fn is identity:
                yield node
            else:
                res = apply_fn(node)
                if res is not None:
                    yield res

            for link in node.links.values():
                terminal = link.parent
                if terminal.id not in visited:
                    node_stack.append(terminal)
                terminal = link.child
                if terminal.id not in visited:
                    node_stack.append(terminal)

    traversal_methods['dfs'] = "depth_first_traversal"
    traversal_methods['depth_first_traversal'] = "depth_first_traversal"

    def breadth_first_traversal(self, from_node=None, apply_fn=identity, visited=None):
        '''
        Make a breadth-first traversal of the glycan graph. Children are explored in descending bond-order.

        When selecting an iteration strategy, this strategy is specified as "bfs".

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
        node_queue = deque([self.root if from_node is None else from_node])
        visited = set() if visited is None else visited
        while len(node_queue) > 0:
            node = node_queue.popleft()
            visited.add(node.id)

            if apply_fn is identity:
                yield node
            else:
                res = apply_fn(node)
                if res is not None:
                    yield res
            # node_queue.extend(terminal for link in node.links.values()
            #                   for terminal in link if terminal.id not in visited)
            for link in node.links.values():
                terminal = link.parent
                if terminal.id not in visited:
                    node_queue.append(terminal)
                terminal = link.child
                if terminal.id not in visited:
                    node_queue.append(terminal)

    traversal_methods['bfs'] = "breadth_first_traversal"
    traversal_methods['breadth_first_traversal'] = "breadth_first_traversal"

    def indexed_traversal(self, from_node=None, apply_fn=identity, visited=None):
        """Traverse the glycan structure along :attr:`index`.

        This is substantially faster than other traversal methods for complete traversals
        at the cost of a) requiring a call to :meth:`reindex` to populate :attr:`index`
        if it has not been called, and b) is not automatically updated if the structure is
        modified.

        When selecting an iteration strategy, this strategy is specified as "index".

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

        See Also
        --------
        reindex
        """
        if not self.index:
            self._build_node_index()
        if from_node is None and apply_fn is identity:
            for node in self.index:
                yield node
        else:
            i = 0
            n = len(self.index)
            if from_node is not None:
                while i < n:
                    node = self.index[i]
                    if node == from_node:
                        break
                    i += 1
            while i < n:
                node = self.index[i]
                if apply_fn is identity:
                    yield node
                else:
                    value = apply_fn(node)
                    if value is not None:
                        yield value

                i += 1

    traversal_methods['index'] = "indexed_traversal"

    def _get_traversal_method(self, method):
        """An internal helper method used to resolve traversal
        methods by name or alias.

        Parameters
        ----------
        method : :class:`str` or :class:`Callable`
            If a :class:`str`, it is looked up in the class-level
            :attr:`traversal_methods` dictionary and the name of
            the appropriate method is retrieved with :func:`getattr`
            and returned.

            If a :class:`Callable`, the function's first parameter is
            bound to `self` and returned.

        Returns
        -------
        :class:`Callable`
        """
        if method == 'dfs':
            return self.depth_first_traversal
        elif method == 'bfs':
            return self.breadth_first_traversal
        elif isinstance(method, Callable):
            return partial(method, self)
        try:
            traversal = self.traversal_methods[method]
        except KeyError:
            raise KeyError("Unknown traversal method: {}".format(method))
        traversal = getattr(self, traversal)
        return traversal

    def __iter__(self):
        return self.depth_first_traversal()

    def iternodes(self, from_node=None, apply_fn=identity, method='dfs', visited=None):
        '''
        Generic iterator over nodes dispatching to a strategy given by `method`, defaulting
        to :meth:`~.depth_first_traversal`.

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

    def iterlinks(self, apply_fn=identity, substituents=False, method='dfs', visited=None):
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

    def canonicalize(self, canonicalizer=None, **kwargs):
        """Canonicalize this glycan, sorting the order in which its links from the
        same monosaccharide are traversed.

        This currently uses the the :title-reference:`GlycoCT` canonicalization
        algorithm.

        Parameters
        ----------
        canonicalizer : subclass of :class:`~.CanonicalizerBase`, optional
            The canonicalization algorithm to use
        **kwargs
            Forwarded to the canonicalizer

        Returns
        -------
        :class:`~.Glycan`
            This glycan, reordered in place.
        """
        from glypy.algorithms.canonicalize import canonicalize
        return canonicalize(self, canonicalizer=canonicalizer, **kwargs)

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
        :class:`~.Monosaccharide`
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

    def ambiguous_links(self):
        """Locate all links which are :class:`.AmbiguousLink` objects

        Returns
        -------
        list
            list of ambiguous links
        """
        ambiguous_links = []

        # find all ambiguous links
        for i, link in self.iterlinks():
            if link.is_ambiguous():
                ambiguous_links.append(link)
        return ambiguous_links

    def iterconfiguration(self):
        '''Iterate over all valid configurations of ambiguous linkages.

        During calculation, the :class:`~.AmbiguousLink` objects may be mutated, but
        by the time a new configuration is yielded all changes should be reversed. If an
        error occurs during configuration adjustment, it may not be possible to restore the
        object to its original state.

        Yields
        ------
        :class:`tuple` of (:class:`~.AmbiguousLink`, :class:`~.Monosaccharide`,
                           :class:`~.Monosaccharide`, :class:`int`, :class:`int`)
            The ambiguous link, the parent chosen, the child chosen, the parent linkage site chose, and the child
            linkage site chosen

        Examples
        --------
        >>> from glypy.io import glyspace
        >>> structure_record = glyspace.get("G81339YK")
        >>> structure = structure_record.structure_
        >>> configurations = []
        >>> for config_list in structure.iterconfiguration():
        ...     instance = structure.clone()
        ...     for link, conf in config_list:
        ...         link = instance.get_link(link.id)
        ...         parent = instance.get(conf[0].id)
        ...         child = instance.get(conf[1].id)
        ...         link.reconfigure(parent, child, conf[2], conf[3])
        ...     configurations.append(instance)
        >>> len(configurations)
        4

        See Also
        --------
        :meth:`~.AmbiguousLink.iterconfiguration`, :meth:`~.AmbiguousLink.reconfigure`

        '''
        ambiguous_links = self.ambiguous_links()

        combos = itertools.product(*[
            link.iterconfiguration() for link in ambiguous_links])

        for configs in combos:
            n_configs = len(configs)
            # count the number of ambiguous links with unknown linkages
            has_unknowns = []
            for i, conf in enumerate(configs):
                if conf.parent_position == UnknownPosition or conf.child_position == UnknownPosition:
                    has_unknowns.append(i)
            if not has_unknowns:
                if len(set([(c.parent, c.parent_position) for c in configs])) < n_configs:
                    continue
                yield tuple(zip(ambiguous_links, configs))
            else:
                # mask the ambiguous links so the possible attachment sites can
                # be calculated
                for link in ambiguous_links:
                    link.break_link(refund=True)

                parent_map = {
                    parent.id: MonosaccharideOccupancy.build(parent)
                    for parent, _, _, _ in configs
                }

                valid = True

                # place all concrete links
                for parent, child, parent_site, child_site in configs:
                    parent_occupancy = parent_map[parent.id]
                    if parent_site == UnknownPosition:
                        continue
                    if parent_site not in parent_occupancy.open_sites:
                        valid = False
                        break
                    # can't add more concrete links if the open sites are allocated
                    # for unknown localizations
                    elif len(parent_occupancy.open_sites) - parent_occupancy.unknown_sites <= 0:
                        valid = False
                        break
                    else:
                        # allocate the requested location so no new connections
                        # can be made there
                        ix = parent_occupancy.open_sites.index(parent_site)
                        parent_occupancy.open_sites.pop(ix)
                        parent_occupancy.occupied_sites.append(parent_site)
                if valid:
                    # place all unknown links
                    for parent, child, parent_site, child_site in configs:
                        parent_occupancy = parent_map[parent.id]
                        if parent_site != UnknownPosition:
                            continue
                        # can't add more unknown links if the open sites are allocated
                        # for unknown localizations
                        if len(parent_occupancy.open_sites) - parent_occupancy.unknown_sites <= 0:
                            valid = False
                            break
                        else:
                            parent_occupancy.unknown_sites += 1

                # re-attach ambiguous links so the state of the glycan returns to
                # its original state
                for link in ambiguous_links:
                    link.apply()

                # if configuration is valid, yield it
                if valid:
                    yield tuple(zip(ambiguous_links, configs))

    def has_undefined_linkages(self):
        """Check if this structure has undefined or ambiguous connectivity
        between its nodes.

        Returns
        -------
        bool:
            If any of its links are :class:`~.AmbiguousLink` instances, or
            have unknown positions (-1).
        """
        ambiguous_links = bool(self.ambiguous_links())
        if ambiguous_links:
            return ambiguous_links
        for node in self:
            if node.has_undefined_linkages():
                return True
        return False

    def label_branches(self):
        '''
        Labels each branch point with an alphabetical symbol. Also computes and stores
        each branch's length and stores it in :attr:`branch_lengths`. Sets :attr:`branch_lengths`
        of `self` and :attr:`Link.label` for each link attached to `self`. Also populates
        :attr:`branch_parent_map`.

        Branch symbols are increasing alphabetical characters. The root branch is denoted '-', though
        glycans having an :attr:`root` with multiple children will not have any actual branches with
        that label.

        :attr:`Link.label` updates use the current branch symbol, and the index of that link along
        that branch.

        .. note:: Labeling always uses a depth-first traversal of nodes.

        '''
        last_branch_label = MAIN_BRANCH_SYM
        self.branch_lengths = defaultdict(int)
        branch_parent_map = {}

        def parent_link_symbol(node):
            try:
                label = node.links[node.parents()[0][0]][0].label
                if label is None:
                    return MAIN_BRANCH_SYM
                else:
                    return label[0]
            except IndexError:
                return MAIN_BRANCH_SYM

        for node in self.depth_first_traversal():
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
                    last_branch_label = chrinc(
                        last_branch_label) if last_branch_label != MAIN_BRANCH_SYM else 'a'
                    new_label_key = last_branch_label
                    branch_parent_map[new_label_key] = last_label_key
                    self.branch_lengths[new_label_key] = count + 1
                    label = "{}{}".format(
                        new_label_key, self.branch_lengths[new_label_key])
                    link.label = label

        # Update parent branch lengths
        longest = self.branch_lengths[MAIN_BRANCH_SYM]
        for branch in sorted(list(self.branch_lengths.keys()), reverse=True):
            if branch == MAIN_BRANCH_SYM:
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

    def order(self, deep=False):
        '''
        The number of nodes in the graph. :meth:`__len__` is an alias of this

        Returns
        -------
        int
        '''
        count = 0
        for node in self.depth_first_traversal():
            count += 1
        return count

    def __len__(self):
        return self.order()

    @classmethod
    def register_serializer(cls, name, method):
        """Add `method` as `name` to the set of serializers to pick from
        in :meth:`serialize`

        Parameters
        ----------
        name : str
            The name of the serializer
        method : Callable
            A callable object that when called with a :class:`Glycan` returns a :class:`str`
        """
        cls._serializers[name] = method

    @classmethod
    def available_serializers(cls):
        """Get the list of available serialization formats

        Returns
        -------
        :class:`list` of :class:`str`
        """
        return list(cls._serializers.keys())

    def serialize(self, name='glycoct'):
        """Convert the structure to text.

        The serialization format is given by a :meth:`available_serializers`.

        Parameters
        ----------
        name : str, optional
            The name of the serialization format (the default is 'glycoct')

        Returns
        -------
        str
        """
        return self._serializers[name](self)

    def __repr__(self):
        """An alias for :meth:`serialize`.

        Returns
        -------
        str
        """
        return self.serialize()

    def mass(self, average=False, charge=0, mass_data=None, method='dfs'):
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
        if charge == 0:
            return sum(
                node.mass(average=average, charge=0, mass_data=mass_data) for node in self.iternodes(method=method))
        else:
            return self.total_composition().calc_mass(average=average, charge=charge, mass_data=mass_data)

    def total_composition(self, method='dfs'):
        '''
        Computes the sum of the composition of all |Monosaccharide| objects in ``self``

        Returns
        -------
        :class:`~glypy.composition.Composition`
        '''

        return sum((node.total_composition() for node in self.iternodes(method=method)), Composition())

    def clone(self, index_method='dfs', visited=None, cls=None):
        '''
        Create a copy of `self`, indexed using `index_method`, a *traversal method*  or |None|.

        Parameters
        ----------
        index_method: :class:`str`
            The indexing method to use when constructing the index of the copied
            structure
        visited: :class:`set`, optional
            A set of nodes to omit traversing through during the copying processing
        cls: :class:`type`
            A subclass of :class:`Glycan`, defaulting to :attr:`__class__`

        Returns
        -------
        :class:`~glypy.structure.glycan.Glycan`
        '''
        if cls is None:
            cls = self.__class__
        clone_root = graph_clone(self.root, visited=visited)
        duplicate = cls(root=clone_root, index_method=index_method)

        return duplicate

    def exact_ordering_equality(self, other):
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
        :term:`Exact Matching`
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
        :term:`Topological Matching`
        '''
        if other is None:
            return False
        elif not isinstance(other, Glycan):
            return False
        return self.root.topological_equality(other.root)

    def __eq__(self, other):
        """Test for exact ordering equality

        Parameters
        ----------
        other : :class:`Glycan`

        Returns
        -------
        :class:`bool`

        See Also
        --------
        :meth:`exact_ordering_equality`
        :term:`Exact Matching`
        """
        return self.exact_ordering_equality(other)

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        """Hashes the structure from the GlycoCT text representation

        See Also
        --------
        :meth:`serialize`
        """
        return hash(self.serialize("glycoct"))

    def substructures(self, max_cleavages=1, min_cleavages=1, inplace=False):
        '''
        Generate disjoint subtrees from this glycan by removing one or
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

        The formal grammar for fragment names in Backus-Naur Form:

        .. code:: xml

            <full-name>                ::= <fragment-name>|<fragment-name-list>
            <fragment-name>            ::= <glycosidic-fragment-name>|<crossring-fragment-name>
            <fragment-name-list>       ::= <fragment-name>"-"<fragment-name-list>|<fragment-name>
            <glycosidic-fragment-name> ::= <branch-identifier><fragment-type><index>
            <crossring-fragment-name>  ::= <ring-coordinates><fragment-type><branch-identifier><index>
            <fragment-type>            ::= "A" | "B" | "C" | "X" | "Y" | "Z"
            <ring-coordinate>          ::= <integer>,<integer>
            <index>                    ::= <integer>
            <integer>                  ::= <digit>|<integer><digit>
            <digit>                    ::= "0" | "1" | "2" | "3" | "4" | "5" | "6" | "7" | "8" | "9"
            <branch-identifier>        ::= <letter>|<letter><digit>|""
            <letter>                   ::= "a" | "b" | "c" | "d" | "e" | "f" | "g" | "h" | "i" | "j" |
                                           "k" | "l" | "m" | "n" | "o" | "p" | "q" | "r" | "s" | "t" |
                                           "u" | "v" | "w" | "x" | "y" | "z"

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
        r"""Iteratively generate all subtrees from glycosidic bond cleavages, creating all
        :math:`2{L \choose n}` subtrees.

        Parameters
        ----------
        n_links : int
            Number of links to break simultaneously

        Yields
        ------
        Subtree
        """
        if len(self.link_index) == 0:
            self._build_link_index()
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
        if len(self.link_index) == 0:
            self._build_link_index()

        links = list(self.link_index)
        # origin_mass = self.mass()
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
                        parent_tree._build_node_index()
                        child_tree = Glycan(child, index_method=None)
                        child_tree._build_node_index()
                        subtrees.append(parent_tree)
                        subtrees.append(child_tree)

                    unique_subtrees = []
                    for tree in subtrees:
                        ids = {n.id for n in tree.index}
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
                            subtree._build_link_index()
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
                            include_nodes = {n.id for n in subtree.indexed_traversal()}
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
                    while(any(ring.is_attached() for ring in breaks)):
                        for ring in breaks:
                            ring.release()
                    for ring in breaks:
                        ring.release()
                    # assert round(self.mass(), 4) == round(origin_mass, 4)

    def fragments(self, kind="BY", max_cleavages=1, average=False, charge=0, mass_data=None,
                  traversal_method='dfs'):
        '''
        Generate carbohydrate backbone fragments from this glycan by examining the disjoint subtrees
        created by removing one or more monosaccharide-monosaccharide bond.


        Parameters
        ----------
        kind: :class:`Iterable`
            Any :class:`Iterable` of characters corresponding to A/B/C/X/Y/Z
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
        :class:`GlycanFragment`

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
                                                     charge=charge, mass_data=mass_data,
                                                     traversal_method=traversal_method):
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


class NamedGlycan(Glycan):
    def __init__(self, name=None, *args, **kwargs):
        self.name = name
        super(NamedGlycan, self).__init__(*args, **kwargs)

    def clone(self, index_method='dfs', visited=None, cls=None):
        if cls is None:
            cls = NamedGlycan
        inst = super(NamedGlycan, self).clone(index_method=index_method, visited=visited, cls=cls)
        inst.name = self.name
        return inst

    def __repr__(self):
        rep = super(NamedGlycan, self).__repr__()
        return "%s\n%s" % (self.name, rep)
