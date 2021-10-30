# coding: utf-8
'''
A parser for :title-reference:`GlycoCT{condensed}` format.

:title-reference:`GlycoCT{condensed}` is a multi-line format for representing
glycan structures and compositions published in [1]. The format is intended to
be human-readable, easily compressed, and includes a canonicalization algorithm
to ensure that there is only a single representation for a glycan structure.

:title-reference:`GlycoCT{condensed}` can represent glycan structures with ambiguous
or repeating sub-units. The specification includes additional section directives with
support for stochastic sub-units as well as disjoint subgraphs, though these have not
been implemented in :mod:`glypy`.


References
----------
[1] Herget, S., Ranzinger, R., Maass, K., & Lieth, C.-W. V. D. (2008).
    GlycoCT-a unifying sequence format for carbohydrates.
    Carbohydrate Research, 343(12), 2162–2171.
    https://doi.org/10.1016/j.carres.2008.03.011

'''

import re
import warnings
from collections import defaultdict, Counter, deque, namedtuple, OrderedDict
from functools import cmp_to_key

try:
    from collections.abc import Iterator
except ImportError:
    from collections import Iterator

from glypy.utils import (
    opener, StringIO, root as rootp, tree as treep,
    make_counter, invert_dict, uid,
    RootProtocolNotSupportedError)
from glypy.utils.multimap import OrderedMultiMap
from glypy.structure import monosaccharide, substituent, glycan, Modification, constants, UnknownPosition, NoPosition
from glypy.structure.link import Link, AmbiguousLink
from .format_constants_map import (anomer_map, superclass_map,
                                   link_replacement_composition_map,
                                   modification_map, linkage_type_map)
from .file_utils import ParserError
from .tree_builder_utils import (
    decorate_tree,
    undecorate_tree,
    find_root,
    try_int,
    StructurePrecisionEnum,
    AbstractGraphEntryEnum, NodeCollection)
from glypy.composition import Composition

try:
    range = xrange
except NameError:
    pass


__id = id

Glycan = glycan.Glycan
Monosaccharide = monosaccharide.Monosaccharide
Substituent = substituent.Substituent

Configuration = constants.Configuration
Stem = constants.Stem


START = "!START"
REPINNER = "!REPINNER"
UNDINNER = "!UNDINNER"
RES = "RES"
LIN = "LIN"
REP = "REP"
ALT = "ALT"
UND = "UND"
ISO = "ISO"
NON = "NON"

TERMINAL_STATES = {
    RES,
    LIN,
    ISO,
    NON
}

subsituent_start = "s"
base_start = "b"
repeat_start = "r"
alternative_start = "a"

#: Pattern for parsing the lines of the RES section corresponding
#: to individual |Monosaccharide| residues
res_pattern = re.compile(
    r'''
    (?P<anomer>[abxo])?
    (?P<conf_stem>(?:-[dlx][a-z]+)+)?-?
    (?P<superclass>[A-Z]+)-?
    (?P<indices>[0-9x]+:[0-9x]+)
    (?P<modifications>(\|[0-9x,]+:[0-9a-z]+)+)?
    ''', re.VERBOSE)

#: Pattern for parsing the potentially repeated |Configuration| and |Stem|
#: regions of the lines of the RES section.
conf_stem_pattern = re.compile(r'(?P<config>[dlx])(?P<stem>[a-z]+)')

#: Pattern for parsing modifications found on monosaccharide residue
#: lines in the RES section
modification_pattern = re.compile(r"\|?([0-9,x]+):([^\|;\n]+)")


#: Pattern for parsing |Link| lines found in the LIN section
link_pattern = re.compile(
    r'''(?P<doc_index>\d+)?:
    (?P<parent_residue_index>\d+)
    (?P<parent_atom_replaced>[odhnx])
    \((?P<parent_attachment_position>-?[0-9\-\|]+)[\+\-]
        (?P<child_attachment_position>-?[0-9\-\|]+)\)
    (?P<child_residue_index>\d+)
    (?P<child_atom_replaced>[odhnx])
        ''', re.VERBOSE)


#: Special truncation of the :data:`link_pattern` which is used on
#: REP header sections
internal_link_pattern = re.compile(
    r'''(?P<parent_residue_index>\d+)
    (?P<parent_atom_replaced>[odhnx])
    \((?P<parent_attachment_position>-?[0-9\-\|]+)[\+\-]
        (?P<child_attachment_position>-?[0-9\-\|]+)\)
    (?P<child_residue_index>\d+)
    (?P<child_atom_replaced>[odhnx])
    ''',
    re.VERBOSE)

#: Pattern for interpreting the REP# instance header section
rep_header_pattern = re.compile(
    r'''REP(?P<repeat_index>\d+):
    (?P<internal_linkage>.+)
    =(?P<lower_multitude>-?\d+)-(?P<higher_multitude>-?\d+)''', re.VERBOSE)

repeat_line_pattern = re.compile(r"^(?P<graph_index>\d+)r:r(?P<repeat_index>\d+)")

und_header_pattern = re.compile(r'''UND(?P<und_index>\d+):
    (?P<major>\d+(\.\d*)?):
    (?P<minor>\d+(\.\d*)?)
    ''', re.VERBOSE)

und_link_pattern = re.compile(r'''
    (?P<parent_atom_replaced>[odhnx])
    \((?P<parent_attachment_position>-?[0-9\-\|]+)[\+\-]
        (?P<child_attachment_position>-?[0-9\-\|]+)\)
    (?P<child_atom_replaced>[odhnx])
    ''', re.VERBOSE)


class GlycoCTError(ParserError):
    """Base error for GlycoCT-based parsing exceptions.
    """
    pass


class GlycoCTSectionUnsupported(GlycoCTError):
    '''Indicates that the GlycoCT parser has encountered a section that
    it does not know how to parse.
    '''
    pass


class DeferredRetrieval(object):
    """Callback object to invoke a :class:`GlycoCTGraph` instance's
    :meth:`get_node` method with a set of stored parameters
    at a later time.

    Attributes
    ----------
    direction : AbstractGraphDirectionEnum
        The direction with which to retieve the node
    graph : GlycoCTGraph
        The node graph to call :meth:`GlycoCTGraph.get_node` with
    id : object
        The id of the node to retrieve
    """
    def __init__(self, graph, id, direction=None):
        self.graph = graph
        self.id = id
        self.direction = direction

    def __call__(self):
        return self.graph.get_node(self.id, self.direction)

    def __repr__(self):
        return "{s.__class__.__name__}({s.graph}, {s.id}, {s.direction})".format(
            s=self)


class GlycoCTGraph(object):
    """A graph to store nodes from parsing GlycoCT
    text in.

    Implements a Mapping interface

    Attributes
    ----------
    graph : dict
        Mapping from node id to node-like objects
    id: tuple
        A pair of (class name, ID)
    """
    def __init__(self, graph=None):
        if graph is None:
            graph = dict()
        self.graph = graph
        self.id = (self.__class__, uid())

    def __repr__(self):
        return "{self.__class__.__name__}({self.graph})".format(self=self)

    def __contains__(self, k):
        return k in self.graph

    def __getitem__(self, k):
        return self.get_node(k)

    def __setitem__(self, k, v):
        self.put_node(k, v)

    def __len__(self):
        return len(self.graph)

    def keys(self):
        '''See :meth:`~.Mapping.keys`
        '''
        return self.graph.keys()

    def values(self):
        '''See :meth:`~.Mapping.values`
        '''
        return self.graph.values()

    def items(self):
        '''See :meth:`~.Mapping.items`
        '''
        return self.graph.items()

    def clear(self):
        '''See :meth:`~.MutableMapping.clear`
        '''
        self.graph.clear()

    def __iter__(self):
        return iter(self.graph)

    def get_node(self, id, direction=None):
        """Get a node by its id value.

        Parameters
        ----------
        id : object
            The node's id. Will be case as an :class:`int`
        direction : object, optional
            Included for compatibility, ignored.

        Returns
        -------
        :class:`~.Monosaccharide` or :class:`~.Substituent
        """
        id = int(id)
        return self.graph[id]

    def put_node(self, id, value):
        """Store a node for the given id value.

        Parameters
        ----------
        id : object
            The node's id. Will be case as an :class:`int`
        value : :class:`~.Monosaccharide` or :class:`~.Substituent
            The node to store
        """
        id = int(id)
        self.graph[id] = value

    def form_link(self, parent, child, parent_position, child_position, parent_loss,
                  child_loss, parent_linkage_type=None, child_linkage_type=None, id=None):
        """Form a :class:`~.Link` between `parent` and `child` with the specified parameters.

        If more than one position is passed for `parent_position` or `child_position`, an
        :class:`~.AmbiguousLink` will be created instead.

        Parameters
        ----------
        parent : :class:`~.Monosaccharide` or :class:`~.Substituent
            The parent node in the bond
        child : :class:`~.Monosaccharide` or :class:`~.Substituent
            The child node in the bond
        parent_position : list
            The set of possible positions on the parent node to attach to.
        child_position : list
            The set of possible positions on the child node to attach to.
        parent_loss : str or :class:`~.Compositition`
            The composition lost from the parent node
        child_loss : str or :class:`~.Compositition`
            The composition lost from the child node
        parent_linkage_type : :class:`~.EnumValue`, optional
            A :class:`~.LinkageType` entry describing how the linkage is formed on the parent
        child_linkage_type : :class:`~.EnumValue`, optional
            A :class:`~.LinkageType` entry describing how the linkage is formed on the child
        id : object, optional
            The within-graph unique identifier of the :class:`~.Link` object

        Returns
        -------
        :class:`~.Link`
        """
        if parent.node_type is Substituent.node_type and\
                child.node_type is Monosaccharide.node_type:
            warnings.warn(
                "A monosaccharide with a substituent parent has been detected. "
                "These structures are not fully supported and may not traverse as expected "
                "by default.", stacklevel=7)

        if len(parent_position) > 1 or len(child_position) > 1:
            link_obj = AmbiguousLink(
                parent, child, parent_position=list(map(int, parent_position)),
                child_position=list(map(int, child_position)), parent_loss=parent_loss,
                child_loss=child_loss, id=id, parent_linkage_type=parent_linkage_type,
                child_linkage_type=child_linkage_type)
            link_obj.find_open_position()
        else:
            link_obj = Link(
                parent, child, parent_position=int(parent_position[0]),
                child_position=int(child_position[0]), parent_loss=parent_loss,
                child_loss=child_loss, parent_linkage_type=parent_linkage_type,
                child_linkage_type=child_linkage_type)
        return link_obj

    def deferred_retrieval(self, id, direction=None):
        """Construct a :class:`DeferredRetrieval` instance to carry out the
        :meth:`get_node` at a later time.

        Parameters
        ----------
        id : object
            The node's id. Will be case as an :class:`int`
        direction : object, optional
            Included for compatibility, ignored.

        Returns
        -------
        :class:`DeferredRetrieval`
        """
        return DeferredRetrieval(self, id, direction)

    def __root__(self):
        return self.find_root_nodes()[0]

    def find_root_nodes(self):
        """Find "root" nodes within the graph.

        Returns
        -------
        list
        """
        roots = []
        for _index, node in self.items():
            try:
                if node.parents():
                    continue
            except AttributeError:
                if not isinstance(node, GlycoCTGraph):
                    raise
                else:
                    if rootp(node).parents():
                        continue
            roots.append(node)
        if not roots:
            roots.append(sorted(self.items())[0][1])
        return roots

    def visit(self, node, visited=None, fn=None):
        """Visit `node`, calling `fn` on it, and then call :meth:`visit`
        on each connected node from `node` that had not previously been
        visited (tracked in `visited`).

        If a `node` is actually a :class:`GlycoCTGraph`, it will be
        traversed.

        Parameters
        ----------
        node : :class:`~.Monosaccharide` or :class:`~.Substituent`
            The node to visit.
        visited : :class:`set`, optional
            The set of previously visited node ids. If not provided, an empty set
            will be used.
        fn : :class:`callable`, optional
            The function to call on each node.

        Returns
        -------
        :class:`set`:
            The visited nodes.
        """
        if visited is None:
            visited = set()
        if isinstance(node, GlycoCTGraph):
            for node in node.find_root_nodes():
                self.visit(node, visited, fn)
            return visited

        visited.add(node.id)
        if fn is not None:
            fn(node, visited)
        for _position, link in node.links.items():
            ref = link.to(node)
            if ref.id in visited:
                continue
            else:
                self.visit(ref, visited, fn)
        try:
            for _position, link in node.substituent_links.items():
                ref = link.to(node)
                if ref.id in visited:
                    continue
                else:
                    self.visit(ref, visited, fn)
        except AttributeError:
            pass
        return visited

    def is_fully_connected(self):
        """Check that the graph is fully connected, meaning it has only
        one root node.

        Returns
        -------
        bool
        """
        roots = self.find_root_nodes()
        visited = self.visit(roots[0])
        return len(visited) >= len(self)


class GlycoCTGraphStack(GlycoCTGraph):
    """Represent a stack of :class:`GlycoCTGraph` instances,
    which may be nested inside another graph.

    Attributes
    ----------
    stack: list
        The stack of :class:`GlycoCTGraph` instances in the current state
    history: list
        The historical sequence of :class:`GlycoCTGraph` instances, added to
        but never removed.
    """
    def __init__(self, stack=None, parent=None):  # pylint: disable=super-init-not-called
        if stack is None:
            stack = deque([GlycoCTSubgraph(parent=parent)])
        else:  # pragma: no cover
            _stack = deque()
            _parent = parent
            for level in stack:
                _stack.append(GlycoCTSubgraph(level, parent=_parent))
                _parent = _stack[-1]
            stack = _stack
        self.stack = stack
        self.history = list(stack)
        self.id = (self.__class__, uid())

    @property
    def graph(self):
        """The top :class:`GlycoCTGraph` on the stack

        Returns
        -------
        :class:`GlycoCTGraph`
        """
        return self.stack[-1]

    @property
    def parent(self):
        """The graph that contains this one.

        Returns
        -------
        :class:`GlycoCTGraph`
        """
        return self.stack[0].parent

    @parent.setter
    def parent(self, value):
        self.stack[0].parent = value

    def get_node(self, id, direction=None):
        for level in reversed(self.history):
            try:
                return level.get_node(id, direction)
            except KeyError:
                continue
        raise KeyError(id)

    def deferred_retrieval(self, id, direction=None):
        return DeferredRetrieval(self, id, direction)

    def add(self, subgraph, parent=None):  # pragma: no cover
        """Add `subgraph` to :attr:`stack`, setting the
        graph's parent.

        Parameters
        ----------
        subgraph: :class:`GlycoCTGraph`
            The graph to add to the stack
        parent: :class:`GlycoCTGraph`, optional
            The parent to set for `subgraph`. If :const:`None`, this
            defaults to `self`.

        """
        if subgraph.parent is None:
            if parent is None:
                parent = self
            subgraph.parent = parent
        self.push_level(subgraph)

    def push_level(self, subgraph):
        """Push `subgraph` onto the stack, and the history.

        Parameters
        ----------
        subgraph : :class:`GlycoCTGraph`
            The graph to add.
        """
        self.stack.append(subgraph)
        self.history.append(subgraph)

    def pop_level(self):
        """Pop the last item from :attr:`stack`

        Returns
        -------
        :class:`GlycoCTGraph`
        """
        return self.stack.pop()

    def clear(self):
        self.stack = deque([GlycoCTSubgraph(parent=None)])
        self.history = list(self.stack)

    def find_subgraph_containing(self, id):
        """Find the first subgraph which contains a node with
        the query `id`

        Parameters
        ----------
        id : tuple
            The node id to find.

        Returns
        -------
        :class:`GlycoCTGraph`

        Raises
        ------
        KeyError:
            If the `id` is not found.
        """
        for level in reversed(self.history):
            if id in level:
                return level
        raise KeyError(id)

    def __repr__(self):
        return "{self.__class__.__name__}({self.history})".format(self=self)

    def _count_subgraph_nodes(self):
        layer_relations = defaultdict(list)
        root_layer = self.stack[0]
        for layer in self.stack:
            if layer.parent is None:
                continue
            else:
                key = layer.parent.id
            layer_relations[key].append(layer)

        # repeated subgraphs add to the node count on a layer, but
        # the root layer is also the parent of all other subgraphs
        # which do not have nodes. This is a hacky solution to
        # count just those with nodes in the root layer.
        first_children = layer_relations.pop(root_layer.id, [])
        adjustment = 0
        for child in first_children:
            if isinstance(child, RepeatedGlycoCTSubgraph):
                adjustment += 1
        return sum(map(len, layer_relations.values())) + adjustment

    def __len__(self):
        return sum(map(len, self.stack)) - self._count_subgraph_nodes()

    def find_root_nodes(self):
        return self.stack[0].find_root_nodes()

    def is_fully_connected(self):
        roots = self.find_root_nodes()
        visited = self.visit(roots[0])
        return len(visited) >= len(self)


class GlycoCTSubgraph(GlycoCTGraph):
    """A :class:`GlycoCTGraph` which has an :attr:`parent` for use with
    :class:`GlycoCTGraphStack`.

    Attributes
    ----------
    parent: :class:`GlycoCTGraph`
        The graph that contains this one.

    """
    def __init__(self, graph=None, parent=None):
        super(GlycoCTSubgraph, self).__init__(graph)
        assert not isinstance(parent, dict)
        self.parent = parent

    def postprocess(self):
        """Apply arbitrary post-processing before finalizing the subgraph.

        To be overriden by subclasses. Base implementation is a no-op.
        """
        pass


RepeatedMultitude = namedtuple("RepeatedMultitude", "lower upper")


class RepeatedGlycoCTSubgraph(GlycoCTSubgraph):
    """Implements the machinery for representing a repeated subgraph in GlycoCT.

    Attributes
    ----------
    graph_index: int
    repeast_index: int
        The ``i``th repeating subgraph in the graph.
    internal_linkage: object
        The linkage connecting two repetitions of the subgraph
    external_linkage: object
        The linkage connecting from the final repetition and the
        outside nodes.
    multitude: :class:`RepeatedMultitude`
        Holds the lower and upper range of multiplicities this subgraph may be repeated
        to.
    repetitions: :class:`~.OrderedDict`
        The repetitions of this subgraph, materialized during :meth:`postprocess`
    postponed: :class:`~.deque`
        A queue of post-processing callbacks.

    """
    def __init__(self, graph_index, repeat_index, internal_linkage=None,
                 external_linkage=None, multitude=None, graph=None,
                 parent=None):
        super(RepeatedGlycoCTSubgraph, self).__init__(graph, parent)

        if multitude is None:
            multitude = RepeatedMultitude(-1, -1)
        else:
            multitude = RepeatedMultitude(*multitude)

        self.graph_index = graph_index
        self.repeat_index = repeat_index
        self.internal_linkage = internal_linkage
        self.external_linkage = external_linkage
        self.multitude = multitude
        self.original_graph = None
        self.repetitions = OrderedDict()
        self.postponed = deque()

    def __repr__(self):
        rep = super(RepeatedGlycoCTSubgraph, self).__repr__()
        rep = "%s, %d)" % (rep[:-1], self.repeat_index)
        return rep

    @property
    def terminal_node_index(self):
        """The index of the node that will connect to either the
        external or internal target node.

        Returns
        -------
        int
        """
        return int(self.external_linkage['parent_residue_index'])

    @property
    def last_repeat_index(self):
        """The last repeat's index in :attr:`repetitions`

        Returns
        -------
        int
        """
        sub_unit_indices = sorted(self.repetitions.keys())
        terminal_unit_ix = sub_unit_indices[-1]
        return terminal_unit_ix

    @property
    def terminal_node(self):
        """Retrieves the terminal residue from the subgraph, where
        outgoing connections start from

        Returns
        -------
        MoleculeBase
        """
        terminal_unit_ix = self.last_repeat_index
        parent_graph = self.repetitions[terminal_unit_ix]
        parent = parent_graph.get((terminal_unit_ix, self.terminal_node_index))
        return parent

    @property
    def origin_node_index(self):
        return int(self.external_linkage['child_residue_index'])

    @property
    def first_repeat_index(self):
        """The first repeat's index in :attr:`repetitions`

        Returns
        -------
        int
        """
        sub_unit_indices = sorted(self.repetitions.keys())
        return sub_unit_indices[0]

    @property
    def origin_node(self):
        """Retrieves the root residue from the subgraph, where
        the incoming external connection ends at.

        Returns
        -------
        MoleculeBase
        """
        root_unit_ix = self.first_repeat_index
        child_graph = self.repetitions[root_unit_ix]
        try:
            child = child_graph.get((root_unit_ix, self.origin_node_index))
        except IndexError:
            # the root of this subgraph is nested, coming from an expanded
            # inner subgraph. The index labels recorded in the section header
            # are not accurate anymore. Attempt to recover by returning the
            # root node of this subgraph
            if isinstance(child_graph.root.id[1], tuple):
                return child_graph.root
            else:
                raise

        return child

    def postpone(self, f, args):
        """Queue a callable `f` with `args`
        to be run during postprocessing.

        Parameters
        ----------
        f : :class:`Callable`
            Some callable object
        args : :class:`Iterable`
            The arguments to call `f` with
        """
        self.postponed.append((f, args))

    def complete_postponed_tasks(self):
        """Call all of the postponed tasks.
        """
        i = 0
        while self.postponed:
            i += 1
            f, args = self.postponed.popleft()
            f(*args)

    def __contains__(self, k):
        if self.original_graph is not None:
            return k in self.original_graph
        else:
            return k in self.graph

    @property
    def structure_precision(self):
        if -1 in self.multitude:
            return StructurePrecisionEnum.unknown
        elif self.multitude.lower == self.multitude.upper:
            return StructurePrecisionEnum.exact
        return StructurePrecisionEnum.ranging

    def _build_minigraph(self, graph):
        sub_unit_indices = sorted(map(try_int, graph.keys()))
        for key, node in list(graph.items()):
            if isinstance(node, GlycoCTGraph):
                node.postprocess()
                graph[key] = node.origin_node

        try:
            root = find_root(graph[sub_unit_indices[0]])
        except RootProtocolNotSupportedError:  # pragma: no cover
            raise GlycoCTError("Could not locate subgraph root")
        # glycan_graph = Glycan(root, index_method=None).clone()
        glycan_graph = NodeCollection.from_node(root).clone()
        return glycan_graph

    def _duplicate_graph(self, graph):
        glycan_graph = self._build_minigraph(graph)
        duplicate_graph = {}
        for k, v in graph.items():
            if isinstance(v, GlycoCTSubgraph):
                duplicate_graph[k] = v
            else:
                try:
                    duplicate_graph[k] = glycan_graph.get(v.id)
                except IndexError:
                    # if the subgraph intersperses repeats and
                    # residues, they won't be connected in the
                    # Glycan object, so we must read them back
                    # from the original graph
                    duplicate_graph[k] = graph[k]
        return duplicate_graph

    def postprocess(self, n=None):  # pylint: disable=arguments-differ
        if n is None:
            if self.multitude.upper != -1:
                n = self.multitude.upper
            elif self.multitude.lower != -1:
                n = self.multitude.lower
            else:
                n = 1

        if self.structure_precision is not StructurePrecisionEnum.unknown:
            if not (self.multitude.lower <= n <= self.multitude.upper):  # pragma: no cover
                raise GlycoCTError("{} is not within the range of {}".format(n, self.multitude))

        self.original_graph = self._duplicate_graph(self.graph)
        glycan_graph = self._build_minigraph(self.graph)

        graph = OrderedDict({1: glycan_graph.clone(index_method=None)})
        decorate_tree(graph[1], 1)

        parent_residue_index = self.terminal_node_index
        parent_atom_replaced = link_replacement_composition_map[self.internal_linkage["parent_atom_replaced"]]
        parent_linkage_type = linkage_type_map[self.internal_linkage["parent_atom_replaced"]]
        parent_attachment_position = self.internal_linkage["parent_attachment_position"]

        child_residue_index = self.origin_node_index
        child_atom_replaced = link_replacement_composition_map[self.internal_linkage["child_atom_replaced"]]
        child_linkage_type = linkage_type_map[self.internal_linkage["child_atom_replaced"]]
        child_attachment_position = self.internal_linkage["child_attachment_position"]

        op_stack = []

        for i in range(2, n + 1):
            graph[i] = glycan_graph.clone(index_method=None)
            graph[i] = decorate_tree(graph[i], i)
            parent_graph = graph[i - 1]
            child_graph = graph[i]
            parent_node = parent_graph.get((i - 1, parent_residue_index))
            child_node = child_graph.get((i, child_residue_index))
            op_stack.append(
                (self.form_link, [parent_node, child_node],
                 dict(parent_position=parent_attachment_position,
                      child_position=child_attachment_position,
                      parent_loss=parent_atom_replaced,
                      child_loss=child_atom_replaced,
                      parent_linkage_type=parent_linkage_type,
                      child_linkage_type=child_linkage_type)))

        for op in op_stack:
            f, args, kwargs = op
            f(*args, **kwargs)

        self.repetitions = graph
        self.complete_postponed_tasks()

    def handle_incoming_link(self, parent_getter, parent_position, parent_loss,
                             child_position, child_loss, id=None):
        child = self.origin_node
        parent = parent_getter()
        if parent_loss == Composition("H"):
            child_loss = Composition("OH")

        self.form_link(
            parent, child, parent_position=parent_position, child_position=child_position,
            parent_loss=parent_loss, child_loss=child_loss, id=id)

    def handle_outgoing_link(self, child_getter, parent_position, parent_loss,
                             child_position, child_loss, id=None):
        child = child_getter()
        parent = self.terminal_node
        if isinstance(child, RepeatedGlycoCTSubgraph):
            child.handle_incoming_link(
                lambda: parent, parent_position=parent_position,
                child_position=child_position, parent_loss=parent_loss, child_loss=child_loss, id=id)
        else:
            self.form_link(
                parent, child, parent_position=parent_position, child_position=child_position,
                parent_loss=parent_loss, child_loss=child_loss, id=id)

    def handle_abstract_subgraph_link(self, parent_getter, child_getter, parent_position,
                                      parent_loss, child_position, child_loss, id=None):
        parent = parent_getter()
        child = child_getter()
        self.form_link(
            parent, child, parent_position=parent_position,
            child_position=child_position, parent_loss=parent_loss,
            child_loss=child_loss, id=id)

    def get_node(self, id, direction=None):
        id = int(id)
        if direction == AbstractGraphEntryEnum.parent:
            child_graph = self.repetitions[self.first_repeat_index]
            child = child_graph.get((self.first_repeat_index, id))
            return child
        elif direction == AbstractGraphEntryEnum.child:
            parent_graph = self.repetitions[self.last_repeat_index]
            parent = parent_graph.get((self.last_repeat_index, id))
            return parent
        elif direction == AbstractGraphEntryEnum.internal or direction is None:
            return self.graph[id]
        else:
            raise Exception("Unknown direction %s" % direction)

    def __root__(self):  # pragma: no cover
        root_node = find_root(self.repetitions[self.first_repeat_index])
        return root_node

    def find_root_nodes(self):
        if not self.repetitions:
            return super(RepeatedGlycoCTSubgraph, self).find_root_nodes()
        roots = []
        for _rep_index, node_set in self.repetitions.items():
            for node in node_set:
                try:
                    if node.parents():
                        continue
                except AttributeError:
                    if not isinstance(node, GlycoCTGraph):
                        raise
                    else:
                        if rootp(node).parents():
                            continue
                roots.append(node)
        if not roots:
            roots.append(sorted(self.items())[0][1])
        return roots

    def prepare_glycan(self):
        subglycan_start = self.repetitions[self.first_repeat_index]
        subglycan_start.deindex()
        return Glycan(subglycan_start.root, index_method=None)


UndeterminedProbability = namedtuple("UndeterminedProbability", "major minor")


LinkageSpecification = namedtuple(
    "LinkageSpecification", [
        "id", "parent_residue_index", "parent_atom_replaced", "parent_attachment_position",
        "child_residue_index", "child_atom_replaced", "child_attachment_position",
        "parent_linkage_type", "child_linkage_type"
    ])


class UndeterminedGlycoCTSubgraph(GlycoCTSubgraph):
    def __init__(self, und_index, probability=None, parent_ids=None,
                 subtree_linkages=None, graph=None, parent=None):
        super(UndeterminedGlycoCTSubgraph, self).__init__(graph, parent)
        if probability is None:
            probability = UndeterminedProbability(100., 100.)
        if parent_ids is None:
            parent_ids = []
        if subtree_linkages is None:
            subtree_linkages = []
        self.parent_ids = parent_ids
        self.subtree_linkages = subtree_linkages

    def __root__(self):
        return self[sorted(self.keys())[0]]

    def postprocess(self):
        parents = [
            self.parent.get_node(parent_id) for parent_id in self.parent_ids]
        try:
            child = [rootp(self)]
        except RootProtocolNotSupportedError:  # pragma: no cover
            raise GlycoCTError("Could not locate subgraph root")

        linkage = self.subtree_linkages[0]
        parent_loss = linkage["parent_atom_replaced"]
        parent_position = linkage["parent_attachment_position"]
        child_loss = linkage["child_atom_replaced"]
        child_position = linkage["child_attachment_position"]

        link_obj = AmbiguousLink(
            parents, child, parent_position=list(map(int, parent_position)),
            child_position=list(map(int, child_position)), parent_loss=parent_loss,
            child_loss=child_loss)
        try:
            link_obj.find_open_position()
        except ValueError:  # pragma: no cover
            if link_obj.child_position == 1 and Modification.Acidic in link_obj.child.modifications[1]:
                link_obj.child_position = 2
                ix = link_obj.child_position_choices.index(1)
                link_obj.child_position_choices.pop(ix)
                link_obj.child_position_choices.insert(ix, 2)
                link_obj.apply()
                link_obj.find_open_position()
            else:
                raise


def _build_graph(glycoct_str):  # pragma: no cover
    rep = StringIO(glycoct_str)
    inst = GlycoCTReader(rep, completes=False)
    return next(inst)


def extract_composition(parser):
    from glypy.structure.glycan_composition import (
        GlycanComposition, MonosaccharideResidue, SubstituentResidue)
    store = GlycanComposition()
    # remove links between layers in the stack
    for layer in list(parser.stack)[1:]:
        node = rootp(layer)
        for _position, link in list(node.links.items()):
            if link.is_child(node):
                link.break_link(refund=True)

    for layer in parser.stack:
        for node in layer.find_root_nodes():
            if isinstance(node, Monosaccharide):
                node = MonosaccharideResidue.from_monosaccharide(node)
            elif isinstance(node, Substituent):
                node = SubstituentResidue(node.name)
            else:
                raise ValueError(node)
            store[node] += 1
    return store


class GlycoCTReader(GlycoCTGraphStack, Iterator):
    """Parse :title-reference:`GlycoCT{condensed}` text data into |Glycan| objects.

    The parser implements the :class:`Iterator` interface, yielding successive glycans
    from a text stream separated by empty lines.

    The parser can understand fully specified and partially ambiguous structures.
    When :attr:`allow_repeats` is |True| and a ``REP`` section is encountered, it
    will be expanded to its minimum multiplicity, or 1 if the minimum is unknown.
    ``UND`` sections will be connected to the main graph by :class:`~.AmbiguousLink`
    instead of :class:`~.Link` objects.

    Attributes
    ----------
    allow_repeats : :class:`bool`
        Whether or not to permit ``REP`` sections. Defaults to |True|
    completes : :class:`bool`
        Whether or not to translate the built graph into a |Glycan| object. Defaults
        to |True|
    handle : file-like
        The text file being read from
    in_repeat : :class:`bool`
        Indicates the parser is currently parsing a ``REP`` section's sub-graph
    in_undetermined : bool
        Indicates the parser is currently parsing a ``UND`` section's sub-graph
    postponed : list
        Holds all the deferred operations for the top-most graph as :class:`callable`
        objects
    root : :class:`Monosaccharide`
        The root node of the produced graph
    state : str
        The current state of the parser's state machine
    structure_class : type
        The |Glycan| sub-class to produce
    repeats : dict
        Maps RES section index to :class:`RepeatedGlycoCTSubgraph`
    undetermineds : dict
        Maps UND section index to :class:`UndeterminedGlycoCTSubgraph`
    """

    @classmethod
    def loads(cls, glycoct_str, structure_class=Glycan, allow_repeats=True):
        '''Parse results from |str|'''
        rep = StringIO(glycoct_str)
        return cls(rep, structure_class=structure_class, allow_repeats=allow_repeats)

    def __init__(self, stream, structure_class=Glycan, allow_repeats=True, completes=True):
        super(GlycoCTReader, self).__init__()

        self._state = None
        self.state = START

        self.completes = completes

        self.handle = opener(stream, "r")
        self.in_repeat = False
        self.in_undetermined = False
        self.repeats = {}
        self.undetermineds = {}
        self.postponed = []

        self.root = None
        self._iter = None

        self.allow_repeats = allow_repeats
        self.structure_class = structure_class

        self._index = 0
        self._source_line = 0
        self._segment_iterator = None

        self._output_queue = deque()

    def _read(self):
        for line in self.handle:
            self._source_line += 1
            for segment in re.split(r"\s|;", line):
                if not segment.strip():
                    continue
                self._current_segment = segment
                yield segment

    def _reset(self):
        self.clear()
        self.root = None
        self.postponed = []
        self.repeats.clear()
        self.undetermineds.clear()
        self.in_repeat = False
        self._index += 1

    def reset(self):
        if self.completes:
            self._reset()

    def __iter__(self):
        '''
        Calls :meth:`parse` and stores it for reuse with :meth:`__next__`
        '''
        if self._iter is None:
            self._iter = self.parse()
        return self._iter

    def next(self):
        '''
        Calls :meth:`parse` if the internal iterator has not been instantiated
        '''
        if self._iter is None:
            iter(self)
        return next(self._iter)

    #: Alias for next. Supports Py3 Iterator interface
    __next__ = next

    @property
    def state(self):
        return self._state

    @state.setter
    def state(self, value):
        self._state = value

    def _parse_modifications(self, residue_dict):
        mods = residue_dict["modifications"]
        modifications = OrderedMultiMap()
        if mods is not None:
            for modp, mod in modification_pattern.findall(mods):
                positions = modp.split(",")
                if len(positions) > 1:
                    warnings.warn("Multi-site Modifications are not fully supported")
                for p in positions:
                    modifications[try_int(p)] = modification_map[mod]
        is_reduced = "aldi" in modifications[1]
        if is_reduced:
            modifications.pop(1, "aldi")
            is_reduced = monosaccharide.ReducedEnd()
        else:
            is_reduced = None
        return modifications, is_reduced

    def _parse_conf_stem(self, residue_dict):
        conf_stem = residue_dict["conf_stem"]
        if conf_stem is not None:
            config, stem = zip(*conf_stem_pattern.findall(conf_stem))
        else:
            config = ('x',)
            stem = ('x',)
        stem_ = tuple(Stem[s] for s in stem)
        configuration_ = tuple(Configuration[c] for c in config)
        return stem_, configuration_

    def handle_residue_line(self, line):
        '''
        Handle a base line, creates an instance of |Monosaccharide|
        and adds it to :attr:`graph` at the given index.

        Called by :meth:`parse`
        '''
        _, ix, residue_str = re.split(r"^(\d+)b", line, maxsplit=1)
        residue_dict = res_pattern.search(residue_str).groupdict()

        modifications, is_reduced = self._parse_modifications(residue_dict)

        stem_, configuration_ = self._parse_conf_stem(residue_dict)

        ring_start_, ring_end_ = [
            (try_int(i) if i != 'x' else UnknownPosition) for i in residue_dict["indices"].split(":")]

        anomer_ = anomer_map[residue_dict['anomer']]
        super_class_ = superclass_map[residue_dict['superclass']]
        ix = int(ix)
        residue = monosaccharide.Monosaccharide(
            fast=True, stem=stem_, modifications=modifications,
            reduced=is_reduced, configuration=configuration_,
            ring_start=ring_start_, ring_end=ring_end_, anomer=anomer_,
            superclass=super_class_, id=ix)

        self.put_node(ix, residue)
        if self.root is None:
            self.root = residue

    def handle_residue_substituent(self, line):
        '''
        Handle a substituent line, creates an instance of |Substituent|
        and adds it to :attr:`graph` at the given index. The |Substituent| object is not yet linked
        to a |Monosaccharide| instance.

        Called by :meth:`parse`

        '''
        _, ix, subsituent_str = re.split(r"^(\d+)s:", line, maxsplit=1)
        sub = Substituent(subsituent_str.strip())

        self[int(ix)] = sub

    def handle_blank(self):
        self._complete_structure()

    def enter_res(self):
        if self.state in (START, REPINNER, UNDINNER):
            pass
        elif self.state in TERMINAL_STATES:
            self.in_repeat = False
            self._complete_structure()
        else:
            raise GlycoCTError("Invalid State Transition at line %d" % self._source_line)

        self.state = RES

    def enter_lin(self):
        if self.state != RES:
            raise GlycoCTError("LIN before RES at line %d" % self._source_line)
        self.state = LIN

    def enter_rep(self):
        if not self.allow_repeats:
            raise GlycoCTSectionUnsupported(
                "Repeat are not allowed (set allow_repeats=True to allow them) at line %d" % self._source_line)
        self.state = REP
        self.in_repeat = True

    def enter_und(self):
        self.state = UND
        self.in_undetermined = True

    def parse_link(self, line):
        link_dict = link_pattern.search(line)
        if link_dict is not None:
            link_dict = link_dict.groupdict()
        else:
            raise GlycoCTError("Could not interpret link", line)
        id = link_dict['doc_index']
        parent_residue_index = int(link_dict['parent_residue_index'])
        child_residue_index = int(link_dict['child_residue_index'])

        parent_atom_replaced = link_replacement_composition_map[link_dict["parent_atom_replaced"]]
        parent_attachment_position = list(map(int, link_dict["parent_attachment_position"].split("|")))
        try:
            parent_linkage_type = linkage_type_map[link_dict['parent_atom_replaced']]
        except KeyError:
            parent_linkage_type = constants.LinkageType.x

        child_atom_replaced = link_replacement_composition_map[link_dict["child_atom_replaced"]]
        child_attachment_position = list(map(int, link_dict["child_attachment_position"].split("|")))
        try:
            child_linkage_type = linkage_type_map[link_dict['child_atom_replaced']]
        except KeyError:
            child_linkage_type = constants.LinkageType.x

        return LinkageSpecification(
            id, parent_residue_index, parent_atom_replaced, parent_attachment_position,
            child_residue_index, child_atom_replaced, child_attachment_position,
            parent_linkage_type, child_linkage_type)

    def handle_linkage(self, line):
        '''
        Handle a linkage line, creates an instance of |Link| and
        attaches it to the two referenced nodes in :attr:`graph`. The parent node is always
        an instance of |Monosaccharide|, and the child node
        may either be an instance of |Monosaccharide| or
        |Substituent| or |Monosaccharide|.

        Called by :meth:`parse`

        See also |Link| for more information on the impact of instantiating
        a |Link| object.
        '''
        id, parent_residue_index, parent_atom_replaced, parent_attachment_position,\
            child_residue_index, child_atom_replaced, child_attachment_position,\
            parent_linkage_type, child_linkage_type = self.parse_link(line)

        parent = self.get_node(parent_residue_index)
        child = self.get_node(child_residue_index)

        is_parent_repeat = isinstance(parent, RepeatedGlycoCTSubgraph)
        is_child_repeat = isinstance(child, RepeatedGlycoCTSubgraph)

        if is_parent_repeat and is_child_repeat:
            inner = max([parent, child], key=lambda x: x.repeat_index)
            if child == inner:

                def child_getter():
                    return child.origin_node

                def parent_getter():
                    return parent.get_node(
                        parent.terminal_node_index,
                        AbstractGraphEntryEnum.internal)
            else:

                def child_getter():
                    return child.get_node(
                        child.origin_node_index,
                        AbstractGraphEntryEnum.internal)

                def parent_getter():
                    return parent.terminal_node
            if parent_atom_replaced == child_atom_replaced == Composition('H'):
                parent_atom_replaced = Composition('H')
                child_atom_replaced = Composition('OH')
            inner.postpone(
                inner.handle_abstract_subgraph_link,
                (parent_getter,
                 child_getter,
                 parent_attachment_position, parent_atom_replaced,
                 child_attachment_position, child_atom_replaced, id)
            )
        elif is_parent_repeat:
            parent.postpone(parent.handle_outgoing_link, (
                self.deferred_retrieval(child_residue_index),
                parent_attachment_position, parent_atom_replaced,
                child_attachment_position, child_atom_replaced, id)
            )
        elif is_child_repeat:
            child.postpone(child.handle_incoming_link, (
                self.deferred_retrieval(parent_residue_index),
                parent_attachment_position, parent_atom_replaced,
                child_attachment_position, child_atom_replaced, id)
            )
        else:
            self.form_link(
                parent, child,
                parent_position=parent_attachment_position, child_position=child_attachment_position,
                parent_loss=parent_atom_replaced, child_loss=child_atom_replaced, id=id,
                parent_linkage_type=parent_linkage_type, child_linkage_type=child_linkage_type)

    def handle_repeat_stub(self, line):
        if not self.allow_repeats:
            raise GlycoCTSectionUnsupported(
                "Repeat are not allowed (set allow_repeats=True to allow them)")

        match = repeat_line_pattern.search(line).groupdict()

        graph_index = try_int(match['graph_index'])
        repeat_index = try_int(match["repeat_index"])

        repeat = RepeatedGlycoCTSubgraph(
            int(graph_index), int(repeat_index), parent=self.graph)
        repeat._index = self._index

        self[graph_index] = repeat
        self.repeats[repeat_index] = repeat

        if self.root is None:
            self.root = repeat

    def handle_repeat_inner(self, line):
        if not self.in_repeat:
            raise GlycoCTError(
                "Encountered %r outside of REP at line %d" % (
                    line, self._source_line))
        header_dict = rep_header_pattern.search(line).groupdict()

        repeat_index = int(header_dict['repeat_index'])
        repeat_record = self.repeats[repeat_index]

        self.push_level(repeat_record)

        linkage = internal_link_pattern.search(header_dict['internal_linkage']).groupdict()
        repeat_record.internal_linkage = linkage
        repeat_record.external_linkage = linkage
        repeat_record.multitude = RepeatedMultitude(
            try_int(header_dict['lower_multitude']),
            try_int(header_dict['higher_multitude']))
        self.state = REPINNER

    def handle_und_inner(self, line):
        if not self.in_undetermined:
            raise GlycoCTError("Encountered %r outside of UND at line %d" % (
                line, self._source_line))
        header_dict = und_header_pattern.search(line).groupdict()
        parent_line = next(self._segment_iterator)
        subtree_linkage_line = next(self._segment_iterator)

        ids = list(map(int, parent_line.split(":")[1].split("|")))
        subtree_linkages = []
        match = und_link_pattern.search(subtree_linkage_line.split(":")[1])
        if match is None:
            raise GlycoCTError("Could not interpret UND SubtreeLinkage %r at line %d" % (
                subtree_linkage_line, self._source_line))
        else:
            link_dict = match.groupdict()
            link_dict["parent_atom_replaced"] = link_replacement_composition_map[
                link_dict["parent_atom_replaced"]]
            link_dict["parent_attachment_position"] = list(
                map(int, link_dict["parent_attachment_position"].split("|")))
            link_dict["child_atom_replaced"] = link_replacement_composition_map[
                link_dict["child_atom_replaced"]]
            link_dict["child_attachment_position"] = list(
                map(int, link_dict["child_attachment_position"].split("|")))
            subtree_linkages.append(link_dict)

        und_index = int(header_dict['und_index'])
        prob = UndeterminedProbability(float(header_dict['major']), float(header_dict['minor']))
        record = UndeterminedGlycoCTSubgraph(
            und_index, prob, parent_ids=ids,
            subtree_linkages=subtree_linkages, parent=self)
        self.undetermineds[und_index] = record
        self.push_level(record)
        self.state = UNDINNER

    def _complete_structure(self):
        if self.completes:
            result = self.postprocess()
            if result is not None:
                self._output_queue.append(result)
            self.reset()
        else:
            self._output_queue.append(self)

    def postprocess(self):
        '''
        Handle all deferred operations such as binding together and expanding
        repeating units. Removes any distinguishing markers on node ids, and
        constructs a new instance of :attr:`structure_class` from the accumulated
        graph

        Returns
        -------
        Glycan
        '''
        for level in reversed(self.history):
            level.postprocess()

        for postop in self.postponed:
            postop[0](*postop[1:])

        if self.root is None:
            return None

        if self.is_fully_connected():
            try:
                inst = undecorate_tree(
                    self.structure_class(
                        root=rootp(self.root), index_method=None)
                ).reindex()
                return inst
            except RootProtocolNotSupportedError:  # pragma: no cover
                raise GlycoCTError("Could not locate graph root")
        else:
            # warnings.warn("The parsed structure was not fully connected. Producing a Composition")
            return extract_composition(self)

    def parse(self):
        '''
        Returns an iterator that yields each complete :class:`Glycan` instance
        from the underlying text stream.
        '''

        # Create a reference to the segment iterator
        # as late as possible, but bind it to the object
        # state so it can be referenced independent of this
        # outermost loop
        self._segment_iterator = self._read()

        for line in self._segment_iterator:
            if line.strip() == "":
                self.handle_blank()
                while self._output_queue:
                    yield self._output_queue.popleft()

            elif RES == line.strip():
                self.enter_res()
                while self._output_queue:
                    yield self._output_queue.popleft()

            elif LIN == line.strip():
                self.enter_lin()

            elif REP == line.strip():
                self.enter_rep()

            elif UND == line.strip():
                self.enter_und()

            # REP definition block
            elif line.strip()[:3] == REP:
                self.handle_repeat_inner(line)
            elif line.strip()[:3] == UND:
                self.handle_und_inner(line)
            elif ALT == line.strip():
                raise GlycoCTSectionUnsupported(ALT)

            elif re.search(r"^(\d+)b", line) and self.state == RES:
                self.handle_residue_line(line)
            elif re.search(r"^(\d+)s:", line) and self.state == RES:
                self.handle_residue_substituent(line)
            elif re.search(r"^(\d+)r:", line) and self.state == RES:
                self.handle_repeat_stub(line)
            elif re.search(r"^(\d+):(\d+)", line) and self.state == LIN:
                self.handle_linkage(line)
            else:
                raise GlycoCTError("Unknown format error: %s on line %d" % (line, self._source_line))

        if self.root is not None:
            self._complete_structure()
            while self._output_queue:
                yield self._output_queue.popleft()


GlycoCT = GlycoCTReader


def read(stream, structure_class=Glycan, allow_repeats=True):
    '''
    A convenience wrapper for :class:`GlycoCTReader`
    '''
    return GlycoCTReader(stream, structure_class=structure_class, allow_repeats=allow_repeats)


def load(stream, structure_class=Glycan, allow_repeats=True, allow_multiple=True):  # pragma: no cover
    """Read all structures from the provided text stream.

    Parameters
    ----------
    stream : file-like
        The text stream to parse structures from
    structure_class : type, optional
        :class:`~.Glycan` subclass to use
    allow_repeats : bool, optional
        Whether or not to allow ``REP`` sections

    Returns
    -------
    :class:`~.Glycan` or :class:`list` of :class:`~.Glycan`
    """
    g = GlycoCTReader(stream, structure_class=structure_class, allow_repeats=allow_repeats)
    first = next(g)
    if not allow_multiple:
        return first
    second = None
    try:
        second = next(g)
        collection = [first, second]
        collection.extend(g)
        return collection
    except StopIteration:
        return first


def loads(text, structure_class=Glycan, allow_repeats=True, allow_multiple=True):
    """Read all structures from the provided text string.

    Parameters
    ----------
    text : str
        The text to parse structures from
    structure_class : type, optional
        :class:`~.Glycan` subclass to use
    allow_repeats : bool, optional
        Whether or not to allow ``REP`` sections

    Returns
    -------
    :class:`~.Glycan` or :class:`list` of :class:`~.Glycan`
    """

    text_buffer = StringIO(text)
    return load(text_buffer, structure_class, allow_repeats, allow_multiple)


def detect_glycoct(string):
    return string.lstrip()[:3] == "RES"


invert_anomer_map = invert_dict(anomer_map)
invert_superclass_map = invert_dict(superclass_map)


class DictTree(object):
    def __init__(self, state=START, store=None):
        if store is None:
            store = {}
        self.store = defaultdict(dict, store)
        self.state = state

    def __getitem__(self, key):
        for subtree in self.store.values():
            try:
                return subtree[key]
            except KeyError:
                continue
        return self.store[key]

    def __setitem__(self, key, value):
        self.store[self.state][key] = value

    def __len__(self):
        return sum(map(len, self.store))

    def __iter__(self):
        return iter(self.store)

    def keys(self):
        return self.store.keys()

    def items(self):
        vals = list(self.store.values())
        acc = vals[0].items()
        for v in vals[1:]:
            acc = acc + v.items()
        return acc

    def __contains__(self, key):
        return key in self.store

    def get(self, key, subtree=None, default=None):
        if subtree is None:
            subtree = self.state
        return self.store[subtree].get(key, default)


class GlycoCTWriterBase(object):
    """Summary

    Attributes
    ----------
    buffer : file-like
        The buffer to write structures to. If :attr:`nobuffer` is |True|,
        this will be a :class:`~.StringIO` object which will be returned
        on each write.
    dependencies : :class:`defaultdict` of :class:`dict`
        Track the relationships between child nodes and their parents.
        Used during linkage writing
    full : :class:`bool`
        Whether or not to traverse :class:`~.Monosaccharide`-:class:`~.Monosaccharide`
        linkages.
    index_to_residue : :class:`DictTree`
        A state-specific mapping from index to :attr:`~.Monosaccharide.id`.
    lin_accumulator : list
        Accumulator list of :class:`~.Link` objects.
    lin_counter : function
        A stateful counter which when called returns the next
        integer in a sequence used to index entries in the `LIN` section
    nobuffer : bool
        Whether or not the writer was initialized with a write-able buffer
    res_counter : function
        A stateful counter which when called returns the next
        integer in a sequence used to index entries in the `RES` section
    residue_to_index : :class:`DictTree`
        A state-specific mapping from :attr:`~.Monosaccharide.id` to index.
    state : str
        The current state of the writer
    structure : :class:`~.SaccharideCollection`
        The structure currently being written. May be a :class:`~.Monosaccharide`,
        :class:`~.Glycan`, :class:`~.GlycanComposition`.
    und_counter : function
        A stateful counter which when called returns the next
        integer in a sequence used to index `UND` sections.
    """

    def __init__(self, structure=None, buffer=None, full=True):
        self.nobuffer = False
        if buffer is None:
            self.nobuffer = True
            buffer = StringIO()

        self.buffer = buffer
        self.structure = structure
        self.full = full

        self.state = START
        self._initialize_counters()
        self._initialize_index_tree()
        self._initialize_link_trackers()

    def _initialize_counters(self):
        self.res_counter = make_counter()
        self.lin_counter = make_counter()
        self.und_counter = make_counter()

    def _initialize_index_tree(self):
        # Look-ups for mapping RES nodes to objects by section index and id,
        # respectively
        self.index_to_residue = DictTree(self.state)
        self.residue_to_index = DictTree(self.state)

    def _initialize_link_trackers(self):
        # Accumulator for linkage indices and mapping linkage indices to
        # dependent RES indices
        self.lin_accumulator = []
        self.dependencies = defaultdict(dict)

    @property
    def structure(self):
        return self._structure

    @structure.setter
    def structure(self, value):
        if value is None:
            self._structure = value
            return
        try:
            structure = treep(value)
        except TypeError:
            try:
                root = rootp(value)
                structure = Glycan(root, index_method=None)
            except TypeError:
                raise TypeError("Could not extract or construct a tree structure from %r" % value)
        self._structure = structure

    def _reset(self):
        self.state = START
        self._initialize_counters()
        self._initialize_index_tree()
        self._initialize_link_trackers()
        if self.nobuffer:
            self.buffer = StringIO()

    def _glycoct_sigils(self, link):
        '''
        Helper method for determining which GlycoCT symbols and losses to present
        '''
        parent_loss_str, child_loss_str = link._glycoct_sigils()

        return parent_loss_str, child_loss_str

    def handle_link(self, link, ix, parent_ix, child_ix):
        parent_loss_str, child_loss_str = self._glycoct_sigils(link)

        if link.has_ambiguous_linkage():
            rep = "{ix}:{parent_ix}{parent_loss}({parent_position}+{child_position}){child_ix}{child_loss}"
            return rep.format(
                ix=ix,
                parent_ix=parent_ix,
                parent_loss=parent_loss_str,
                parent_position='|'.join(map(str, link.parent_position_choices)),
                child_ix=child_ix,
                child_loss=child_loss_str,
                child_position='|'.join(map(str, link.child_position_choices)))
        else:
            rep = "{ix}:{parent_ix}{parent_loss}({parent_position}+{child_position}){child_ix}{child_loss}"
            return rep.format(
                ix=ix,
                parent_ix=parent_ix,
                parent_loss=parent_loss_str,
                parent_position=link.parent_position,
                child_ix=child_ix,
                child_loss=child_loss_str,
                child_position=link.child_position)

    def handle_substituent(self, substituent):  # pylint: disable=redefined-outer-name
        return "s:{0}".format(substituent.name.replace("_", "-"))

    def _format_monosaccharide(self, monosaccharide):  # pylint: disable=redefined-outer-name
        residue_template = "{ix}b:{anomer}{conf_stem}{superclass}-{ring_start}:{ring_end}{modifications}"

        # This index is reused many times
        monosaccharide_index = self.res_counter()

        # Format individual fields
        anomer = invert_anomer_map[monosaccharide.anomer]
        conf_stem = ''.join("-{0}{1}".format(c.name, s.name)
                            for c, s in zip(monosaccharide.configuration, monosaccharide.stem))
        if None in monosaccharide.configuration and None in monosaccharide.stem:
            conf_stem = ''
        superclass = "-" + invert_superclass_map[monosaccharide.superclass]

        modifications = '|'.join(
            "{0}:{1}".format(k, v.name) for k, v in monosaccharide.modifications.items())
        null_positions = (UnknownPosition, NoPosition)
        modifications = "|" + modifications if modifications != "" else ""
        ring_start = monosaccharide.ring_start if monosaccharide.ring_start not in null_positions else 'x'
        ring_end = monosaccharide.ring_end if monosaccharide.ring_end not in null_positions else 'x'

        # The complete monosaccharide residue line
        residue_str = residue_template.format(ix=monosaccharide_index, anomer=anomer, conf_stem=conf_stem,
                                              superclass=superclass, modifications=modifications,
                                              ring_start=ring_start, ring_end=ring_end)
        return residue_str, monosaccharide_index

    def handle_monosaccharide(self, monosaccharide):  # pylint: disable=redefined-outer-name
        residue_str, monosaccharide_index = self._format_monosaccharide(monosaccharide)
        res = [residue_str]
        lin = []
        visited_subst = dict()
        # Construct the substituent lines
        # and their links
        for _lin_pos, link_obj in monosaccharide.substituent_links.items():
            sub = link_obj.to(monosaccharide)
            if sub.id not in visited_subst:
                sub_index = self.res_counter()
                subst_str = str(sub_index) + self.handle_substituent(sub)
                res.append(subst_str)
                visited_subst[sub.id] = sub_index
            lin.append(
                self.handle_link(
                    link_obj, self.lin_counter(), monosaccharide_index, visited_subst[sub.id]))

        return [res, lin, monosaccharide_index]

    def handle_glycan(self, structure):  # pragma: no cover
        if structure is None:
            raise GlycoCTError("No structure is ready to be written.")

        self.lin_accumulator = []
        self.buffer.write("RES\n")

        visited = set()
        for node in (structure):
            if node.id in visited:
                continue
            visited.add(node.id)
            res, lin, index = self.handle_monosaccharide(node)

            self.lin_accumulator.append((index, lin))
            self.residue_to_index[node.id] = index
            self.index_to_residue[index] = node

            if self.full:
                for _pos, lin in node.links.items():
                    if lin.is_child(node):
                        continue
                    self.dependencies[lin.child.id][node.id] = ((self.lin_counter(), lin))
            for line in res:
                self.buffer.write(line + '\n')

            # If this serialization is not meant to be full
            # do not visit residues beyond the first.
            if not self.full:
                break

        self.buffer.write("LIN\n")
        for res_ix, links in self.lin_accumulator:
            for line in links:
                self.buffer.write(line + '\n')
            residue = self.index_to_residue[res_ix]
            if self.full:
                for _pos, lin in residue.links.items():
                    if lin.is_child(residue):
                        continue
                    child_res = lin.child
                    ix, lin = self.dependencies[child_res.id][residue.id]
                    self.buffer.write(
                        self.handle_link(lin, ix, res_ix, self.residue_to_index[child_res.id]) + "\n")
        return self.buffer

    def begin_underdetermined(self):
        if self.state != UND:
            self.buffer.write("UND\n")
            self.state = UND
            self.index_to_residue.state = UND
            self.residue_to_index.state = UND

    def _format_subtree_linkage(self, linkage_args):
        (parent_link_type, parent_position,
         child_position, child_link_type) = linkage_args
        return "%s(%d+%d)%s" % (parent_link_type, parent_position,
                                child_position, child_link_type)

    def _get_viable_und_parents(self):
        valid_parent_inds = []
        for index, node in self.index_to_residue[START].items():
            if node.node_type == Monosaccharide.node_type:
                valid_parent_inds.append(index)
        return valid_parent_inds

    def handle_und_header(self, major_probability=100.0, minor_probability=100.0, parent_ids=None,
                          subtree_linkage_args=None):
        if parent_ids is None:
            parent_ids = list(self._get_viable_und_parents())

        index = self.und_counter()
        self.buffer.write("UND%d:%0.1f:%0.1f\n" % (index, major_probability, minor_probability))
        self.buffer.write(
            "ParentIDs:%s\n" % ('|'.join(map(str, parent_ids))))
        self.buffer.write("SubtreeLinkageID1:%s\n" % (
            self._format_subtree_linkage(subtree_linkage_args)))

    def dump(self):
        buffer = self.handle_glycan(self.structure)
        if self.nobuffer:
            value = buffer.getvalue()
            self._reset()
            return value
        return buffer

    def write(self, structure):
        self.structure = structure
        self._reset()
        return self.dump()

    def embed(self, writer):
        writer.res_counter = self.res_counter
        writer.lin_counter = self.lin_counter
        writer.und_counter = self.und_counter
        writer.index_to_residue = self.index_to_residue
        writer.buffer = self.buffer
        writer.state = self.state
        return writer

    def _determine_und_linkage_type_glycan_composition(self, glycan_composition):
        if len(glycan_composition) == 1:
            keys = list(glycan_composition)
            key = keys[0]
            if key.node_type == Substituent.node_type:
                return "d", "n"
            else:
                return "o", "d"
        else:
            return "o", "d"

    def add_glycan_composition(self, glycan_composition):
        for m in OrderingComparisonContext(self).sort_residues(glycan_composition):
            for _i in range(glycan_composition[m]):
                self.add_glycan_composition_single({m: 1})

    def add_glycan_composition_single(self, glycan_composition):
        self.begin_underdetermined()
        linkage_types = self._determine_und_linkage_type_glycan_composition(glycan_composition)
        self.handle_und_header(subtree_linkage_args=(linkage_types[0], -1, -1, linkage_types[1]))
        writer = GlycanCompositionGlycoCTWriter(
            glycan_composition, self.buffer)
        self.embed(writer)
        writer.handle_glycan(glycan_composition)


class GlycanCompositionGlycoCTWriter(GlycoCTWriterBase):
    def __init__(self, structure=None, buffer=None, full=True, standardize=False):
        super(GlycanCompositionGlycoCTWriter, self).__init__(structure, buffer, full)
        self.standardize = standardize

    @property
    def structure(self):
        return self._structure

    def _standardize_substituent_linkage(self, link):
        if self.standardize:
            from glypy.io.nomenclature import identity
            from glypy.structure.named_structures import monosaccharides
            try:
                if link.child.name == 'n_acetyl':
                    if identity.is_a(link.parent, monosaccharides.HexNAc):
                        if link.parent_position == -1:
                            link.parent_position = 2
            except AttributeError:
                pass

    @structure.setter
    def structure(self, value):
        from glypy.structure.glycan_composition import GlycanComposition

        if value is None:
            self._structure = value
            return
        if isinstance(value, (GlycanComposition, dict)):
            value = GlycanComposition(value)
        else:
            try:
                structure = treep(value)
            except TypeError:
                try:
                    root = rootp(value)
                    structure = Glycan(root, index_method=None)
                except TypeError:
                    raise TypeError("Could not extract or construct a tree structure from %r" % value)
            value = GlycanComposition.from_glycan(structure)
        self._structure = value

    def _unspool(self, mapping):
        sorter = OrderingComparisonContext(self)
        order = sorter.sort_residues(mapping.keys(), reverse=True)
        for key in order:
            count = mapping[key]
            if count < 1:
                continue
            for _i in range(count):
                yield key

    def _write_und_subgraph(self, substituent):  # pylint: disable=redefined-outer-name
        self.handle_und_header(subtree_linkage_args=('o', -1, 1, 'n'))
        self.buffer.write("RES\n")
        sub_index = self.res_counter()
        subst_str = str(sub_index) + self.handle_substituent(substituent)
        self.buffer.write("%s\n" % subst_str)

    def handle_link(self, link, ix, parent_ix, child_ix):
        self._standardize_substituent_linkage(link)
        return super(GlycanCompositionGlycoCTWriter, self).handle_link(link, ix, parent_ix, child_ix)

    def handle_glycan(self, structure):
        if structure is None:
            raise GlycoCTError("No structure is ready to be written.")

        from glypy.structure.glycan_composition import (
            SubstituentResidue, MolecularComposition)

        self.lin_accumulator = []

        disconnected_substituents = []
        molecules = []

        nodes = []

        for node in self._unspool(structure):
            if isinstance(node, SubstituentResidue):
                disconnected_substituents.append(node)
                continue
            elif isinstance(node, MolecularComposition):
                molecules.append(node)
                continue
            else:
                nodes.append(node)

        if nodes:
            self.buffer.write("RES\n")
            for node in nodes:
                res, lin, index = self.handle_monosaccharide(node)

                self.lin_accumulator.append((index, lin))
                self.residue_to_index[node.id] = index
                self.index_to_residue[index] = node

                if self.full:
                    for _pos, lin in node.links.items():
                        if lin.is_child(node):
                            continue
                        self.dependencies[lin.child.id][node.id] = ((self.lin_counter(), lin))
                for line in res:
                    self.buffer.write(line + '\n')

                # If this serialization is not meant to be full
                # do not visit residues beyond the first.
                if not self.full:
                    break
            # if self.lin_accumulator:
            if any(linkages for res_ix, linkages in self.lin_accumulator):
                self.buffer.write("LIN\n")
                for res_ix, links in self.lin_accumulator:
                    for line in links:
                        self.buffer.write(line + '\n')
                    residue = self.index_to_residue[res_ix]
                    if self.full:
                        for _pos, lin in residue.links.items():
                            if lin.is_child(residue):
                                continue
                            child_res = lin.child
                            ix, lin = self.dependencies[child_res.id][residue.id]
                            self.buffer.write(
                                self.handle_link(lin, ix, res_ix, self.residue_to_index[child_res.id]) + "\n")
            if disconnected_substituents:
                self.begin_underdetermined()
                for node in disconnected_substituents:
                    self._write_und_subgraph(node)
        elif disconnected_substituents:
            for subst in disconnected_substituents:
                self.buffer.write("RES\n")
                sub_index = self.res_counter()
                subst_str = str(sub_index) + self.handle_substituent(subst)
                self.buffer.write("%s\n" % subst_str)

        if molecules:
            raise TypeError("Cannot serialize MolecularComposition to GlycoCT")

        return self.buffer


def all_node_depth(node, visited=None):
    if visited is None:
        visited = set()
    if node.id in visited:  # pragma: no cover
        return 0
    visited.add(node.id)
    depth_count = 1
    children = list(node.children())
    try:
        children += list(node.substituents())
    except AttributeError:
        pass
    if children:
        depth_count += max(all_node_depth(ch, visited) for p, ch in children)
    return depth_count


class OrderingComparisonContext(object):
    def __init__(self, parent):
        self.parent = parent
        if isinstance(self.structure, Glycan) and not self.structure.has_index():
            self.structure.reindex()
        self.branch_to_terminal_count = self.build_branch_to_terminal_count()

    @property
    def structure(self):
        return self.parent.structure

    def get_branch_from_link_label(self, link):
        return link.label[0]

    def build_branch_to_terminal_count(self):
        counter = Counter()
        try:
            for key in sorted(self.structure.branch_parent_map.keys(), reverse=True):
                parent = self.structure.branch_parent_map[key]
                counter[parent] += counter[key] + 1
        except AttributeError:
            pass
        return counter

    def _residue_diff(self, res_a, res_b):
        n_child_residues_a = all_node_depth(res_a)
        n_child_residues_b = all_node_depth(res_b)
        diff_child_res = n_child_residues_a - n_child_residues_b

        if diff_child_res != 0:
            if diff_child_res < 0:
                diff_child_res = -1
            else:
                diff_child_res = 1

        try:
            branch_length_a = max((all_node_depth(cr) for p, cr in res_a.children()))
        except ValueError:
            branch_length_a = 0
        try:
            branch_length_b = max((all_node_depth(cr) for p, cr in res_b.children()))
        except ValueError:
            branch_length_b = 0

        diff_longest_branch = branch_length_a - branch_length_b

        if diff_longest_branch != 0:
            if diff_longest_branch < 0:
                diff_longest_branch = -1
            else:
                diff_longest_branch = 1

        n_branches_from_a = 0
        n_branches_from_b = 0
        for link in res_a.links.values():
            if link.is_parent(res_a):
                branch_label = self.get_branch_from_link_label(link)
                n_branches_from_a = max(n_branches_from_a, self.branch_to_terminal_count[branch_label])

        for link in res_b.links.values():
            if link.is_parent(res_b):
                branch_label = self.get_branch_from_link_label(link)
                n_branches_from_b = max(n_branches_from_b, self.branch_to_terminal_count[branch_label])
        diff_n_branches_from = n_branches_from_a - n_branches_from_b

        if diff_n_branches_from != 0:
            if diff_n_branches_from < 0:
                diff_n_branches_from = -1
            else:
                diff_n_branches_from = 1

        if res_a == res_b:
            subtree_diff = 0
        else:
            subtree_a = GlycoCTWriter(Glycan.subtree_from(self.structure, res_a)).dump()
            subtree_b = GlycoCTWriter(Glycan.subtree_from(self.structure, res_b)).dump()
            subtree_diff = (subtree_b > subtree_a) - (subtree_b < subtree_a)
            # cmp(subtree_b, subtree_a)
        return (diff_child_res, diff_longest_branch, diff_n_branches_from, subtree_diff, subtree_a, subtree_b)

    def _compare_residue_ordering(self, res_a, res_b):
        n_child_residues_a = all_node_depth(res_a)
        n_child_residues_b = all_node_depth(res_b)
        diff_child_res = n_child_residues_a - n_child_residues_b

        if diff_child_res != 0:
            if diff_child_res < 0:
                return -1
            else:
                return 1

        try:
            branch_length_a = max((all_node_depth(cr) for p, cr in res_a.children()))
        except ValueError:
            branch_length_a = 0
        try:
            branch_length_b = max((all_node_depth(cr) for p, cr in res_b.children()))
        except ValueError:
            branch_length_b = 0

        diff_longest_branch = branch_length_a - branch_length_b

        if diff_longest_branch != 0:
            if diff_longest_branch < 0:
                return -1
            else:
                return 1

        n_branches_from_a = 0
        n_branches_from_b = 0
        for link in res_a.links.values():
            if link.is_parent(res_a):
                branch_label = self.get_branch_from_link_label(link)
                n_branches_from_a = max(n_branches_from_a, self.branch_to_terminal_count[branch_label])

        for link in res_b.links.values():
            if link.is_parent(res_b):
                branch_label = self.get_branch_from_link_label(link)
                n_branches_from_b = max(n_branches_from_b, self.branch_to_terminal_count[branch_label])
        diff_n_branches_from = n_branches_from_a - n_branches_from_b

        if diff_n_branches_from != 0:
            if diff_n_branches_from < 0:
                return -1
            else:
                return 1

        if res_a == res_b:
            return 0

        subtree_a = GlycoCTWriter(Glycan.subtree_from(self.structure, res_a)).dump()
        subtree_b = GlycoCTWriter(Glycan.subtree_from(self.structure, res_b)).dump()
        return (subtree_b > subtree_a) - (subtree_b < subtree_a)

    def compare_residue_ordering(self, res_a, res_b):
        ordered = self._compare_residue_ordering(res_a, res_b)
        return ordered

    def _link_diff(self, link_a, link_b):  # pragma: no cover
        parent_pos_a = link_a.parent_position
        parent_pos_b = link_b.parent_position
        try:
            diff_parent = parent_pos_a - parent_pos_b
        except TypeError as e:
            print(parent_pos_a, parent_pos_b, link_a, link_b)
            raise e

        if diff_parent != 0:
            if diff_parent < 0:
                diff_parent = -1
            else:
                diff_parent = 1

        child_pos_a = link_a.child_position
        child_pos_b = link_b.child_position
        diff_child = child_pos_a - child_pos_b

        if diff_child != 0:
            if diff_child < 0:
                diff_child = -1
            else:
                diff_child = 1

        sigils_a = link_a._glycoct_sigils()
        sigils_b = link_b._glycoct_sigils()

        diff_sig0 = 0
        if sigils_a[0] != sigils_b[0]:
            diff_sig0 = ord(sigils_a[0]) - ord(sigils_b[0])
            if diff_sig0 < 0:
                diff_sig0 = -1
            else:
                diff_sig0 = 1

        diff_sig1 = 0
        if sigils_a[1] != sigils_b[1]:
            diff_sig1 = ord(sigils_a[1]) - ord(sigils_b[1])
            if diff_sig1 < 0:
                diff_sig1 = -1
            else:
                diff_sig1 = 1

        child_a = link_a.child
        child_b = link_b.child
        ordered = self.compare_residue_ordering(child_a, child_b)
        return (diff_parent, diff_child, diff_sig0, diff_sig1, ordered)

    def _compare_link_ordering(self, link_a, link_b):
        # Ignoring # of links for now since it is difficult
        # to compute
        parent_pos_a = link_a.parent_position
        parent_pos_b = link_b.parent_position
        try:
            diff_parent = parent_pos_a - parent_pos_b
        except TypeError as e:
            print(parent_pos_a, parent_pos_b, link_a, link_b)
            raise e

        if diff_parent != 0:
            if diff_parent < 0:
                return -1
            else:
                return 1

        child_pos_a = link_a.child_position
        child_pos_b = link_b.child_position
        diff_child = child_pos_a - child_pos_b

        if diff_child != 0:
            if diff_child < 0:
                return -1
            else:
                return 1

        sigils_a = link_a._glycoct_sigils()
        sigils_b = link_b._glycoct_sigils()

        if sigils_a[0] != sigils_b[0]:
            diff_sig0 = ord(sigils_a[0]) - ord(sigils_b[0])
            if diff_sig0 < 0:
                return -1
            else:
                return 1

        if sigils_a[1] != sigils_b[1]:
            diff_sig1 = ord(sigils_a[1]) - ord(sigils_b[1])
            if diff_sig1 < 0:
                return -1
            else:
                return 1

        child_a = link_a.child
        child_b = link_b.child
        ordered = self.compare_residue_ordering(child_a, child_b)
        return ordered

    def compare_link_ordering(self, link_a, link_b):
        ordered = self._compare_link_ordering(link_a, link_b)
        return ordered

    def sort_links(self, links, reverse=False):
        return sorted(links, key=cmp_to_key(self.compare_link_ordering),
                      reverse=reverse)

    def sort_residues(self, residues, reverse=False):
        return sorted(residues, key=cmp_to_key(self.compare_residue_ordering),
                      reverse=reverse)


class SubtreeJourney(object):
    def __init__(self, links):
        self.links = links

    def __repr__(self):
        template = "{self.__class__.__name__}({self.links!r})"
        return template.format(self=self)

    def __iter__(self):
        return iter(self.links)

    def __len__(self):
        return len(self.links)

    def __getitem__(self, i):
        return self.links[i]


class OrderedSubtreeTraverser(object):
    def __init__(self, ordering_context):
        self.ordering_context = ordering_context

    def visit_link(self, link, ignore_ambiguous=False):
        if isinstance(link, AmbiguousLink) and len(link.parent_choices) and not ignore_ambiguous:
            subtree = OrderedSubtreeTraverser(self.ordering_context)
            return [SubtreeJourney(subtree.visit_subtree(link, True))]
        node = link.child
        if node.node_type is Monosaccharide.node_type:
            next_links = self.visit_monosaccharide(node)
        else:
            next_links = self.visit_substituent(node)
        return next_links

    def visit_monosaccharide(self, monosaccharide):
        link_collection = list(monosaccharide.substituent_links.values())
        link_collection.extend(
            [cl for p, cl in monosaccharide.children(links=True)])
        links = self.ordering_context.sort_links(link_collection)
        return links[::-1]

    def visit_substituent(self, substituent):
        links = self.ordering_context.sort_links(
            [cl for p, cl in substituent.children(links=True)])
        return links[::-1]

    def visit_subtree(self, link, ignore_ambiguous=False):
        visited = set()
        link_queue = deque([link])
        journey = deque()

        while link_queue:
            link = link_queue.popleft()
            # Explicitly add before skipping to avoid double-writing
            # residues, but including multiple-link cases
            journey.append(link)
            if isinstance(link, SubtreeJourney):
                continue
            if link.child.id in visited:
                continue
            visited.add(link.child.id)
            link_queue.extend(
                self.visit_link(
                    link, ignore_ambiguous=ignore_ambiguous))
            ignore_ambiguous &= False
        return journey


class OrderRespectingGlycoCTWriter(GlycoCTWriterBase):
    def __init__(self, structure, buffer=None, full=True):
        super(OrderRespectingGlycoCTWriter, self).__init__(structure, buffer, full)

        self.ordering_context = OrderingComparisonContext(self)
        self.link_queue = deque()

    def handle_monosaccharide(self, monosaccharide):
        residue_str, monosaccharide_index = self._format_monosaccharide(monosaccharide)

        self.index_to_residue[monosaccharide_index] = monosaccharide
        self.residue_to_index[monosaccharide.id] = monosaccharide_index

        link_collection = list(monosaccharide.substituent_links.values())
        if self.full:
            link_collection.extend([cl for p, cl in monosaccharide.children(links=True)])

        links = self.ordering_context.sort_links(link_collection)
        self.link_queue.extendleft(links[::-1])
        return residue_str

    def handle_substituent(self, substituent):
        substituent_index = self.res_counter()

        self.index_to_residue[substituent_index] = substituent
        self.residue_to_index[substituent.id] = substituent_index

        subst_str = "%ss:%s" % (substituent_index, substituent.name.replace("_", "-"))

        links = self.ordering_context.sort_links([cl for p, cl in substituent.children(links=True)])
        self.link_queue.extendleft(links[::-1])
        return subst_str

    def handle_glycan(self, structure):
        if structure is None:
            raise GlycoCTError("No structure is ready to be written.")

        self.process_graph(structure.root)
        return self.buffer

    def process_graph(self, root, visited=None, link_queue=None):
        if visited is None:
            visited = set()

        if link_queue is None:
            link_queue = self.link_queue

        links_in_order = []
        underdetermined_subtrees = []

        self.buffer.write("RES\n")
        if root.node_type is Monosaccharide.node_type:
            res_str = self.handle_monosaccharide(root)
            self.buffer.write(res_str + "\n")
        else:
            res_str = self.handle_substituent(root)
            self.buffer.write(res_str + "\n")

        while link_queue:
            link = link_queue.popleft()
            # Explicitly add before skipping to avoid double-writing
            # residues, but including multiple-link cases
            links_in_order.append(link)
            if link.child.id in visited:
                continue
            visited.add(link.child.id)
            if link.child.node_type is Monosaccharide.node_type:
                line = self.handle_monosaccharide(link.child)
            else:
                line = self.handle_substituent(link.child)
            self.buffer.write(line + "\n")

        self.buffer.write("LIN\n")
        for link in links_in_order:
            if not link.is_substituent_link() and not self.full:
                continue
            parent_ix = self.residue_to_index[link.parent.id]
            child_ix = self.residue_to_index[link.child.id]
            line = self.handle_link(
                link, self.lin_counter(), parent_ix, child_ix)
            self.buffer.write(line + "\n")

        return visited, underdetermined_subtrees


class UNDOrderRespectingGlycoCTWriter(OrderRespectingGlycoCTWriter):
    def gobble_subtree(self, link):
        traverser = OrderedSubtreeTraverser(self.ordering_context)
        return SubtreeJourney(traverser.visit_subtree(link, True))

    def process_graph(self, root, visited=None, link_queue=None):
        if visited is None:
            visited = set()

        if link_queue is None:
            link_queue = self.link_queue

        links_in_order = []
        underdetermined_subtrees = []

        self.buffer.write("RES\n")
        if root.node_type is Monosaccharide.node_type:
            res_str = self.handle_monosaccharide(root)
            self.buffer.write(res_str + "\n")
        else:
            res_str = self.handle_substituent(root)
            self.buffer.write(res_str + "\n")

        while link_queue:
            link = link_queue.popleft()
            if isinstance(link, AmbiguousLink) and len(link.parent_choices) > 1:
                underdetermined_subtrees.append(self.gobble_subtree(link))
                continue
            # Explicitly add before skipping to avoid double-writing
            # residues, but including multiple-link cases
            links_in_order.append(link)
            if link.child.id in visited:
                continue
            visited.add(link.child.id)
            if link.child.node_type is Monosaccharide.node_type:
                line = self.handle_monosaccharide(link.child)
            else:
                line = self.handle_substituent(link.child)
            self.buffer.write(line + "\n")

        self.buffer.write("LIN\n")
        for link in links_in_order:
            if not link.is_substituent_link() and not self.full:
                continue
            parent_ix = self.residue_to_index[link.parent.id]
            child_ix = self.residue_to_index[link.child.id]
            line = self.handle_link(
                link, self.lin_counter(), parent_ix, child_ix)
            self.buffer.write(line + "\n")

        return visited, underdetermined_subtrees

    def handle_glycan(self, structure):
        if structure is None:
            raise GlycoCTError("No structure is ready to be written.")

        visited, underdetermined_subtrees = self.process_graph(structure.root)

        if underdetermined_subtrees:
            self.begin_underdetermined()
            for i, und in enumerate(underdetermined_subtrees, 1):
                und = deque(und)
                self.buffer.write("UND%d:100.0:100.0\n" % i)
                base_link = und.popleft()
                parent_ids = "|".join([str(self.residue_to_index[p.id]) for p in base_link.parent_choices])
                self.buffer.write("ParentIDs:%s\n" % parent_ids)
                subtree_linkage = self.handle_link(base_link, '-', '', '')[2:]
                self.buffer.write("SubtreeLinkageID1:%s\n" % subtree_linkage)

                outer_link_deque = self.link_queue
                self.link_queue = deque([])

                child = base_link.child
                self.process_graph(child, visited=visited, link_queue=und)

                self.link_queue = outer_link_deque

        return self.buffer


GlycoCTWriter = UNDOrderRespectingGlycoCTWriter


def dump(structure, buffer=None):
    '''
    Serialize the |Glycan| into :title-reference:`GlycoCT{condensed}`, using
    `buffer` to store the result. If `buffer` is |None|, then the
    function will operate on a newly created :class:`StringIO` object.

    Parameters
    ----------
    structure: |Glycan|
        The structure to serialize
    buffer: file-like or None
        The stream to write the serialized structure to. If |None|, uses an instance
        of :class:`StringIO`

    Returns
    -------
    file-like or str if ``buffer`` is :const:`None`
    '''
    from glypy import GlycanComposition
    if isinstance(structure, GlycanComposition):
        return GlycanCompositionGlycoCTWriter(structure, buffer).dump()
    return GlycoCTWriter(structure, buffer).dump()


def dumps(structure):
    '''
    Serialize the |Glycan| into :title-reference:`GlycoCT{condensed}`, returning
    the text as a string.

    Parameters
    ----------
    structure: |Glycan|
        The structure to serialize

    Returns
    -------
    str
    '''
    from glypy import GlycanComposition
    if isinstance(structure, GlycanComposition):
        return GlycanCompositionGlycoCTWriter(structure, None).dump()
    return GlycoCTWriter(structure, None).dump()


def _postprocessed_single_monosaccharide(monosaccharide, convert=True):
    if convert:
        monostring = GlycoCTWriterBase(monosaccharide, None, full=False).dump()
    else:
        monostring = monosaccharide
    monostring = monostring.replace("\n", " ")
    if monostring.endswith("LIN "):
        monostring = monostring.replace(" LIN ", "")
    else:
        monostring = monostring.strip()
    return monostring


Monosaccharide.register_serializer("glycoct", _postprocessed_single_monosaccharide)
Glycan.register_serializer("glycoct", dumps)


def canonicalize(structure):
    return loads(dumps(structure))
