'''
A GlycoCT (Condensed) parser.
Supports RES, LIN, and un-nested REP sections.

.. code-block:: python

    >>>from glypy.io import glycoct
    >>>glycoct.loads("""RES
    1b:o-dman-HEX-0:0|1:aldi
    2b:a-lido-HEX-1:5|6:a
    3s:sulfate
    LIN
    1:1o(3+1)2d
    2:2o(2+1)3n
    3:2o(4+1)4d
    """)
    >>>
'''

import re
import warnings
from collections import defaultdict, Counter, deque, namedtuple, OrderedDict
from functools import cmp_to_key
import zlib

from glypy.utils import (
    opener, StringIO, root as rootp, tree as treep,
    make_counter, invert_dict, uid,
    RootProtocolNotSupportedError)
from glypy.utils.multimap import OrderedMultiMap
from glypy.structure import monosaccharide, substituent, glycan, Modification
from glypy.structure.link import Link, AmbiguousLink
from .format_constants_map import (anomer_map, superclass_map,
                                   link_replacement_composition_map, modification_map)
from .file_utils import ParserError
from .tree_builder_utils import (
    decorate_tree,
    undecorate_tree,
    find_root,
    try_int,
    StructurePrecisionEnum,
    AbstractGraphEntryEnum)
from glypy.composition import Composition

try:
    range = xrange
except NameError:
    pass


__id = id

Glycan = glycan.Glycan
Monosaccharide = monosaccharide.Monosaccharide
Substituent = substituent.Substituent

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
    '''
    (?P<anomer>[abxo])?
    (?P<conf_stem>(?:-[dlx][a-z]+)+)?-?
    (?P<superclass>[A-Z]+)-?
    (?P<indices>[0-9x]+:[0-9x]+)
    (?P<modifications>(\|[0-9x]+:[0-9a-z]+)+)?
    ''', re.VERBOSE)

#: Pattern for parsing the potentially repeated |Configuration| and |Stem|
#: regions of the lines of the RES section.
conf_stem_pattern = re.compile(r'(?P<config>[dlx])(?P<stem>[a-z]+)')

#: Pattern for parsing modifications found on monosaccharide residue
#: lines in the RES section
modification_pattern = re.compile(r"\|?(\d+):([^\|;\n]+)")


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

repeat_line_pattern = re.compile("^(?P<graph_index>\d+)r:r(?P<repeat_index>\d+)")

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
    pass


class GlycoCTSectionUnsupported(GlycoCTError):
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
        return self.graph.keys()

    def values(self):
        return self.graph.values()

    def items(self):
        return self.graph.items()

    def clear(self):
        self.graph.clear()

    def __iter__(self):
        return iter(self.graph)

    def get_node(self, id, direction=None):
        id = int(id)
        return self.graph[id]

    def put_node(self, id, value):
        id = int(id)
        self.graph[id] = value

    def form_link(self, parent, child, parent_position, child_position, parent_loss,
                  child_loss, id=None):
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
                child_loss=child_loss, id=id)
            link_obj.find_open_position()
        else:
            link_obj = Link(
                parent, child, parent_position=int(parent_position[0]),
                child_position=int(child_position[0]), parent_loss=parent_loss,
                child_loss=child_loss)
        return link_obj

    def deferred_retrieval(self, id, direction=None):
        return DeferredRetrieval(self, id, direction)

    def __root__(self):
        return self.find_root_nodes()[0]

    def find_root_nodes(self):
        roots = []
        for index, node in self.items():
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

    def visit(self, node, visited=None):
        if visited is None:
            visited = set()
        if isinstance(node, GlycoCTGraph):
            return self.visit(rootp(node), visited)

        visited.add(node.id)
        for position, link in node.links.items():
            ref = link.to(node)
            if ref.id in visited:
                continue
            else:
                self.visit(ref, visited)
        try:
            for position, link in node.substituent_links.items():
                ref = link.to(node)
                if ref.id in visited:
                    continue
                else:
                    self.visit(ref, visited)
        except AttributeError:
            pass
        return visited

    def is_fully_connected(self):
        roots = self.find_root_nodes()
        visited = self.visit(roots[0])
        return len(visited) >= len(self)


class NodeCollection(object):

    @classmethod
    def from_node(cls, root):
        nodes = []
        for node in monosaccharide.traverse(root):
            nodes.append(node)
            try:
                for _, subst in node.substituents():
                    nodes.append(subst)
            except AttributeError:
                pass
        return cls(nodes)

    def __init__(self, nodes):
        self.nodes = list(nodes)
        self.root = self.nodes[0]

    def __iter__(self):
        return monosaccharide.traverse(self.root)

    def deindex(self):
        base = uid()
        for i, node in enumerate(self.nodes):
            node.id += base
            node.id *= -1

    def reindex(self):
        for i, node in enumerate(self.nodes):
            node.id = i + 1

    def clone(self, *args, **kwargs):
        root = monosaccharide.graph_clone(self.root)
        return self.from_node(root)

    def get(self, key):
        for node in self.nodes:
            if node.id == key:
                return node
        raise IndexError(key)

    def __root__(self):
        return self.root


class GlycoCTGraphStack(GlycoCTGraph):
    def __init__(self, stack=None, parent=None):
        if stack is None:
            stack = deque([GlycoCTSubgraph(parent=parent)])
        else:
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
        return self.stack[-1]

    @property
    def parent(self):
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

    def add(self, subgraph, parent=None):
        if subgraph.parent is None:
            if parent is None:
                parent = self
            subgraph.parent = parent
        self.push_level(subgraph)

    def push_level(self, subgraph):
        self.stack.append(subgraph)
        self.history.append(subgraph)

    def pop_level(self):
        return self.stack.pop()

    def clear(self):
        self.stack = deque([GlycoCTSubgraph(parent=None)])
        self.history = list(self.stack)

    def find_subgraph_containing(self, id):
        for level in reversed(self.history):
            if id in level:
                return level
        else:
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
    def __init__(self, graph=None, parent=None):
        super(GlycoCTSubgraph, self).__init__(graph)
        assert not isinstance(parent, dict)
        self.parent = parent

    def postprocess(self):
        pass


RepeatedMultitude = namedtuple("RepeatedMultitude", "lower upper")


class RepeatedGlycoCTSubgraph(GlycoCTSubgraph):
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
        return int(self.external_linkage['parent_residue_index'])

    @property
    def last_repeat_index(self):
        sub_unit_indices = sorted(self.repetitions.keys())
        terminal_unit_ix = sub_unit_indices[-1]
        return terminal_unit_ix

    @property
    def terminal_node(self):
        """Retrieves the terminal residue from the subgraph, where
        outgoing connections start from

        Connections between repeats of this subgraph will start from
        this node, as will connections between the final repeat of this
        subgraph and its parent graph

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
        sub_unit_indices = sorted(self.repetitions.keys())
        return sub_unit_indices[0]

    @property
    def origin_node(self):
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
        self.postponed.append((f, args))

    def complete_postponed_tasks(self):
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
        except RootProtocolNotSupportedError:
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

    def postprocess(self, n=None):
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
        parent_attachment_position = self.internal_linkage["parent_attachment_position"]

        child_residue_index = self.origin_node_index
        child_atom_replaced = link_replacement_composition_map[self.internal_linkage["child_atom_replaced"]]
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
                      parent_loss=parent_atom_replaced, child_loss=child_atom_replaced)))

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

    def prepare_glycan(self):
        glycan = self.repetitions[self.first_repeat_index]
        glycan.deindex()
        return Glycan(glycan.root, index_method=None)


UndeterminedProbability = namedtuple("UndeterminedProbability", "major minor")


class UnderdeterminedRecord(GlycoCTSubgraph):
    def __init__(self, und_index, probability=None, parent_ids=None,
                 subtree_linkages=None, graph=None, parent=None):
        super(UnderdeterminedRecord, self).__init__(graph, parent)
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
        except RootProtocolNotSupportedError:
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
        except ValueError:
            if link_obj.child_position == 1 and Modification.Acidic in link_obj.child.modifications[1]:
                link_obj.child_position = 2
                ix = link_obj.child_position_choices.index(1)
                link_obj.child_position_choices.pop(ix)
                link_obj.child_position_choices.insert(ix, 2)
                link_obj.apply()
                link_obj.find_open_position()
            else:
                raise


def _build_graph(glycoct_str):
    rep = StringIO(glycoct_str)
    inst = GlycoCTReader(rep, completes=False)
    return next(inst)


def extract_composition(parser):
    from glypy.structure.glycan_composition import (
        GlycanComposition, MonosaccharideResidue, SubstituentResidue)
    store = GlycanComposition()
    for layer in list(parser.stack)[1:]:
        node = rootp(layer)
        for position, link in list(node.links.items()):
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


class GlycoCTReader(GlycoCTGraphStack):
    '''
    Simple State-Machine parser for condensed GlycoCT representations. Yields
    |Glycan| instances.
    '''

    @classmethod
    def loads(cls, glycoct_str, structure_class=Glycan, allow_repeats=True):
        '''Parse results from |str|'''
        rep = StringIO(glycoct_str)
        return cls(rep, structure_class=structure_class, allow_repeats=allow_repeats)

    def __init__(self, stream, structure_class=Glycan, allow_repeats=True, completes=True):
        '''
        Creates a parser of condensed GlycoCT.

        Parameters
        ----------
        stream: basestring or file-like
            A path to a file or a file-like object to be processed
        '''
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
        # was_last_line_blank = False
        for line in self.handle:
            self._source_line += 1
            # if line.strip() == "":
            #     if not was_last_line_blank:
            #         self._current_segment = line
            #         yield line
            #         was_last_line_blank = True
            # else:
            #     was_last_line_blank = False
            #     for segment in re.split(r"\s|;", line):
            #         if "" == segment.strip():
            #             continue
            #         self._current_segment = segment
            #         yield segment
            for segment in re.split(r"\s|;", line):
                if "" == segment.strip():
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

    def handle_residue_line(self, line):
        '''
        Handle a base line, creates an instance of |Monosaccharide|
        and adds it to :attr:`graph` at the given index.

        Called by :meth:`parse`
        '''
        _, ix, residue_str = re.split(r"^(\d+)b", line, maxsplit=1)
        residue_dict = res_pattern.search(residue_str).groupdict()

        mods = residue_dict.pop("modifications")
        modifications = OrderedMultiMap()
        if mods is not None:
            for p, mod in modification_pattern.findall(mods):
                modifications[try_int(p)] = modification_map[mod]

        residue_dict["modifications"] = modifications
        is_reduced = "aldi" in modifications[1]
        if is_reduced:
            modifications.pop(1, "aldi")
            residue_dict['reduced'] = True

        conf_stem = residue_dict.pop("conf_stem")
        if conf_stem is not None:
            config, stem = zip(*conf_stem_pattern.findall(conf_stem))
        else:
            config = ('x',)
            stem = ('x',)
        residue_dict['stem'] = stem
        residue_dict['configuration'] = config

        residue_dict["ring_start"], residue_dict["ring_end"] = list(map(
            try_int, residue_dict.pop("indices").split(":")))

        residue_dict['anomer'] = anomer_map[residue_dict['anomer']]
        residue_dict['superclass'] = superclass_map[residue_dict['superclass']]
        residue = monosaccharide.Monosaccharide(**residue_dict)

        self[int(ix)] = residue

        residue.id = int(ix)
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

        child_atom_replaced = link_replacement_composition_map[link_dict["child_atom_replaced"]]
        child_attachment_position = list(map(int, link_dict["child_attachment_position"].split("|")))

        return (id, parent_residue_index, parent_atom_replaced, parent_attachment_position,
                child_residue_index, child_atom_replaced, child_attachment_position)

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
            child_residue_index, child_atom_replaced, child_attachment_position = self.parse_link(line)

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
                parent_loss=parent_atom_replaced, child_loss=child_atom_replaced, id=id)

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
        record = UnderdeterminedRecord(
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
        # if True:
            try:
                inst = undecorate_tree(
                    self.structure_class(
                        root=rootp(self.root), index_method=None)
                ).reindex()
                return inst
            except RootProtocolNotSupportedError:
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


def load(stream, structure_class=Glycan, allow_repeats=True, allow_multiple=True):
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


def loads(glycoct_str, structure_class=Glycan, allow_repeats=True, allow_multiple=True):
    '''
    A convenience wrapper for :meth:`GlycoCTReader.loads`

    As additional convenience, this function does not return an
    iterator over glycans, and returns a single instance if only
    one is present, or a list of instances otherwise.
    '''

    g = GlycoCTReader.loads(glycoct_str, structure_class=structure_class, allow_repeats=allow_repeats)
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


def detect_glycoct(string):
    return string.lstrip()[:3] == "RES"


invert_anomer_map = invert_dict(anomer_map)
invert_superclass_map = invert_dict(superclass_map)


class GlycoCTWriterBase(object):
    def __init__(self, structure=None, buffer=None, full=True):
        self.nobuffer = False
        if buffer is None:
            buffer = StringIO()
            self.nobuffer = True

        self.buffer = buffer
        self.structure = structure
        self.full = full

        self.res_counter = make_counter()
        self.lin_counter = make_counter()
        self.und_counter = make_counter()

        # Look-ups for mapping RES nodes to objects by section index and id,
        # respectively
        self.index_to_residue = {}
        self.residue_to_index = {}

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
        self.res_counter = make_counter()
        self.lin_counter = make_counter()
        self.und_counter = make_counter()

        self.index_to_residue = {}
        self.residue_to_index = {}

        self.lin_accumulator = []
        self.dependencies = defaultdict(dict)

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

    def handle_substituent(self, substituent):
        return "s:{0}".format(substituent.name.replace("_", "-"))

    def handle_monosaccharide(self, monosaccharide):
        residue_template = "{ix}b:{anomer}{conf_stem}{superclass}-{ring_start}:{ring_end}{modifications}"

        # This index is reused many times
        monosaccharide_index = self.res_counter()

        # Format individual fields
        anomer = invert_anomer_map[monosaccharide.anomer]
        conf_stem = ''.join("-{0}{1}".format(c.name, s.name)
                            for c, s in zip(monosaccharide.configuration, monosaccharide.stem))
        if None in monosaccharide.configuration or None in monosaccharide.stem:
            conf_stem = ''
        superclass = "-" + invert_superclass_map[monosaccharide.superclass]

        modifications = '|'.join(
            "{0}:{1}".format(k, v.name) for k, v in monosaccharide.modifications.items())

        modifications = "|" + modifications if modifications != "" else ""
        ring_start = monosaccharide.ring_start if monosaccharide.ring_start is not None else 'x'
        ring_end = monosaccharide.ring_end if monosaccharide.ring_end is not None else 'x'

        # The complete monosaccharide residue line
        residue_str = residue_template.format(ix=monosaccharide_index, anomer=anomer, conf_stem=conf_stem,
                                              superclass=superclass, modifications=modifications,
                                              ring_start=ring_start, ring_end=ring_end)
        res = [residue_str]
        lin = []
        visited_subst = dict()
        # Construct the substituent lines
        # and their links
        for lin_pos, link_obj in monosaccharide.substituent_links.items():
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

    def handle_glycan(self):
        if self.structure is None:
            raise GlycoCTError("No structure is ready to be written.")

        self.buffer.write("RES\n")

        visited = set()
        for node in (self.structure):
            if node.id in visited:
                continue
            visited.add(node.id)
            res, lin, index = self.handle_monosaccharide(node)

            self.lin_accumulator.append((index, lin))
            self.residue_to_index[node.id] = index
            self.index_to_residue[index] = node

            if self.full:
                for pos, lin in node.links.items():
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
                for pos, lin in residue.links.items():
                    if lin.is_child(residue):
                        continue
                    child_res = lin.child
                    ix, lin = self.dependencies[child_res.id][residue.id]
                    self.buffer.write(
                        self.handle_link(lin, ix, res_ix, self.residue_to_index[child_res.id]) + "\n")
        return self.buffer

    def dump(self):
        buffer = self.handle_glycan()
        if self.nobuffer:
            value = buffer.getvalue()
            self._reset()
            return value
        return buffer

    def write(self, structure):
        self.structure = structure
        self._reset()
        return self.dump()


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
            for i in range(count):
                yield key

    def _write_und_subgraph(self, substituent):
        index = self.und_counter()
        self.buffer.write("UND%d:100.0:100.0\n" % (index,))
        self.buffer.write(
            "ParentIDs:%s\n" % ('|'.join(map(str, self.index_to_residue.keys()))))
        self.buffer.write("SubtreeLinkageID1:o(-1+1)n\nRES\n")
        sub_index = self.res_counter()
        subst_str = str(sub_index) + self.handle_substituent(substituent)
        self.buffer.write("%s\n" % subst_str)

    def handle_link(self, link, ix, parent_ix, child_ix):
        self._standardize_substituent_linkage(link)
        return super(GlycanCompositionGlycoCTWriter, self).handle_link(link, ix, parent_ix, child_ix)

    def handle_glycan(self):
        if self.structure is None:
            raise GlycoCTError("No structure is ready to be written.")

        from glypy.structure.glycan_composition import (
            SubstituentResidue, MolecularComposition)

        self.buffer.write("RES\n")

        disconnected_substituents = []
        molecules = []

        for node in self._unspool(self.structure):
            if isinstance(node, SubstituentResidue):
                disconnected_substituents.append(node)
                continue
            elif isinstance(node, MolecularComposition):
                molecules.append(node)
                continue
            res, lin, index = self.handle_monosaccharide(node)

            self.lin_accumulator.append((index, lin))
            self.residue_to_index[node.id] = index
            self.index_to_residue[index] = node

            if self.full:
                for pos, lin in node.links.items():
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
                for pos, lin in residue.links.items():
                    if lin.is_child(residue):
                        continue
                    child_res = lin.child
                    ix, lin = self.dependencies[child_res.id][residue.id]
                    self.buffer.write(
                        self.handle_link(lin, ix, res_ix, self.residue_to_index[child_res.id]) + "\n")
        if disconnected_substituents:
            self.buffer.write("UND\n")
        for node in disconnected_substituents:
            self._write_und_subgraph(node)

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

    def _link_diff(self, link_a, link_b):
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


class OrderRespectingGlycoCTWriter(GlycoCTWriterBase):
    def __init__(self, structure, buffer=None, full=True):
        self.nobuffer = False
        if buffer is None:
            buffer = StringIO()
            self.nobuffer = True

        self.ordering_context = OrderingComparisonContext(self)

        self.buffer = buffer
        self.structure = structure
        self.full = full

        self.link_queue = deque()
        self.res_counter = make_counter()
        self.lin_counter = make_counter()

        # Look-ups for mapping RES nodes to objects by section index and id,
        # respectively
        self.index_to_residue = {}
        self.residue_to_index = {}

        # Accumulator for linkage indices and mapping linkage indices to
        # dependent RES indices
        self.lin_accumulator = []
        self.dependencies = defaultdict(dict)

    def handle_monosaccharide(self, monosaccharide):
        residue_template = "{ix}b:{anomer}{conf_stem}{superclass}-{ring_start}:{ring_end}{modifications}"

        # This index is reused many times
        monosaccharide_index = self.res_counter()

        self.index_to_residue[monosaccharide_index] = monosaccharide
        self.residue_to_index[monosaccharide.id] = monosaccharide_index

        # Format individual fields
        anomer = invert_anomer_map[monosaccharide.anomer]
        conf_stem = ''.join("-{0}{1}".format(c.name, s.name)
                            for c, s in zip(monosaccharide.configuration, monosaccharide.stem))
        superclass = "-" + invert_superclass_map[monosaccharide.superclass]

        modifications = '|'.join(
            "{0}:{1}".format(k, v.name) for k, v in monosaccharide.modifications.items())

        modifications = "|" + modifications if modifications != "" else ""
        ring_start = monosaccharide.ring_start if monosaccharide.ring_start is not None else 'x'
        ring_end = monosaccharide.ring_end if monosaccharide.ring_end is not None else 'x'

        # The complete monosaccharide residue line
        residue_str = residue_template.format(ix=monosaccharide_index, anomer=anomer, conf_stem=conf_stem,
                                              superclass=superclass, modifications=modifications,
                                              ring_start=ring_start, ring_end=ring_end)

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

    def handle_glycan(self):
        if self.structure is None:
            raise GlycoCTError("No structure is ready to be written.")

        self.buffer.write("RES\n")

        visited = set()
        if self.structure.root.node_type is Monosaccharide.node_type:
            res_str = self.handle_monosaccharide(self.structure.root)
            self.buffer.write(res_str + "\n")
        else:
            res_str = self.handle_substituent(self.structure.root)
            self.buffer.write(res_str + "\n")
        links_in_order = []
        while self.link_queue:
            link = self.link_queue.popleft()
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
            line = self.handle_link(link, self.lin_counter(), parent_ix, child_ix)
            self.buffer.write(line + "\n")

        return self.buffer


GlycoCTWriter = OrderRespectingGlycoCTWriter


def dump(structure, buffer=None):
    '''
    Serialize the |Glycan| graph object into condensed GlycoCT, using
    `buffer` to store the result. If `buffer` is |None|, then the
    function will operate on a newly created :class:`~glypy.utils.StringIO` object.

    Parameters
    ----------
    structure: |Glycan|
        The structure to serialize
    buffer: file-like or None
        The stream to write the serialized structure to. If |None|, uses an instance
        of `StringIO`

    Returns
    -------
    file-like or str if ``buffer`` is :const:`None`

    '''
    from glypy import GlycanComposition
    if isinstance(structure, GlycanComposition):
        return GlycanCompositionGlycoCTWriter(structure, buffer).dump()
    return GlycoCTWriter(structure, buffer).dump()


def dumps(structure, full=True):
    from glypy import GlycanComposition
    if isinstance(structure, GlycanComposition):
        return GlycanCompositionGlycoCTWriter(structure, None).dump()
    return GlycoCTWriter(structure, None, full=full).dump()


def _postprocessed_single_monosaccharide(monosaccharide, convert=True):
    if convert:
        monostring = dumps(monosaccharide, full=False)
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


class DistinctGlycanSet(object):

    def __init__(self, structures=None):
        if structures is None:
            structures = []
        self.raw_data_buffer = set()
        for structure in structures:
            self.add(structure)

    def add(self, structure):
        key = self._transform_text(
            self._structure_to_text(structure))
        if key in self.raw_data_buffer:
            return key
        self.raw_data_buffer.add(key)
        return key

    def _structure_to_text(self, structure):
        return dumps(structure)

    def _text_to_structure(self, text):
        return loads(text)

    def _transform_text(self, text):
        return zlib.compress(text)

    def _untransform_text(self, compressed):
        return zlib.decompress(compressed)

    def pop(self):
        text = self._untransform_text(self.raw_data_buffer.pop())
        return self._text_to_structure(text)

    def __len__(self):
        return len(self.raw_data_buffer)

    def __iter__(self):
        for compressed in self.raw_data_buffer:
            yield self._text_to_structure(
                self._untransform_text(compressed))
