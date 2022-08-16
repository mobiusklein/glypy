from collections import defaultdict
from weakref import ref

from six import string_types as basestring

from glypy.utils import uid
from glypy.composition import Composition

from . import atomic_data


class Graph(object):
    def __init__(self, vertices=None, edges=None):
        self.vertices = vertices or {}
        self.edges = edges or set()
        self.adjacency_list = defaultdict(set)

    def __repr__(self):
        template = "{self.__class__.__name__}({self.vertices}, {self.edges})"
        return template.format(self=self)

    def add_vertex(self, vertex):
        self.vertices[vertex.id] = vertex
        vertex.graph = ref(self)
        return self

    def __getitem__(self, i):
        return self.vertices[i]

    def __iter__(self):
        keys = sorted(self.vertices.keys())
        for key in keys:
            yield self.vertices[key]

    def __len__(self):
        return len(self.vertices)

    def add_edge(self, edge):
        self.edges.add(edge)
        self.adjacency_list[edge.vertex_ids[0]].add((edge.vertex_ids[1], edge))
        self.adjacency_list[edge.vertex_ids[1]].add((edge.vertex_ids[0], edge))
        return self

    def neighbors_of(self, vertex_id):
        neighbors = [self[x[0]] for x in self.adjacency_list[vertex_id]]
        return neighbors

    def edges_of(self, vertex_id):
        edges = [x[1] for x in self.adjacency_list[vertex_id]]
        return edges


class Vertex(object):
    def __init__(self, id=None, _graph=None, **node_data):
        if id is None:
            id = uid()
        self.id = id
        self.data = node_data
        self._graph = None
        self.graph = _graph

    @property
    def graph(self):
        if self._graph is not None:
            return self._graph()
        return None

    @graph.setter
    def graph(self, value):
        self._graph = value

    def __repr__(self):
        if self.data:
            data = ', {self.data}'.format(self=self)
        else:
            data = ''
        template = "{self.__class__.__name__}({self.id}{data})"
        return template.format(self=self, data=data)

    def neighbors(self):
        graph = self.graph
        if graph is None:
            return []
        neighbors = graph.neighbors_of(self.id)
        return neighbors


class Edge(object):
    def __init__(self, vertex_ids):
        self.vertex_ids = tuple(vertex_ids)

    def __repr__(self):
        template = '{self.__class__.__name__}({self.vertex_ids})'
        return template.format(self=self)

    def __eq__(self, other):
        return self.vertex_ids == other.vertex_ids

    def __hash__(self):
        return hash(self.vertex_ids)


class Atom(Vertex):
    def __init__(self, element, id=None, _graph=None, **node_data):
        if isinstance(element, basestring):
            element = atomic_data.element_information[element]
        super(Atom, self).__init__(id, _graph=_graph, **node_data)
        self.element = element

    def __repr__(self):
        if self.data:
            data = ', {self.data}'.format(self=self)
        else:
            data = ''
        template = "{self.__class__.__name__}({self.element.symbol!r}, {self.id}{data})"
        return template.format(self=self, data=data)

    @property
    def valence(self):
        return self.element.valence

    @property
    def unpaired_electrons(self):
        theoretical = self.valence
        graph = self.graph
        bonds = graph.edges_of(self.id)
        for bond in bonds:
            theoretical -= bond.order
        return theoretical


class Bond(Edge):
    def __init__(self, vertex_ids, order=1):
        super(Bond, self).__init__(vertex_ids)
        self.order = order

    def __repr__(self):
        template = '{self.__class__.__name__}({self.vertex_ids}, {self.order})'
        return template.format(self=self)


class MolecularGraph(Graph):

    def total_composition(self):
        composition = Composition()
        unpaired_electrons = 0
        for node in self:
            composition[node.element.symbol] += 1
            unpaired_electrons += node.unpaired_electrons
        composition['H'] = unpaired_electrons
        return composition

    def mass(self):
        return self.total_composition().mass
