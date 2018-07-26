from collections import OrderedDict

from .node_type import NodeTypeSpec
from .utils import base52


class WURCSWriter(object):
    version = '2.0'

    def __init__(self, glycan):
        self.glycan = glycan
        self.node_type_map = OrderedDict()
        self.node_index_to_node_type = OrderedDict()
        self.index_to_glyph = OrderedDict()
        self.id_to_index = OrderedDict()
        self.extract_node_types()

    def extract_node_types(self):
        node_types = OrderedDict()
        node_index_to_node_type = OrderedDict()
        index_to_glyph = OrderedDict()
        id_to_index = OrderedDict()
        for i, node in enumerate(self.glycan, 1):
            node_type = NodeTypeSpec.from_monosaccharide(node)
            index_to_glyph[i] = base52(i - 1)
            id_to_index[node.id] = i
            node_index_to_node_type[i] = node_type
            if node_type not in node_types:
                node_types[node_type] = len(node_types) + 1
        self.node_type_map = node_types
        self.node_index_to_node_type = node_index_to_node_type
        self.index_to_glyph = index_to_glyph
        self.id_to_index = id_to_index

    def format_version(self):
        return "WURCS=%s" % (self.version, )

    def format_count_section(self):
        count_nodes = len(list(self.glycan.iternodes()))
        count_links = len(list(self.glycan.iterlinks()))
        count_section = "%s,%s,%s" % (len(self.node_type_map), count_nodes, count_links)
        return count_section

    def format_node_types(self):
        return ''.join('[%s]' % (str(s),) for s in self.node_type_map.keys())

    def format_node_type_index(self):
        node_type_sequence = []
        for index, node_type in self.node_index_to_node_type.items():
            node_type_sequence.append(self.node_type_map[node_type])
        return '-'.join(map(str, node_type_sequence))

    def format_links(self):
        links = []
        for _, link in self.glycan.iterlinks():
            parent_index = self.id_to_index[link.parent.id]
            child_index = self.id_to_index[link.child.id]
            parent_glyph = self.index_to_glyph[parent_index]
            child_glyph = self.index_to_glyph[child_index]
            parent_position = link.parent_position
            if parent_position == -1:
                parent_position = '?'
            child_position = link.child_position
            if child_position == -1:
                child_position = '?'
            link_spec = '%s%s-%s%s' % (parent_glyph, parent_position, child_glyph, child_position)
            links.append(link_spec)
        return '_'.join(links)

    def write(self):
        self.extract_node_types()
        sections = (self.format_version(), self.format_count_section(), self.format_node_types(),
                    self.format_node_type_index(), self.format_links())
        return '/'.join(sections)


def dumps(glycan):
    return WURCSWriter(glycan).write()
