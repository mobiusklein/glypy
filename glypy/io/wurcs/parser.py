import re
import string
from collections import defaultdict, namedtuple

from glypy.composition import Composition
from glypy.structure import substituent as _substituent, glycan, link as _link
from glypy.io.tree_builder_utils import try_int
from .carbon_descriptors import CarbonDescriptors
from .substituent_conversion import alin_to_substituent


_NodeTypeSpec = namedtuple("NodeTypeSpec", [
    "carbon_descriptor",
    "substituents"
])


class NodeTypeSpec(_NodeTypeSpec):
    @classmethod
    def parse(cls, text, version_number):
        if version_number < 2.0:
            raise TypeError("Cannot parse type version earlier than 2.0")
        parts = text.split("_")
        skeleton_anomer = parts[0]
        substituents = []
        if len(parts) > 1:
            ring_specification_or_substituents = parts[1]
            if re.match(r"[0-9\-1?]+-[0-9\-1?]+", ring_specification_or_substituents):
                ring_specification = ring_specification_or_substituents
            else:
                ring_specification = '?-?'
                substituents.append(ring_specification_or_substituents)
        substituents.extend(parts[2:])
        try:
            skeleton, anomer = skeleton_anomer.split("-")
        except Exception:
            skeleton = skeleton_anomer
            anomer = (-1, None)
        ring_start, ring_end = map(try_int, ring_specification.split("-"))
        anomeric_position, anomer_type = anomer
        return cls(CarbonDescriptors(skeleton, anomer_type, anomeric_position, ring_start, ring_end),
                   cls.translate_substituents(substituents))

    @classmethod
    def translate_substituents(self, substituents):
        result = []
        for map_spec in substituents:
            position = map_spec[0]
            if position == "?":
                position = -1
                i = 1
            else:
                i = 0
                position = ''
                while map_spec[i].isdigit():
                    position += map_spec[i]
                    i += 1
                position = int(position)
            map_code = map_spec[1:]
            record = alin_to_substituent(map_code)
            result.append((position, record.name))
        return result

    def to_base_type(self):
        base = self.carbon_descriptor.to_base_type()
        child_loss = Composition("H")
        for position, subst in self.substituents:
            parent_loss = _substituent.attachment_composition_info[subst]
            base.add_substituent(subst, position, parent_loss=parent_loss, child_loss=child_loss)
        return base


def base52(x):
    code = []
    if x == 0:
        return string.ascii_letters[0]
    while x > 0:
        i = x % 52
        code.append(i)
        x //= 52
    code = code[::-1]
    n = len(code)
    if n == 1:
        return ''.join([string.ascii_letters[c] for j, c in enumerate(code)])
    else:
        return ''.join([string.ascii_letters[c - 1 if j != n - 1 else c] for j, c in enumerate(code)])


class WURCSParser(object):
    def __init__(self, line, structure_class=glycan.Glycan):
        self.line = line
        self.structure_class = structure_class
        self.version = self.parse_version()
        self.node_type_count = None
        self.node_count = None
        self.edge_count = None
        self.node_type_map = {}
        self.node_index_to_node = {}
        self.glyph_to_node_index = {}

    def extract_sections(self):
        version_section, count_section, rest = self.line.split("/", 2)
        node_type_section, rest = rest.split("]/")
        node_type_section += ']'
        parts = rest.split("/")
        node_index_to_type_section = parts[0]
        node_linkage_section = parts[1]
        return (count_section, node_type_section, node_index_to_type_section, node_linkage_section)

    def parse_version(self, section=None):
        if section is None:
            section = self.line.split("/", 1)[0]
        number = float(section.split("=")[1])
        return number

    def parse_counts(self, section=None):
        if section is None:
            section = self.line.split("/", 2)[1]
        counts = (
            self.node_type_count, self.node_count,
            self.edge_count) = map(int, section.split(","))
        return counts

    def parse_node_type_section(self, section=None):
        if section is None:
            section = self.line.split("/", 2)[2].split("]/")[0] + ']'
        node_types = [s[:-1] for s in section.split("[")[1:]]
        for i, node_type in enumerate(node_types, 1):
            self.node_type_map[i] = NodeTypeSpec.parse(node_type, self.version)
        return self.node_type_map

    def parse_node_index_to_type_section(self, section=None):
        if section is None:
            section = self.extract_sections()[2]
        for i, index in enumerate(map(int, section.split('-'))):
            alpha = base52(i)
            mono = self.node_type_map[index].to_base_type()
            mono.id = i
            self.node_index_to_node[i] = mono
            self.glyph_to_node_index[alpha] = i
        return self.node_index_to_node

    parse_connection = re.compile(r"([a-zA-Z]+)([0-9]+|\?)")

    def parse_connectivity_map(self, section=None):
        if section is None:
            section = self.extract_sections()[3]
        links = section.split("_")
        for link in links:
            has_ambiguity = "|" in link
            parent_link_def, child_link_def = link.split("-")

            parent_spec = self.parse_connection.findall(parent_link_def)
            child_spec = self.parse_connection.findall(child_link_def)

            parent_glyph = parent_spec[0][0]
            child_glyph = child_spec[0][0]
            if self.glyph_to_node_index[child_glyph] < self.glyph_to_node_index[parent_glyph]:
                parent_spec, child_spec = child_spec, parent_spec

            if has_ambiguity:
                parent_nodes = []
                parent_positions = []
                for parent_glyph, parent_position in parent_spec:
                    parent_positions.append(try_int(parent_position) or -1)
                    parent_nodes.append(
                        self.node_index_to_node[self.glyph_to_node_index[parent_glyph]])

                child_nodes = []
                child_positions = []
                for child_glyph, child_position in child_spec:
                    child_positions.append(try_int(child_position) or -1)
                    child_nodes.append(
                        self.node_index_to_node[self.glyph_to_node_index[child_glyph]])
                bond = _link.AmbiguousLink(
                    parent_nodes, child_nodes, parent_position=parent_positions,
                    child_position=child_positions, parent_loss=Composition("H"),
                    child_loss=Composition("OH"))
                bond.find_open_position()
            else:
                parent_glyph, parent_position = parent_spec[0]
                child_glyph, child_position = child_spec[0]
                parent_position = try_int(parent_position) or -1
                child_position = try_int(child_position) or -1

                parent = self.node_index_to_node[self.glyph_to_node_index[parent_glyph]]
                child = self.node_index_to_node[self.glyph_to_node_index[child_glyph]]
                bond = _link.Link(
                    parent, child, parent_position=parent_position, child_position=child_position,
                    parent_loss=Composition("H"), child_loss=Composition("OH"))

    def parse(self):
        (count_section, node_type_section, node_index_to_type_section, node_linkage_section) = self.extract_sections()
        self.parse_counts(count_section)
        self.parse_node_type_section(node_type_section)
        self.parse_node_index_to_type_section(node_index_to_type_section)
        self.parse_connectivity_map(node_linkage_section)
        return self.structure_class(root=self.node_index_to_node[0], index_method='dfs', canonicalize=True)


def loads(text, structure_class=glycan.Glycan):
    parser = WURCSParser(text, structure_class=structure_class)
    structure = parser.parse()
    return structure
