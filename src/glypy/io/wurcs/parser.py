import re
try:
    from urllib import unquote
except ImportError:
    from urllib.parse import unquote

from glypy.composition import Composition
from glypy.structure import glycan, link as _link, glycan_composition
from glypy.io.tree_builder_utils import try_int

from .node_type import NodeTypeSpec
from .utils import base52, WURCSFeatureNotSupported


class WURCSParser(object):
    def __init__(self, line, structure_class=glycan.Glycan):
        self.line = unquote(line)
        self.structure_class = structure_class
        self.version = self.parse_version()
        self.node_type_count = None
        self.node_count = None
        self.edge_count = None
        self.node_type_map = {}
        self.node_index_to_node = {}
        self.glyph_to_node_index = {}
        self.has_uncertain_linkages = False

    def extract_sections(self):
        version_section, count_section, rest = self.line.split("/", 2)
        node_type_section, rest = rest.split("]/")
        node_type_section += ']'
        parts = rest.split("/", 1)
        node_index_to_type_section = parts[0]
        if len(parts) == 2:
            rest = parts[1]
        else:
            rest = ''
        node_linkage_section = rest
        return (count_section, node_type_section, node_index_to_type_section, node_linkage_section)

    def parse_version(self, section=None):
        if section is None:
            section = self.line.split("/", 1)[0]
        number = float(section.split("=")[1])
        return number

    def parse_counts(self, section=None):
        if section is None:
            section = self.line.split("/", 2)[1]
        if "+" in section:
            self.has_uncertain_linkages = True
        counts = (
            self.node_type_count, self.node_count,
            self.edge_count) = list(map(lambda x: int(x.replace("+", "")), section.split(",")))
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
            mono = self.node_type_map[index].to_monosaccharide()
            mono.id = i
            self.node_index_to_node[i] = mono
            self.glyph_to_node_index[alpha] = i
        return self.node_index_to_node

    parse_connection = re.compile(r"([a-zA-Z]+)([0-9]+|\?)")

    def parse_connectivity_map(self, section=None):
        if section is None:
            section = self.extract_sections()[3]
        if "{" in section or "}" in section:
            links = section.split("_")
            if len(links) > 1 and len(set(links)) == 1:
                # This is a composition, everybody is ambiguously linked to everybody
                return False
            raise WURCSFeatureNotSupported("Braced Undefined Linkages are not supported")

        links = section.split("_")
        for link in links:
            has_ambiguity = "|" in link
            has_bridge = "*" in link
            if has_bridge:
                raise WURCSFeatureNotSupported("Bridging MAPs are not supported.")

            parent_link_def, child_link_def = link.split("-", 1)

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
        return True

    def _to_composition(self):
        gc = glycan_composition.GlycanComposition()
        for node in self.node_index_to_node.values():
            gc[glycan_composition.MonosaccharideResidue.from_monosaccharide(node)] += 1
        return gc

    def parse(self):
        (count_section, node_type_section, node_index_to_type_section, node_linkage_section) = self.extract_sections()
        self.parse_counts(count_section)
        self.parse_node_type_section(node_type_section)
        self.parse_node_index_to_type_section(node_index_to_type_section)
        if node_linkage_section:
            if self.parse_connectivity_map(node_linkage_section):
                return self.structure_class(root=self.node_index_to_node[0], index_method='dfs', canonicalize=True)
            return self._to_composition()
        else:
            return self._to_composition()


def loads(text, structure_class=glycan.Glycan):
    """Parse a WURCS-encoded glycan structure from `text` into a :class:`~.Glycan`
    or :class:`~.GlycanComposition`.

    Parameters
    ----------
    text : str
        The WURCS string to parse
    structure_class : :class:`type`, optional
        The class to use to wrap the :class:`~.Monosaccharide` graph (the default is :class:`~.Glycan`)

    Returns
    -------
    :class:`~.Glycan` or :class:`~.GlycanComposition`
        The parsed result
    """
    parser = WURCSParser(text, structure_class=structure_class)
    structure = parser.parse()
    return structure
