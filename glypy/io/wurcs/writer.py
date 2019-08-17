from collections import OrderedDict

from glypy.structure import Glycan, Monosaccharide
from glypy.structure.glycan_composition import GlycanComposition
from glypy.utils import tree

from .node_type import NodeTypeSpec
from .utils import base52


class WURCSWriter(object):
    """Implementation of WURCS encoding process.

    Includes the steps for creating each section of the WURCS encoding,
    and populating them from a saccharide structure, composition, or
    monosaccharide.

    """

    version = '2.0'

    def __init__(self, glycan):
        self.glycan = glycan
        self.node_type_map = OrderedDict()
        self.node_index_to_node_type = OrderedDict()
        self.index_to_glyph = dict()
        self.id_to_index = dict()
        self.extract_node_types()

    def extract_node_types(self):
        node_types = OrderedDict()
        node_index_to_node_type = OrderedDict()
        index_to_glyph = dict()
        id_to_index = dict()
        for i, node in enumerate(self._iter_monosaccharides(), 1):
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

    def _iter_monosaccharides(self):
        if isinstance(self.glycan, GlycanComposition):
            for key, value in self.glycan.items():
                for i in range(value):
                    yield key
        else:
            for x in self.glycan.iternodes():
                yield x

    def _iter_links(self):
        if isinstance(self.glycan, GlycanComposition):
            pass
        else:
            for x in self.glycan.iterlinks():
                yield x

    def format_count_section(self):
        count_nodes = len(list(self._iter_monosaccharides()))
        count_links = len(list(self._iter_links()))
        count_section = "%s,%s,%s" % (len(self.node_type_map), count_nodes, count_links)
        if isinstance(self.glycan, GlycanComposition):
            count_section += '+'
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
        if isinstance(self.glycan, GlycanComposition):
            return ""
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
        sections = [self.format_version(), self.format_count_section(), self.format_node_types(),
                    self.format_node_type_index(), ]
        if not isinstance(self.glycan, GlycanComposition):
            sections.append(self.format_links())
        return '/'.join(sections)


def dumps(glycan):
    """Encode a saccharide object as a WURCS 2.0 string.

    .. note::
        The WURCS canonicalization has not been implemented yet, so the generated
        string may differ from other sources. However, the :mod:`~.canonicalize` module
        can be used to standardize structures prior to encoding them.

    Parameters
    ----------
    glycan : :class:`~.Glycan`, :class:`~.GlycanComposition`, or :class:`~.Monosaccharide`
        The structure to encode

    Returns
    -------
    :class:`str`
        The structure encoded as a string.
    """
    if not isinstance(glycan, GlycanComposition):
        try:
            glycan = tree(glycan)
        except TypeError:
            if isinstance(glycan, Monosaccharide):
                nts = NodeTypeSpec.from_monosaccharide(glycan)
                return nts.to_res()
            else:
                raise
    return WURCSWriter(glycan).write()


Glycan.register_serializer('wurcs', dumps)
Monosaccharide.register_serializer('wurcs', dumps)
