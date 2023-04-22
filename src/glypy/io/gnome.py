"""
GNOme
-----

Tools to traverse the
`GNOme Glycan Naming and Subsumption Ontology <https://gnome.glyomics.org/>`_
"""
import re
import bisect
import warnings

from urllib.request import urlopen
from typing import Dict, DefaultDict, List, Any, Optional, Tuple, Deque
from dataclasses import dataclass, field

from lxml import etree

from glypy.structure.glycan_composition import (
    GlycanComposition, FrozenMonosaccharideResidue, HashableGlycanComposition)

from glypy.algorithms.similarity import monosaccharide_similarity
from glypy.io import glyspace
from glypy.utils import enum


class SubsumptionLevel(enum.Enum):
    molecular_weight = 1
    basecomposition = 2
    composition = 3
    topology = 4
    saccharide = 5

SubsumptionLevel.molecular_weight.add_name('molecular weight')
SubsumptionLevel.basecomposition.add_name("base composition")
SubsumptionLevel.basecomposition.add_name("base_composition")


DEFAULT_URI = "http://purl.obolibrary.org/obo/gno.owl"


generic_residues = {
    "dHex": ['Fuc'],
    "Hex": ["Gal", "Glc", "Man", "Ido",],
    "HexNAc": ['GlcNAc', 'GalNAc', 'ManNAc',],
    "Sia": ["NeuAc", "NeuGc", "Kdn"],
    "Pent": ["Xyl",],
    "HexA": ["GalA", "IdoA", "ManA", "GlcA", ],
    "HexN": ['GalN', 'ManN', 'GlcN', ],
    "Xxx": []
}

modifications = {"P", "S", "X", "Me",}


def parse_gnome_glycan(text: str) -> Dict[str, int]:
    mapping = {}
    for name, count in re.findall("([^0-9]+)(\d+)", text):
        if name.endswith("+aldi"):
            name = name[:-5]
        ct = mapping.get(name, 0)
        ct += int(count)
        mapping[name] = ct
    for generic_name, alts in generic_residues.items():
        if generic_name in mapping:
            for alt in alts:
                if alt in mapping:
                    mapping[generic_name] -= mapping[alt]
                    if mapping[generic_name] == 0:
                        mapping.pop(generic_name)
    return mapping


def _local_name(element: etree.Element) -> str:
    """Strip namespace from the XML element's name"""
    tag = element.tag
    if tag and tag[0] == "{":
        return tag.rpartition("}")[2]
    return tag


def _local_name_and_namespace(element: etree.Element) -> Tuple[str, str]:
    tag = element.tag
    if tag and tag[0] == "{":
        parts = tag.rpartition("}")
        return parts[2], parts[0][1:]
    return tag, None


OBO = "http://purl.obolibrary.org/obo/"
OBOinOWL = "http://www.geneontology.org/formats/oboInOwl#"
RDF = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"


class _GNOmeOWLXMLParser:
    property_names: Dict[str, str]
    entity_names: Dict[str, Tuple[str, Dict]]
    object_classes: Dict[str, Dict[str, Any]]
    subclasses: DefaultDict[str, List[str]]
    subsumption_levels: DefaultDict[str, List[str]]

    namespace_map = {
        "obo": "http://purl.obolibrary.org/obo/",
        "dcterms": "http://purl.org/dc/terms/",
        "owl": "http://www.w3.org/2002/07/owl#",
        "oboInOwl": "http://www.geneontology.org/formats/oboInOwl#",
        "dc": "http://purl.org/dc/elements/1.1/",
        "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
        "rdfs": "http://www.w3.org/2000/01/rdf-schema#",
    }

    def __init__(self):
        self.property_names = {}
        self.entity_names = {}
        self.object_classes = {}
        self.subclasses = DefaultDict(list)
        self.subsumption_levels = DefaultDict(list)

    def register_property(self, element: etree.Element):
        state = self.element_as_dict(element)
        self.property_names[state["about"].rpartition("/")[2]] = state["label"]

    def register_entity(self, element: etree.Element):
        state = self.element_as_dict(element)
        self.entity_names[state["about"].rpartition(
            "/")[2]] = (state["label"], state)

    def register_object_class(self, element: etree.Element):
        state = self.element_as_dict(element)
        name = state["about"].rpartition("/")[2]
        self.object_classes[name] = state
        if "subClassOf" in state:
            subclassing = state["subClassOf"]
            if isinstance(subclassing, list):
                for cls in subclassing:
                    self.subclasses[cls.split("/")[-1]].append(name)
            else:
                self.subclasses[subclassing.split("/")[-1]].append(name)
        if "has_subsumption_category" in state:
            self.subsumption_levels[state["has_subsumption_category"]].append(
                name)
        return state

    def element_as_dict(self, element: etree.Element):
        name = _local_name(element)
        info = dict()
        is_resource = False
        for k, v in element.attrib.items():
            if k and k[0] == "{":
                parts = k.rpartition("}")
                k = parts[2]
                ns = parts[0][1:]
                if ns == RDF and k == "resource":
                    is_resource = True
                    if v.startswith(OBO):
                        v_part = v.rpartition("}")[2]
                        if v_part in self.property_names:
                            v = self.property_names[v_part]
                        if v_part in self.entity_names:
                            v = self.entity_names[v_part][0]
            info[k] = v

        if element.text:
            stext = element.text.strip()
            if stext:
                if info:
                    info[name] = stext
                else:
                    return stext

        for child in element:
            name, ns = _local_name_and_namespace(child)
            if (ns == OBO or ns == OBOinOWL) and name in self.property_names:
                name = self.property_names[name]
            elif (ns == RDF or ns == OBO) and name in self.entity_names:
                name = self.entity_names[name][0]

            if name not in info:
                info[name] = self.element_as_dict(child)
            else:
                current = info[name]
                if isinstance(current, list):
                    current.append(self.element_as_dict(child))
                else:
                    info[name] = [current, self.element_as_dict(child)]
        if is_resource and len(info):
            uri = info["resource"]
            if uri.startswith(OBO):
                name = uri.rpartition("/")[2]
                if name in self.entity_names:
                    return self.entity_names[name][0]
                elif name in self.object_classes:
                    return name
            return uri
        return info

    @classmethod
    def from_element_tree(cls, tree: etree.ElementTree) -> "_GNOmeOWLXMLParser":
        self = cls()
        root = tree.getroot()
        for prop in root.iterfind("./owl:AnnotationProperty", self.namespace_map):
            self.register_property(prop)
        for ni in root.iterfind("./owl:NamedIndividual", self.namespace_map):
            self.register_entity(ni)
        for cls in root.iterfind("./owl:Class", self.namespace_map):
            self.register_object_class(cls)
        return self

    @classmethod
    def parse(cls, uri) -> "_GNOmeOWLXMLParser":
        tree = etree.parse(uri)
        return cls.from_element_tree(tree)


def strip_uri(identifier: str) -> str:
    if not identifier:
        return identifier
    if identifier.startswith("http"):
        return identifier.rsplit('/', 1)[-1]
    return identifier


@dataclass
class SubsumptionNode:
    """
    A single node in the GNOme subsumption graph, representing one or more glycans.

    Attributes
    ----------
    id : str
        A unique identifier in GNOme. Often matches the GlyTouCan accession number but
        not universally.
    definition : str
        A brief description of this node.
    subsumption_category : :class:`SubsumptionLevel`
        One of the five subsumption levels, "molecular weight", "basecomposition",
        "composition", "topology", or "saccharide" describing the level of resolution
        for the glycan in this node.
    glytoucan_id : str
        The GlyTouCan/glySpace accession number for this glycan, if it exists.
    base_composition : str
        The identifier of the "base composition" node this node is under, if any.
    composition : str
        The identifier of the "composition" node this node is under, if any.
    topology : str
        The identifier of the "topology" node this node is under, if any.
    glycan : str
        The monosaccharide composition of the glycan at this node, if any. It
        may contain generic or placeholder monosaccharides which do GNOme does
        not encode exactly or at all.
    """

    id: str
    definition: str
    subsumption_category: SubsumptionLevel
    glytoucan_id: str = field(default=None)
    base_composition: str = field(default=None)
    composition: str = field(default=None)
    topology: str = field(default=None)
    subclass_of: List[str] = field(default_factory=list)
    subclasses: List[str] = field(default_factory=list)
    glycan: Any = field(default=None)

    @classmethod
    def from_state(cls, state: dict):
        d = {}
        d['id'] = state['about'].split("/")[-1]
        d['definition'] = state['definition']
        d['subsumption_category'] = SubsumptionLevel[state['has_subsumption_category']]
        d['glytoucan_id'] = state.get('has_glytoucan_id')
        d['base_composition'] = strip_uri(state.get('has_basecomposition'))
        d['composition'] = strip_uri(state.get('has_composition'))
        d['topology'] = strip_uri(state.get('has_topology'))
        d['glycan'] = state.get('_widget_button_state')
        subclass = d['subclass_of'] = state.get('subClassOf')
        if subclass is None:
            d['subclass_of'] = []
        elif not isinstance(subclass, list):
            d['subclass_of'] = [subclass]
        return cls(**d)

    def glytoucan(self):
        if self.glytoucan_id:
            return glyspace.get(self.glytoucan_id)

    def glycan_composition(self) -> Optional[HashableGlycanComposition]:
        if self.glycan is not None:
            residue_counts = parse_gnome_glycan(self.glycan)
            gc = {}
            for k, v in residue_counts.items():
                k = FrozenMonosaccharideResidue.from_iupac_lite(k)
                if k.mass() == 0:
                    return None
                gc[k] = v
            return HashableGlycanComposition(gc)
        return None


@dataclass
class MolecularMassNode(SubsumptionNode):
    """
    A :class:`SubsumptionNode` that represents a specific mass bucket.

    This subclass is orderable by its :attr:`mass` attribute for convenient
    searching and sorting.

    Attributes
    ----------
    mass : float
        The averaged mass value of all glycans contained in this class.
    """
    mass: float = field(init=False, default=None)

    def __post_init__(self):
        self.mass = float(self.definition.split(" ")[-2])

    def __lt__(self, other):
        return self.mass < float(other)

    def __gt__(self, other):
        return self.mass > float(other)

    def __float__(self):
        return self.mass


class GNOme:
    """
    An interface for the `GNOme <https://gnome.glyomics.org/>`_ glycan subsumption graph.
    """
    terms: Dict[str, SubsumptionNode]
    subsumption_levels: Dict[str, List[str]]
    mass_index: List[SubsumptionNode]

    def __init__(self, terms: Dict[str, SubsumptionNode],
                 subsumption_levels: Dict[str, List[str]]) -> None:
        self.terms = terms
        self.subsumption_levels = subsumption_levels
        self.mass_index = self._make_mass_index()

    @classmethod
    def from_parser(cls, parser: _GNOmeOWLXMLParser):
        nodes: Dict[str, SubsumptionNode] = {}
        for k, node in parser.object_classes.items():
            if 'has_subsumption_category' not in node:
                continue
            if node['has_subsumption_category'] == "molecular weight":
                nodes[k] = MolecularMassNode.from_state(node)
            else:
                nodes[k] = SubsumptionNode.from_state(node)
            nodes[k].subclasses = parser.subclasses[k]
        return cls(nodes, parser.subsumption_levels)

    def _make_mass_index(self):
        mass_index = []
        for acc in self.subsumption_levels['molecular weight']:
            node = self.terms[acc]
            mass_index.append(node)
        mass_index.sort()
        return mass_index

    def resolve_mass(self, mass: float) -> SubsumptionNode:
        i = bisect.bisect_left(self.mass_index, mass)
        lo = self.mass_index[i - 1]
        lo_err = abs(lo.mass - mass)
        hi = self.mass_index[i]
        hi_err = abs(hi.mass - mass)
        if hi_err < lo_err:
            term = hi
        elif hi_err > lo_err:
            term = lo
        else:
            raise ValueError(
                "Ambiguous duplicate masses (%0.2f, %0.2f)" % (lo.mass, hi.mass))
        return term

    def resolve_base_composition(self, glycan_composition: GlycanComposition,
                                 node: SubsumptionNode) -> Optional[SubsumptionNode]:
        for acc in node.subclasses:
            subnode = self.terms[acc]
            if subnode.glycan:
                node_gc = parse_gnome_glycan(subnode.glycan)
                if self._match_base_composition(glycan_composition, node_gc):
                    return subnode
        return None

    def _match_monosaccharide_to_str(self, mono, mono_str) -> bool:
        if mono_str == 'Sia':
            return any(self._match_monosaccharide_to_str(mono, alt)
                       for alt in generic_residues['Sia'])
        try:
            base_mono = FrozenMonosaccharideResidue.from_iupac_lite(mono_str)
            return self._match_monosaccharide(mono, base_mono)
        except Exception as err:
            warnings.warn(f"Failed to convert {mono_str}: {err}")
            return False

    def _match_monosaccharide(self, mono_a, mono_b) -> bool:
        a, b = monosaccharide_similarity(mono_a, mono_b)
        return a == b

    def _match_base_composition(self, glycan_composition: GlycanComposition, base: dict):
        base = dict(base)
        for mono, ct in glycan_composition.items():
            for base_mono, _base_ct in base.items():
                if self._match_monosaccharide_to_str(mono, base_mono):
                    base[base_mono] -= ct
                    if base[base_mono] < 0:
                        return False
                    break
            else:
                return False
        return all(ct == 0 for ct in base.values())

    def resolve(self, glycan_composition: GlycanComposition,
                subsumption_level: SubsumptionLevel=None) -> Optional[SubsumptionNode]:
        """
        Resolve a glycan composition against the GNOme subsumption hierarchy
        to a specific level of resolution is desired.

        """
        if subsumption_level is None:
            subsumption_level = SubsumptionLevel
        node = self.resolve_mass(glycan_composition.mass())
        if node is None:
            return None
        if subsumption_level == SubsumptionLevel.molecular_weight:
            return node
        node = self.resolve_base_composition(glycan_composition, node)
        if node is None:
            return None
        if subsumption_level == SubsumptionLevel.basecomposition:
            return node
        node_queue = Deque([node])
        root = node
        possible_matches = Deque()
        while node_queue:
            node = node_queue.popleft()
            next_node = self.resolve_base_composition(glycan_composition, node)
            if next_node is None:
                possible_matches.appendleft(node)
            elif not subsumption_level or (next_node.subsumption_category <= subsumption_level):
                node_queue.append(next_node)
        if possible_matches:
            return possible_matches[0]
        return root

    @classmethod
    def load(cls, uri: Optional[str]=None):
        """
        Load the subsumption hierarchy from ``uri``, which may be a local path or a URL
        over the network to the OWL2 RDF/XML serialization of the ontology.

        If no value is passed, :const:`gnome.DEFAULT_URI` will be used.

        Because of the XML file's size, it may be preferable to store it locally if used
        frequently, possibly compressed. :mod:`lxml` can usually decompress files
        transparently.

        Parameters
        ----------
        uri : str or file-like, optional
            The location or stream to parse the GNOme OWL2 RDF/XML from. Defaults
            to :const:`gnome.DEFAULT_URI`.

        Returns
        -------
        :class:`GNOme`
        """
        if uri is None:
            uri = DEFAULT_URI
        if isinstance(uri, str) and uri.startswith("http"):
            uri = urlopen(uri)
        parser = _GNOmeOWLXMLParser.parse(uri)
        return cls.from_parser(parser)
