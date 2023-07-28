import re
from typing import List, NamedTuple, Tuple

from glypy.composition import molecular_graph, Composition, formula
from glypy.structure import Substituent

from .utils import WURCSFeatureNotSupported


def tokenize_alin(code: str) -> List[str]:
    STATE_ALIN = 0
    STATE_ISO = 1
    tokens = []
    current_token = ''
    i = 0
    n = len(code)
    connector_symbols = ('/', '$')
    bond_symbols = ('=', '#')
    state = STATE_ALIN
    while i < n:
        c = code[i]
        is_asterisk = c == '*'
        if c.isalpha() or is_asterisk:
            if c.isupper() or is_asterisk:
                if current_token and state == STATE_ALIN:
                    tokens.append(current_token)
                    current_token = ''
                current_token += c
                state = STATE_ALIN
        elif c in ('^', '~'):
            if current_token:
                tokens.append(current_token)
                current_token = ''
            current_token += code[i:i + 2]
            i += 2
            state = STATE_ISO
            continue
        elif c in connector_symbols:
            if current_token:
                tokens.append(current_token)
                current_token = ''
            current_token += c
        elif c in bond_symbols:
            if current_token:
                tokens.append(current_token)
                current_token = ''
            tokens.append(c)
        elif c.isdigit():
            current_token += c

        i += 1
    if current_token:
        tokens.append(current_token)
    return tokens


def parse_alin(code: str) -> molecular_graph.MolecularGraph:
    graph = molecular_graph.MolecularGraph()
    last_vertex = None
    bond_type = 1

    node_id_counter = 0
    for i, token in enumerate(tokenize_alin(code), 1):
        if token[0] in ('/', '$'):
            last_vertex = graph.vertices[int(token[1:])]
        elif token == '=':
            bond_type = 2
        elif token == '#':
            bond_type = 3
        else:
            properties = {}
            if "^" in token:
                isomorphisms = re.findall(r'[\^~][XSREZ]', token)
                token = re.sub(r'[\^~][XSREZ]', '', token)
                properties['stereoisomerism'] = isomorphisms
            if "*" in token:
                properties["is_backbone_carbon"] = True
                token = 'C'
            node_id_counter += 1
            vertex = molecular_graph.Atom(token, node_id_counter, **properties)
            graph.add_vertex(vertex)
            if last_vertex is not None:
                graph.add_edge(
                    molecular_graph.Bond([last_vertex.id, vertex.id], bond_type))
            bond_type = 1
            last_vertex = vertex
    return graph


_attachment_points_priority_pattern = re.compile(r"\*(\d+)")


def alin_to_substituent(alin: str) -> Tuple[Substituent, bool]:
    if alin.startswith("*"):
        name_query = alin.replace("*", "", 1)
    else:
        name_query = alin
    deoxy = False
    try:
        if name_query in deoxy_map_to_substituent:
            record = deoxy_map_to_substituent[name_query]
            deoxy = True
        elif name_query in map_to_substituent:
            record = map_to_substituent[name_query]
        else:
            drop_sites_name_query = name_query.replace('*', '')
            if drop_sites_name_query in deoxy_map_to_substituent:
                record = deoxy_map_to_substituent[drop_sites_name_query]
                deoxy = True
            elif drop_sites_name_query in map_to_substituent:
                record = map_to_substituent[drop_sites_name_query]
            else:
                raise KeyError(name_query)
        name = record.name
    except KeyError:
        if _attachment_points_priority_pattern.search(alin):
            raise WURCSFeatureNotSupported("ALIN with prioritized attachment points not supported")
        raise KeyError(
            "Could not locate name for substituent described by %r" % (alin,)) from None
    graph = parse_alin(alin)
    composition = Composition()
    unpaired_electrons = 2 if alin.startswith("*") else 1
    for node in graph:
        if node.data.get("is_backbone_carbon"):
            unpaired_electrons -= 1
            continue
        composition[node.element.symbol] += 1
        unpaired_electrons += node.unpaired_electrons
    composition["H"] += unpaired_electrons
    return Substituent(name, composition=composition), deoxy


def substituent_to_alin(substituent: Substituent, deoxy: bool) -> str:
    try:
        name = substituent.name
    except AttributeError:
        name = str(substituent)
    record = substituent_to_map[name]
    if not deoxy:
        return record.map_for_bmu
    else:
        return record.deoxy_map


class SubstituentTranslation(NamedTuple):
    name: str
    map_for_bmu: str
    is_carbon_swap: bool
    deoxy_map: str = None

    @property
    def formula(self):
        graph = parse_alin(self.map_for_bmu)
        return formula(graph.total_composition())

    @property
    def internal_name(self):
        subst = Substituent.internalize_name(self.name)
        return subst


_substituent_translation = [
    SubstituentTranslation('acetyl', 'OCC/3=O', None, 'CC/2=O'),
    SubstituentTranslation('amino', 'ON', None, 'N'),
    SubstituentTranslation('anhydro', '', None),
    SubstituentTranslation('bromo', None, None, 'Br'),
    SubstituentTranslation('chloro', None, None, 'Cl'),
    SubstituentTranslation('epoxy', '', None),
    SubstituentTranslation('ethanolamine', 'NCCO', False),
    SubstituentTranslation('ethyl', None, None, 'CC'),
    SubstituentTranslation('fluoro', None, None, 'F'),
    SubstituentTranslation('formyl', "OC=O", None, 'C=O'),
    SubstituentTranslation('glycolyl', None, None, 'CCO/2=O'),
    SubstituentTranslation('hydroxymethyl', 'OCO', None, 'CO',),
    SubstituentTranslation('imino', None, False, '=N'),
    SubstituentTranslation('iodo', None, None, 'I'),
    SubstituentTranslation('lactone', '', None),
    SubstituentTranslation('methyl', 'OC', None, 'C'),
    SubstituentTranslation('n-acetyl', None, None, 'NCC/3=O'),
    SubstituentTranslation('n-alanine', None, None, 'NCC^XC/4N/3=O'),
    SubstituentTranslation('n-dimethyl', None, None, 'NC/2C'),
    SubstituentTranslation('n-formyl', None, None, 'NC=O'),
    SubstituentTranslation('n-glycolyl', None, None, 'NCCO/3=O'),
    SubstituentTranslation('n-methyl', None, None, 'NC'),
    SubstituentTranslation('n-succinate', None, None, 'NCCCCO/6=O/3=O'),
    SubstituentTranslation('succinate', 'OCCCCO/6=O/3=O',
                           None, 'CCCCO/5=O/2=O'),
    SubstituentTranslation('n-sulfate', None, True, 'NSO/3=O/3=O'),
    SubstituentTranslation('n-triflouroacetyl', None,
                           None, 'NCCF/4F/4F/3=O'),
    SubstituentTranslation('nitrate', None, None, 'C=O/2=O'),
    SubstituentTranslation('phosphate', 'OPO/3O/3=O', None),
    SubstituentTranslation('pyruvate', 'OC^XO*/3CO/6=O/3C', None),
    SubstituentTranslation('pyrophosphate', 'P^XOPO/4O/4=O/2O/2=O', None),
    SubstituentTranslation('triphosphate', 'P^XOP^XOPO/6O/6=O/4O/4=O/2O/2=O', None),
    SubstituentTranslation('(r)-lactate', 'OCC^RC/4O/3=O', None),
    # SubstituentTranslation('(r)-pyruvate', '*1OC^RO*2/3CO/6=O/3C', None),
    SubstituentTranslation('(s)-lactate', 'OCC^SC/4O/3=O', None),
    SubstituentTranslation('(s)-pyruvate', None, None),
    SubstituentTranslation('sulfate', 'OSO/3=O/3=O', None),
    SubstituentTranslation('thio', "OS", None, 'S'),
    SubstituentTranslation('amidino', None, None, 'CN/2=N'),
    SubstituentTranslation('n-amidino', None, None, 'NCN/3=N'),
    SubstituentTranslation('(r)-carboxymethyl', '?*', None),
    SubstituentTranslation('(s)-carboxyethyl', 'OC^SCO/4=O/3C', None),
    SubstituentTranslation('carboxyethyl', 'CCO/3=O', None),
    SubstituentTranslation('(r)-carboxyethyl', 'OC^RCO/4=O/3C', None),
    SubstituentTranslation('(s)-carboxymethyl', 'C^SCO/3=O/2C', None),
    SubstituentTranslation('n-methyl-carbamoyl', None, None, 'CNC/2=O'),
    SubstituentTranslation('phospho-ethanolamine', 'P^XOCCN/2O/2=O', True),
    SubstituentTranslation('phospho-ethanolamine', 'OP^XOCCN/3O/3=O', True),
    SubstituentTranslation('diphospho-ethanolamine', 'P^XOP^XOCCN/4O/4=O/2O/2=O', True),
    SubstituentTranslation('phospho-choline', 'OP^XOCCNC/7C/7C/3O/3=O', None),
    SubstituentTranslation('(x)-lactate', 'OCC^XC/4O/3=O', None),
    SubstituentTranslation('(r)-1-hydroxymethyl', '?*', None),
    SubstituentTranslation('(s)-1-hydroxymethyl', '?*', None),
]


substituent_to_map = {
    unit.internal_name: unit for unit in _substituent_translation
}

deoxy_map_to_substituent = {}
map_to_substituent = {}
for unit in _substituent_translation:
    if unit.map_for_bmu is not None:
        map_to_substituent[unit.map_for_bmu] = unit
    if unit.deoxy_map is not None:
        deoxy_map_to_substituent[unit.deoxy_map] = unit


formula_to_substituent = {
    unit.formula: unit for unit in _substituent_translation
    if unit.map_for_bmu
}
