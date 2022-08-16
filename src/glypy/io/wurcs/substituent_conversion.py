import re
from collections import namedtuple

from glypy.composition import molecular_graph, Composition, formula
from glypy.structure import Substituent


def tokenize_alin(code):
    tokens = []
    current_token = ''
    i = 0
    n = len(code)
    connector_symbols = ('=', '/', '#', '$')
    while i < n:
        c = code[i]
        is_asterisk = c == '*'
        if c.isalpha() or is_asterisk:
            if c.isupper() or is_asterisk:
                if current_token:
                    tokens.append(current_token)
                    current_token = ''
                current_token += c
        elif c in ('^', '~'):
            current_token += code[i:i + 2]
            i += 2
            continue
        elif c in connector_symbols:
            if current_token:
                tokens.append(current_token)
                current_token = ''
            current_token += c
        elif c.isdigit():
            current_token += c

        i += 1
    if current_token:
        tokens.append(current_token)
    return tokens


def parse_alin(code):
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
                graph.add_edge(molecular_graph.Bond([last_vertex.id, vertex.id], bond_type))
            bond_type = 1
            last_vertex = vertex
    return graph


def alin_to_substituent(alin):
    if "*" in alin:
        name_query = alin.replace("*", "")
    else:
        name_query = alin
    try:
        record = map_to_substituent[name_query]
        name = record.name
    except KeyError:
        raise KeyError("Could not locate name for substituent described by %r" % (alin,))
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
    return Substituent(name, composition=composition)


def substituent_to_alin(substituent):
    try:
        name = substituent.name
    except AttributeError:
        name = str(substituent)
    record = substituent_to_map[name]
    return record.map_for_bmu


_SubstituentTranslation = namedtuple(
    "SubstituentTranslation", ("name", "map_for_bmu", "is_carbon_swap"))


class SubstituentTranslation(_SubstituentTranslation):
    @property
    def formula(self):
        graph = parse_alin(self.map_for_bmu)
        return formula(graph.total_composition())

    @property
    def internal_name(self):
        subst = Substituent.internalize_name(self.name)
        return subst


_substituent_translation = [
    SubstituentTranslation('acetyl', 'CC/2=O', None),
    SubstituentTranslation('amino', 'N', None),
    SubstituentTranslation('anhydro', '', None),
    SubstituentTranslation('bromo', 'Br', None),
    SubstituentTranslation('chloro', 'Cl', None),
    SubstituentTranslation('epoxy', '', None),
    SubstituentTranslation('ethanolamine', 'NCCO', False),
    SubstituentTranslation('ethyl', 'CC', None),
    SubstituentTranslation('fluoro', 'F', None),
    SubstituentTranslation('formyl', 'C=O', None),
    SubstituentTranslation('glycolyl', 'CCO/2=O', None),
    SubstituentTranslation('hydroxymethyl', 'CO', None),
    SubstituentTranslation('imino', '=N', False),
    SubstituentTranslation('iodo', 'I', None),
    SubstituentTranslation('lactone', '', None),
    SubstituentTranslation('methyl', 'C', None),
    SubstituentTranslation('n-acetyl', 'NCC/3=O', None),
    SubstituentTranslation('n-alanine', 'NCC^XC/4N/3=O', None),
    SubstituentTranslation('n-dimethyl', 'NC/2C', None),
    SubstituentTranslation('n-formyl', 'NC=O', None),
    SubstituentTranslation('n-glycolyl', 'NCCO/3=O', None),
    SubstituentTranslation('n-methyl', 'NC', None),
    SubstituentTranslation('n-succinate', 'NCCCCO/6=O/3=O', None),
    SubstituentTranslation('succinate', 'CCCCO/5=O/2=O', None),
    SubstituentTranslation('n-sulfate', 'NSO/3=O/3=O', True),
    SubstituentTranslation('n-triflouroacetyl', 'NCCF/4F/4F/3=O', None),
    SubstituentTranslation('nitrate', 'C=O/2=O', None),
    SubstituentTranslation('phosphate', 'OPO/3O/3=O', None),
    SubstituentTranslation('pyruvate', None, None),
    SubstituentTranslation('pyrophosphate', 'P^XOPO/4O/4=O/2O/2=O', None),
    SubstituentTranslation('triphosphate', 'P^XOP^XOPO/6O/6=O/4O/4=O/2O/2=O', None),
    SubstituentTranslation('(r)-lactate', 'CC^RC/3O/2=O', None),
    SubstituentTranslation('(r)-pyruvate', None, None),
    SubstituentTranslation('(s)-lactate', 'CC^SC/3O/2=O', None),
    SubstituentTranslation('(s)-pyruvate', None, None),
    SubstituentTranslation('sulfate', 'OSO/3=O/3=O', None),
    SubstituentTranslation('thio', 'S', None),
    SubstituentTranslation('amidino', 'CN/2=N', None),
    SubstituentTranslation('n-amidino', 'NCN/3=N', None),
    SubstituentTranslation('(r)-carboxymethyl', '?*', None),
    SubstituentTranslation('carboxymethyl', 'CCO/3=O', None),
    SubstituentTranslation('(s)-carboxymethyl', '?*', None),
    SubstituentTranslation('(r)-carboxyethyl', 'C^RCO/3=O/2C', None),
    SubstituentTranslation('(s)-carboxyethyl', 'C^SCO/3=O/2C', None),
    SubstituentTranslation('n-methyl-carbamoyl', 'CNC/2=O', None),
    SubstituentTranslation('phospho-ethanolamine', 'P^XOCCN/2O/2=O', True),
    SubstituentTranslation('phospho-ethanolamine', 'OP^XOCCN/3O/3=O', True),
    SubstituentTranslation('diphospho-ethanolamine', 'P^XOP^XOCCN/4O/4=O/2O/2=O', True),
    SubstituentTranslation('phospho-choline', 'P^XOCCN/5N/5N/2O/2=O', None),
    SubstituentTranslation('(x)-lactate', 'CC^XC/3O/2=O', None),
    SubstituentTranslation('(r)-1-hydroxymethyl', '?*', None),
    SubstituentTranslation('(s)-1-hydroxymethyl', '?*', None),
]


substituent_to_map = {
    unit.internal_name: unit for unit in _substituent_translation
}


map_to_substituent = {
    unit.map_for_bmu: unit for unit in _substituent_translation
    if unit.map_for_bmu is not None
}


formula_to_substituent = {
    unit.formula: unit for unit in _substituent_translation
    if unit.map_for_bmu
}
