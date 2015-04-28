import operator
from collections import defaultdict


def monosaccharide_similarity(node, target, include_substituents=True,
                              include_modifications=True, include_children=False,
                              exact=True, visited=None):
    '''
    A heuristic for measuring similarity between monosaccharide instances

    Compares:
        1. ring_start and ring_end
        2. superclass
        3. configuration
        4. stem
        5. anomer
        6. If `include_modifications`, each modification
        7. If `include_substituents`, each substituent
        8. If `include_children`, each child |Monosaccharide|

    Parameters
    ----------
    node: Monosaccharide
        Object to compare with
    target: Monosaccharide
        Object to compare against
    include_substituents: bool
        Include substituents in comparison (Defaults |True|)
    include_modifications: bool
        Include modifications in comparison (Defaults |True|)
    include_children: bool
        Include children in comparison (Defaults |False|)
    exact: bool
        Penalize for having unmatched attachments (Defaults |True|)

    Returns
    -------
    res: int
        Number of actual matching traits
    qs: int
        Number of expected matching traits assuming perfect equality
    '''

    if visited is None:
        visited = set()

    if (node.id, target.id) in visited:
        return 0, 0

    visited.add((node.id, target.id))

    res = 0
    qs = 0
    res += (node.anomer == target.anomer) or (target.anomer.value is None)
    qs += 1
    res += (node.superclass == target.superclass) or (target.superclass.value is None)
    qs += 1
    res += (node.stem == target.stem) or (target.stem[0].value is None)
    qs += 1
    res += (node.configuration == target.configuration) or (target.configuration[0].value is None)
    qs += 1
    if include_modifications:
        node_mods = list(node.modifications.values())
        for pos, mod in target.modifications.items():
            check = (mod in node_mods)
            if check:
                res += 1
                node_mods.pop(node_mods.index(mod))
            qs += 1
        qs += len(node_mods) if exact else 0
    if include_substituents:
        node_subs = list(node for p, node in node.substituents())
        for pos, sub in target.substituents():
            check = (sub in node_subs)
            if check:
                res += 1
                node_subs.pop(node_subs.index(sub))
            qs += 1
        qs += len(node_subs) if exact else 0
    if include_children:
        node_children = list(child for p, child in node.children())
        match_index = dict()
        for p, target_child in target.children():
            for node_child in node_children:
                c_res, c_qs = monosaccharide_similarity(
                    node_child, target_child, include_substituents=include_substituents,
                    include_modifications=include_modifications,
                    include_children=include_children, visited=visited)
                match_index[node_child.id, target_child.id] = (c_res, c_qs)

        assignments = optimal_assignment(match_index, operator.sub)
        for ix in assignments:
            a_res, a_qs = match_index[ix]
            res += a_res
            qs += a_qs

    return res, qs


def optimal_assignment(assignments, score_fn):
    '''
    Given a set of possibly overlapping matches, brute-force find the
    optimal solution. Evaluate each pairing in `assignments` with `score_fn`
    '''
    score_matrix = dict()
    for ids in assignments:
            score_matrix[ids] = score_fn(*assignments[ids])

    best_score = -float('inf')
    best_mapping = {}
    for assignment in build_unique_index_pairs(assignments):
        current_score = 0
        for ix in assignment:
            current_score += score_matrix[ix]
        if current_score > best_score:
            best_score = current_score
            best_mapping = assignment
    return best_mapping


def build_unique_index_pairs(pairs):
    '''
    Generate all unique non-overlapping sets of pairs, given in
    `pairs`
    '''
    depth = 0
    pairings = defaultdict(set)
    for a, b in pairs:
        pairings[a].add(b)
    next_current = [()]
    options_a = list(pairings)
    partial_solutions = set()
    while depth < len(options_a):
        for current in next_current:
            if len(current) > 0:
                current_a, current_b = map(set, zip(*current))
            else:
                current_b = set()
            components = set()
            a = options_a[depth]
            for b in pairings[a] - current_b:
                components.add((current + ((a, b),)))
            partial_solutions.update(components)
        depth += 1
        next_current = partial_solutions
        partial_solutions = set()
    return list(next_current)
