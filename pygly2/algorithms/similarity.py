import operator
from collections import defaultdict


def monosaccharide_similarity(node, target, include_substituents=True,
                              include_modifications=True, include_children=False):
    '''A heuristic for measuring similarity between monosaccharide instances'''
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
            res += check
            if check:
                node_mods.pop(node_mods.index(mod))
            qs += 1
        qs += len(node_mods)
    if include_substituents:
        node_subs = list(node for p, node in node.substituents())
        for pos, sub in target.substituents():
            check = (sub in node_subs)
            res += check
            if check:
                node_subs.pop(node_subs.index(sub))
            qs += 1
        qs += len(node_subs)
    if include_children:
        node_children = list(child for p, child in node.children())
        match_index = dict()
        for p, target_child in target.children():
            for node_child in node_children:
                c_res, c_qs = monosaccharide_similarity(
                    node_child, target_child, include_substituents, include_modifications,
                    include_children)
                match_index[target_child.id, target_child.id] = (c_res, c_qs)
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


# def build_unique_index(options_a, options_b, depth=0):
#     next_current = [()]
#     set_options_b = set(options_b)
#     partial_solutions = set()
#     while depth < len(options_a):
#         for current in next_current:
#             if len(current) > 0:
#                 current_a, current_b = map(set, (zip(*current)))
#             else:
#                 current_b = set()
#             components = set()
#             a = options_a[depth]
#             for b in set_options_b - current_b:
#                 components.add((current + ((a, b),)))
#             partial_solutions.update(components)
#         depth += 1
#         next_current = partial_solutions
#         partial_solutions = set()
#     return list(next_current)


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
