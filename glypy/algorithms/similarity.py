import operator
from collections import defaultdict
import functools

from glypy import Substituent, monosaccharides
from glypy.structure.constants import Modification, Stem


def monosaccharide_similarity(node, target, include_substituents=True,
                              include_modifications=True, include_children=False,
                              exact=True, ignore_reduction=False, visited=None,
                              short_circuit_after=None):
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
    try:
        res += (node.anomer == target.anomer) or (target.anomer.value is None)
        qs += 1
    except AttributeError:
        # Handle Substituents
        res = int(node.total_composition() == target.total_composition())
        q = 1
        return res, q

    res += (node.superclass == target.superclass) or (target.superclass.value is None)
    qs += 1
    res += (node.stem == target.stem) or (target.stem[0].value is None)
    qs += 1
    res += (node.configuration == target.configuration) or (target.configuration[0].value is None)
    qs += 1
    if short_circuit_after is not None and (res - qs) < short_circuit_after:
        return res, qs
    if include_modifications:
        node_mods = list(node.modifications.values())
        node_reduced = False
        target_reduced = False
        for mod in target.modifications.values():
            if mod == 'aldi':
                target_reduced = True
            check = (mod in node_mods)
            if check:
                if mod == 'aldi':
                    node_reduced = True
                res += 1
                node_mods.pop(node_mods.index(mod))
            qs += 1
        if ignore_reduction:
            if target_reduced:
                qs -= 1
            if node_reduced:
                res -= 1
        qs += len(node_mods) if exact else 0
    if short_circuit_after is not None and (res - qs) < short_circuit_after:
        return res, qs
    if include_substituents:
        node_subs = list(node for p, node in node.substituents())
        for pos, sub in target.substituents():
            check = (sub in node_subs)
            if check:
                res += 1
                node_subs.pop(node_subs.index(sub))
            qs += 1
        qs += len(node_subs) if exact else 0
    if short_circuit_after is not None and (res - qs) < short_circuit_after:
        return res, qs
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


def commutative_similarity(node, target, tolerance=0, *args, **kwargs):
    obs, expect = monosaccharide_similarity(node, target, *args, **kwargs)
    if (obs - expect) >= -tolerance:
        return True
    else:
        obs, expect = monosaccharide_similarity(target, node, *args, **kwargs)
        return (obs - expect) >= -tolerance


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


def has_substituent(monosaccharide, substituent):
    # Use the setter property to force the translation
    # of the name string.
    if isinstance(substituent, basestring):
        substituent = Substituent(substituent)
    substituent = substituent._name
    for position, subst in monosaccharide.substituents():
        if substituent == subst._name:
            return True
    return False


def has_modification(monosaccharide, modification):
    for position, mod in monosaccharide.modifications.items():
        if mod == modification:
            return True
    return False


def has_monosaccharide(glycan, monosaccharide, tolerance=0):
    if isinstance(monosaccharide, basestring):
        monosaccharide = monosaccharides[monosaccharide]
    visited = set()
    for node in glycan:
        if commutative_similarity(node, monosaccharide, tolerance=tolerance, visited=visited):
            return node
    return False


def is_reduced(obj):
    try:
        return obj.reducing_end is not None
    except:
        return False


def is_amine(substituent):
    if isinstance(substituent, Substituent):
        name = substituent.name
    else:
        name = substituent
    return name.startswith("n_") or name == "amino"


def is_aminated(monosaccharide):
    for p, substituent in monosaccharide.substituents():
        if is_amine(substituent):
            return True
    return False


has_fucose = functools.partial(has_monosaccharide, monosaccharide=monosaccharides["Fucose"])
has_n_acetyl = functools.partial(has_substituent, substituent=Substituent("n-acetyl"))
is_acidic = functools.partial(has_modification, modification=Modification.Acidic)
is_sulfated = functools.partial(has_substituent, substituent=Substituent("sulfate"))


def is_generic_monosaccharide(monosaccharide):
    return monosaccharide.stem[0] is Stem.x


def is_derivatized(monosaccharide):
    for pos, sub in monosaccharide.substituents():
        if sub._derivatize:
            return True
    return False
