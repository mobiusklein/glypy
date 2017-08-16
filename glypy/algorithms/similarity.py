import operator
from collections import defaultdict
import functools

from six import string_types as basestring

from glypy import Substituent, monosaccharides
from glypy.structure.constants import Modification, Stem


class NodeSimilarityComparator(object):
    '''A heuristic comparison for measuring similarity between monosaccharides.

    Compares:
        1. ring_start and ring_end
        2. superclass
        3. configuration
        4. stem
        5. anomer
        6. If `include_modifications`, each modification
        7. If `include_substituents`, each substituent
        8. If `include_children`, each child |Monosaccharide|

    Attributes
    ----------
    include_substituents: bool
        Include substituents in comparison (Defaults |True|)
    include_modifications: bool
        Include modifications in comparison (Defaults |True|)
    include_children: bool
        Include children in comparison (Defaults |False|)
    exact: bool
        Penalize for having unmatched attachments (Defaults |True|)
    short_circuit_after: None or Number
        Controls whether to quit comparing nodes if the difference
        becomes too large, useful for speeding up pessimistic
        comparisons
    '''
    def __init__(self, include_substituents=True, include_modifications=True,
                 include_children=False, exact=True, ignore_reduction=False,
                 short_circuit_after=None, visited=None):
        if visited is None:
            visited = set()
        self.include_substituents = include_substituents
        self.include_modifications = include_modifications
        self.include_children = include_children
        self.exact = exact
        self.ignore_reduction = ignore_reduction
        self.visited = visited
        self.short_circuit_after = short_circuit_after
    
    def reset(self):
        self.visited.clear()

    @classmethod
    def similarity(cls, node, target, include_substituents=True,
                   include_modifications=True, include_children=False,
                   exact=True, ignore_reduction=False,
                   short_circuit_after=None, visited=None):
        inst = cls(
            include_substituents=include_substituents,
            include_modifications=include_modifications,
            include_children=include_children,
            exact=exact, ignore_reduction=ignore_reduction,
            short_circuit_after=short_circuit_after,
            visited=visited)
        return inst.compare(node, target)
    
    def compare_anomer(self, node, target):
        test = (node.anomer == target.anomer) or (
                target.anomer.value is None)
        reference = 1
        return int(test), reference

    def compare_compositions(self, node, target):
        test = int(node.total_composition() == target.total_composition())
        reference = 1
        return test, reference
    
    def compare_ring_structure(self, node, target):
        test = reference = 0
        test += (node.superclass == target.superclass) or (target.superclass.value is None)
        reference += 1
        test += (node.stem == target.stem) or (target.stem[0].value is None)
        reference += 1
        test += (node.configuration == target.configuration) or (target.configuration[0].value is None)
        reference += 1
        return test, reference
    
    def compare_modifications(self, node, target):
        test = reference = 0
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
                test += 1
                node_mods.pop(node_mods.index(mod))
            reference += 1
        if self.ignore_reduction:
            if target_reduced:
                reference -= 1
            if node_reduced:
                test -= 1
        reference += len(node_mods) if self.exact else 0
        return test, reference
    
    def compare_substituents(self, node, target):
        test = reference = 0
        node_subs = list(node for p, node in node.substituents())
        for pos, sub in target.substituents():
            check = (sub in node_subs)
            if check:
                test += 1
                node_subs.pop(node_subs.index(sub))
            reference += 1
        reference += len(node_subs) if self.exact else 0
        return test, reference
    
    def _build_child_pair_score_map(self, node, target):
        '''Compute all pair-wise similarities between children
        of ``node`` and ``target``
        '''
        # TODO: support exactness penalty here, maybe recursively?
        node_children = list(child for p, child in node.children())
        match_index = dict()
        for p, target_child in target.children():
            for node_child in node_children:
                c_res, c_qs = self.compare(node_child, target_child)
                match_index[node_child.id, target_child.id] = (c_res, c_qs)
        return match_index

    def compare_children(self, node, target):
        test = reference = 0
        match_index = self._build_child_pair_score_map(node, target)

        assignments = self.optimal_assignment(match_index)
        for ix in assignments:
            a_test, a_reference = match_index[ix]
            test += a_test
            reference += a_reference
        return test, reference
    
    def _check_short_circuit(self, test, reference):
        return self.short_circuit_after is not None and (test - reference) < self.short_circuit_after

    def compare(self, node, target):
        key = (node.id, target.id)
        if key in self.visited:
            return 0, 0
        self.visited.add(key)
        test = 0
        reference = 0
        try:
            t, r = self.compare_anomer(node, target)
            test += t
            reference += r
        except AttributeError:
            # must be handling substituents
            t, r = self.compare_compositions(node, target)
            test += t
            reference += r
            return test, reference
        t, r = self.compare_ring_structure(node, target)
        test += t
        reference += r
        if self._check_short_circuit(test, reference):
            return test, reference
        if self.include_modifications:
            t, r = self.compare_modifications(node, target)
            test += t
            reference += r
            if self._check_short_circuit(test, reference):
                return test, reference
        if self.include_substituents:
            t, r = self.compare_substituents(node, target)
            test += t
            reference += r
            if self._check_short_circuit(test, reference):
                return test, reference
        if self.include_children:
            t, r = self.compare_children(node, target)
            test += t
            reference += r
        return test, reference

    def optimal_assignment(self, assignments):
        '''
        Given a set of possibly overlapping matches, find the
        optimal solution.
        '''
        diff_map = dict()
        for ids in assignments:
            diff_map[ids] = operator.sub(*assignments[ids])

        index_pair_sets = self.build_unique_index_pairs(assignments)

        def ordering_fn(index_pairs):
            total = 0
            for ix in index_pairs:
                total += max(assignments[ix])
            return total

        ordered_index_pairs_sets = sorted(index_pair_sets, key=ordering_fn, reverse=True)

        best_score = -float('inf')
        best_mapping = {}
        # prefer solutions which have a higher maximum value (more points of comparison)
        for assignment in ordered_index_pairs_sets:
            current_score = 0
            for ix in assignment:
                current_score += diff_map[ix]
            if current_score > best_score:
                best_score = current_score
                best_mapping = assignment
        return best_mapping

    def build_unique_index_pairs(self, pairs):
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


monosaccharide_similarity = NodeSimilarityComparator.similarity


def commutative_similarity(node, target, tolerance=0, *args, **kwargs):
    obs, expect = monosaccharide_similarity(node, target, *args, **kwargs)
    if (obs - expect) >= -tolerance:
        return True
    else:
        obs, expect = monosaccharide_similarity(target, node, *args, **kwargs)
        return (obs - expect) >= -tolerance


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
    except AttributeError:
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
