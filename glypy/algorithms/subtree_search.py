import itertools
import operator

from collections import deque, defaultdict

from .similarity import monosaccharide_similarity, commutative_similarity
from glypy.utils import make_struct, root, groupby
from glypy.structure import Glycan, Monosaccharide
from glypy.structure.monosaccharide import depth


class TopologicalInclusionMatcher(object):
    def __init__(self, target, reference, substituents=True, tolerance=0, visited=None):
        self.target = target
        self.reference = reference
        self.substituents = substituents
        self.tolerance = tolerance
        self.visited = visited or set()

    @classmethod
    def compare(cls, target, reference, substituents=True, tolerance=0, visited=None):
        inst = cls(
            target, reference, substituents=substituents,
            tolerance=tolerance, visited=visited)
        score = inst.test()
        return score

    def test_similarity(self):
        similar = commutative_similarity(
            self.target, self.reference, self.tolerance,
            include_substituents=self.substituents)
        return similar

    def test_children(self):
        match_index = self._build_child_pair_score_map()
        required_nodes = {node.id for p, node in self.target.children()}
        optimal_assignment, score = self.optimal_assignment(match_index, required_nodes)
        return optimal_assignment, score, bool(required_nodes)

    def test(self):
        key = (self.target.id, self.reference.id)
        if key in self.visited:
            return True
        self.visited.add(key)
        similar = self.test_similarity()
        if similar:
            child_pairs, score, had_children = self.test_children()
            if child_pairs:
                return similar + score
            elif not had_children:
                return similar
            else:
                return 0
        return 0

    def _build_child_pair_score_map(self):
        combinations = dict()
        for a_pos, a_child in self.target.children():
            for b_pos, b_child in self.reference.children():
                inclusion_score = self.compare(
                    a_child, b_child, substituents=self.substituents,
                    tolerance=self.tolerance, visited=self.visited)
                if inclusion_score > 0:
                    combinations[a_child.id, b_child.id] = inclusion_score
        return combinations

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

    def optimal_assignment(self, assignments, required_nodes=None):
        index_pair_sets = self.build_unique_index_pairs(assignments)
        best_score = -float('inf')
        best_mapping = None
        for assignment in index_pair_sets:
            current_score = 0
            has_nodes = set()
            for ix in assignment:
                current_score += assignments[ix]
                has_nodes.add(ix[0])
            if current_score > best_score and has_nodes == required_nodes:
                best_score = current_score
                best_mapping = assignment
        return best_mapping, best_score


topological_inclusion = TopologicalInclusionMatcher.compare


def exact_ordering_inclusion(self, other, substituents=True, tolerance=0, visited=None):
    '''
    A generalization of :meth:`~glypy.structure.monosaccharide.Monosaccharide.exact_ordering_equality` which
    allows for ``self`` to be matched to ``other``, but for ``other`` to include more. Consequently, this method is
    not commutative.
    '''
    if visited is None:
        visited = set()
    if (self.id, other.id) in visited:
        return True
    similar = commutative_similarity(self, other, tolerance, include_substituents=substituents)
    if similar:
        if substituents:
            other_substituents = dict(other.substituents())
            for a_pos, a_sub in self.substituents():
                b_sub = other_substituents.get(a_pos)
                if b_sub is None:  # pragma: no cover
                    return False
                if a_sub != b_sub:  # pragma: no cover
                    return False
        other_mods = dict(other.modifications.items())
        for a_pos, a_mod in self.modifications.items():
            b_mod = other_mods.get(a_pos)
            if b_mod is None:  # pragma: no cover
                return False
            if a_mod != b_mod:  # pragma: no cover
                return False
        other_children = dict(other.children())
        for pos, a_child in self.children():
            b_child = other_children.get(pos)
            if b_child is None:  # pragma: no cover
                return False
            if a_child[0] == b_child[0]:
                if not exact_ordering_inclusion(a_child, b_child, substituents=substituents,
                                                tolerance=tolerance, visited=visited):
                    return False
            else:
                return False
        return True


def subtree_of(subtree, tree, exact=False, tolerance=0):
    '''
    Test to see if `subtree` is included in `tree` anywhere. Returns the
    node id number of the first occurence of `subtree` included in `tree` or |None|
    if it is not found.

    Parameters
    ----------
    subtree: Glycan
        The structure to search for. The search attempts to match the complete structure of subtree.
    tree: Glycan
        The sturcture to search in. The search iterates over each residue in `tree` and calls a comparator
        function, comparing the `subtree` to the substructure rooted at that residue.
    exact: bool
        If |True|, use :func:`exact_ordering_inclusion` to compare nodes. Otherwise use :func:`topological_inclusion`.
        Defaults to |False|.

    Returns
    -------
    |int| or |None| if no match

    '''
    if exact:
        comparator = exact_ordering_inclusion
    else:
        comparator = topological_inclusion
    tree_root = root(subtree)
    for node in tree:
        if comparator(tree_root, node):
            return node.id
    return None


def find_matching_subtree_roots(subtree, tree, exact=False):
    if exact:
        comparator = exact_ordering_inclusion
    else:
        comparator = topological_inclusion
    tree_root = root(subtree)
    matched_nodes = []
    for node in tree:
        if comparator(tree_root, node):
            matched_nodes.append(node)
    return matched_nodes


def walk_with(query, reference, visited=None, comparator=commutative_similarity):
    if visited is None:
        visited = set()
    query_root = root(query)
    reference_root = root(reference)
    node_stack = deque([[reference_root, query_root]])

    while len(node_stack) != 0:
        rnode, qnode = node_stack.pop()
        key = (rnode.id, qnode.id)
        if key in visited:
            continue
        visited.add(key)
        yield rnode, qnode

        for p, rlink in rnode.links.items():
            rparent = rlink.parent
            qparent = None
            rchild = rlink.child
            qchild = None
            if p != -1:
                for qlink in qnode.links[p]:
                    if comparator(qlink.parent, rparent):
                        qparent = qlink.parent
                        break
                if qparent is not None:
                    key = (rparent.id, qparent.id)
                    if key not in visited:
                        node_stack.append((rparent, qparent))
                for qlink in qnode.links[p]:
                    if comparator(qlink.child, rchild):
                        qchild = qlink.child
                        break
                if qchild is not None:
                    key = (rchild.id, qchild.id)
                    if key not in visited:
                        node_stack.append((rchild, qchild))
            else:
                for qlink in qnode.links.values():
                    if comparator(qlink.parent, rparent) and comparator(qlink.child, rchild):
                        if rnode is rparent:
                            if (rchild.id, qlink.child.id) in visited:
                                continue
                        elif rnode is rchild:
                            if (rparent.id, qlink.parent.id) in visited:
                                continue
                        qparent = qlink.parent
                        qchild = qlink.child
                        break
                if qparent is not None:
                    key = (rparent.id, qparent.id)
                    if key not in visited:
                        node_stack.append((rparent, qparent))
                    key = (rchild.id, qchild.id)
                    if key not in visited:
                        node_stack.append((rchild, qchild))


#: Results Container for :func:`maximum_common_subgraph`
MaximumCommonSubtreeResults = make_struct(
    "MaximumCommonSubtreeResults", ("score", "tree", "similarity_matrix"))


def compare_nodes(node_a, node_b, include_substituents=True, exact=False, visited=None):
    '''
    Score the pair of `node_a` and `node_b` and all pairings of their children.

    Returns
    -------
    float:
        The similarity score between the input nodes
    '''
    if visited is None:
        visited = set()
    score = 0.
    if (node_a.id, node_b.id) in visited:  # pragma: no cover
        return score
    visited.add((node_a.id, node_b.id))

    observed, expected = monosaccharide_similarity(
        node_a, node_b, include_substituents=include_substituents, include_children=False)
    if exact:
        if observed == expected:
            score += 1
    else:
        score += observed / float(expected)
    for child_a, child_b in itertools.product((ch_a for p, ch_a in node_a.children()),
                                              (ch_b for p, ch_b in node_b.children())):
        score += compare_nodes(child_a, child_b,
                               include_substituents=include_substituents,
                               exact=exact, visited=visited)
    return score


def maximum_common_subgraph(seq_a, seq_b, exact=True):
    '''
    Find the maximum common subgraph between `seq_a` and `seq_b`.

    Parameters
    ----------
    seq_a: Glycan
    seq_b: Glycan
    exact: bool
        Whether to use exact equality or fuzzy equality


    Returns
    -------
    MaximumCommonSubtreeResults

    Reference
    ---------
    [1] K. F. Aoki, A. Yamaguchi, Y. Okuno, T. Akutsu, N. Ueda, M. Kanehisa, and H. Mamitsuka,
    "Efficient tree-matching methods for accurate carbohydrate database queries."
    Genome Inform. Jan. 2003.
    '''
    solution_matrix = [
        [0. for i in range(len(seq_b))] for j in range(len(seq_a))]
    for i, a_node in enumerate(seq_a):
        for j, b_node in enumerate(seq_b):
            res = compare_nodes(a_node, b_node, exact=exact)
            solution_matrix[i][j] = res
    score, ix_a, ix_b = _find_max_of_matrix(solution_matrix)
    node_a = seq_a[ix_a]
    node_b = seq_b[ix_b]
    return MaximumCommonSubtreeResults(
        score,
        _extract_maximum_common_subgraph(node_a, node_b, exact=exact),
        solution_matrix)


def _find_max_of_matrix(solution_matrix):
    '''
    Given the `solution_matrix`, find the coordinates of the maximum score

    Returns
    -------
    float:
        The maximum similarity score
    int:
        The index of the maximum similarity in `seq_a`
    int:
        The index of the maximum similarity in `seq_b`
    '''
    ix_a = ix_b = score = 0
    for i in range(len(solution_matrix)):
        for j in range(len(solution_matrix[0])):
            if solution_matrix[i][j] > score:
                score = solution_matrix[i][j]
                ix_a = i
                ix_b = j
    return score, ix_a, ix_b


def _coordinates_for(solution_matrix, value):
    solutions = []
    for i in range(len(solution_matrix)):
        for j in range(len(solution_matrix[0])):
            if solution_matrix[i][j] == value:
                solutions.append((i, j))
    return solutions


def _extract_maximum_common_subgraph(node_a, node_b, exact=False):
    '''
    Given a pair of matched starting nodes from two separate glycan structures,
    traverse them together, copying the best matching branches into a new |Glycan|
    object.

    Parameters
    ----------
    node_a: Monosaccharide
    node_b: Monosaccharide
    exact: bool
        Whether or not to take exact matches, or just take the best
        pairing. If `exact` = |True| and there is no exact match, the
        branch will terminate.
    '''
    root_ = node_a.clone()
    node_stack = [(root_, node_a, node_b)]
    b_taken = set()
    index = {node_a.id: root_}
    while len(node_stack) > 0:
        mcs_node, node_a, node_b = node_stack.pop()
        for a_pos, a_child in node_a.children():
            matched_node = None
            score_pairs = {}
            if len(node_b.links) == 0:
                continue
            for b_pos, b_child in node_b.children():
                if b_child.id in b_taken:
                    continue
                observed, expected = monosaccharide_similarity(
                    a_child, b_child, include_children=True)
                if exact and observed == expected:
                    matched_node = b_child
                    break
                else:
                    score_pairs[b_child.id] = (expected - observed, b_child)
            if not exact and len(score_pairs) > 0:
                score, contestant = min(score_pairs.values(), key=lambda x: x[0])
                cont_depth = depth(contestant)
                for diff, node in score_pairs.values():
                    if diff == score:
                        node_depth = depth(node)
                        if cont_depth < node_depth:
                            contestant = node
                            cont_depth = node_depth
                matched_node = contestant

            if matched_node is None:
                continue

            b_taken.add(matched_node.id)
            if a_child.id in index:
                terminal = index[a_child.id]
            else:
                terminal = index[a_child.id] = a_child.clone()
            link = [
                link for link in node_a.links[a_pos] if link.is_child(a_child)][0]
            link.clone(mcs_node, terminal)
            node_stack.append((terminal, a_child, matched_node))

    return Glycan(root_)


def n_saccharide_similarity(self, other, n=2, exact=False):
    """Calculate n-saccharide similarity between two structures

    Parameters
    ----------
    self: Glycan
    other: Glycan
    n: int
        Size of the fragment saccharide to consider. Defaults to *2*
    exact: bool
        Whether to use :meth:`Glycan.exact_ordering_equality` or :meth:`Glycan.topological_equality`.
        Defaults to `exact_ordering_equality`

    Returns
    -------
    score: float
        How similar these structures are at the n-saccharide level. Ranges between 0 and 1.0 where
        1.0 is exactly the same, while 0.0 means no shared n-saccharides.

    Reference
    ---------
    [1] K. F. Aoki, A. Yamaguchi, Y. Okuno, T. Akutsu, N. Ueda, M. Kanehisa, and H. Mamitsuka,
    "Efficient tree-matching methods for accurate carbohydrate database queries."
    Genome Inform. Jan. 2003.
    """
    _len = len
    _id = id

    comparator = Glycan.exact_ordering_equality
    if not exact:
        comparator = Glycan.topological_equality

    self_n_saccharides = list(subtree.tree for subtree in self.substructures(
        max_cleavages=max(self, key=operator.methodcaller('degree')).degree())
                              if _len(subtree.tree) == n)
    other_n_saccharides = list(subtree.tree for subtree in other.substructures(
        max_cleavages=max(other, key=operator.methodcaller('degree')).degree())
                               if _len(subtree.include_nodes) == n)

    n_sacch_max = max(len(self_n_saccharides), len(other_n_saccharides))
    matched = 0.
    paired = set()
    for stree in self_n_saccharides:
        for otree in other_n_saccharides:
            tid = _id(otree)
            if tid in paired:
                continue
            if comparator(stree, otree):
                matched += 1
                paired.add(tid)
    return matched / n_sacch_max


def distinct_fragments(self, other, fragmentation_parameters=None):
    """Compute the set of distinct masses observed between two structures

    Parameters
    ----------
    self: Glycan or list
    other : Glycan or list
        Glycan objects whose fragments will be compared or lists of Fragment
        to compare.
    fragmentation_parameters : dict, optional
        If `self` and `other` are |Glycan| objects, these parameters will be used
        to call :meth:`Glycan.fragments`.

    Returns
    -------
    set: The distinct fragment masses of `self`
    set: The distinct fragment masses of `other`
    """
    self_fragments = []
    other_fragments = []
    if isinstance(self, (tuple, list)):
        self_fragments = self
        other_fragments = other
    else:
        if fragmentation_parameters is None:
            fragmentation_parameters = dict(kind="BY", max_cleavages=1, average=False, charge=0)
        self_fragments = list(self.fragments(**fragmentation_parameters))
        other_fragments = list(other.fragments(**fragmentation_parameters))

    def grouper(fragment):
        return round(fragment.mass, 4)

    self_masses = set(groupby(self_fragments, grouper))
    other_masses = set(groupby(other_fragments, grouper))
    self_unique = self_masses - other_masses
    other_unique = other_masses - self_masses
    return self_unique, other_unique
