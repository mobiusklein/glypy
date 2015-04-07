import logging
import itertools
from .similarity import monosaccharide_similarity
from ..utils import make_struct
from ..structure import Glycan, Monosaccharide
logger = logging.getLogger(__name__)


def substituent_inclusion(self, other):
    taken_b = set()
    b_substituents = list(other.substituents())
    cntr = 0
    for a_pos, a_substituent in self.substituents():
        matched = False
        cntr += 1
        for b_pos, b_substituent in b_substituents:
            if b_pos in taken_b:
                continue
            if b_substituent == a_substituent:
                matched = True
                taken_b.add(b_pos)
                break
        if not matched and cntr > 0:
            return False
    return True


def topological_inclusion(self, other, substituents=True, tolerance=0): 
    '''
    Performs equality testing between two monosaccharides where
    the exact ordering of child links does not have match between
    the input |Monosaccharide|s, so long as an ``a`` is included in ``b``

    Equality testing is done by :func:`similarity.monosaccharide_similarity`.

    Parameters
    ----------
    self: Glycan
        The |Glycan| to test inclusion of
    other: Glycan
        The |Glycan| to test inclusion in
    substituents: bool
        Consider substituents when comparing |Monosaccharide|s. Defaults to |True|
    tolerance: int
        The amount of error allowed when checking for flat similarity.

    Returns
    -------
    bool

    See Also
    --------
    :func:`similarity.monosaccharide_similarity`
    '''
    obs, expect = monosaccharide_similarity(self, other, include_substituents=substituents)
    if (obs - expect) >= -tolerance:
        taken_b = set()
        for a_pos, a_child in self.children():
            matched = False
            for b_pos, b_child in other.children():
                if (b_pos, b_child.id) in taken_b:
                    continue
                if topological_inclusion(a_child, b_child, substituents=substituents, tolerance=tolerance):
                    matched = True
                    taken_b.add((b_pos, b_child.id))
                    break
            if not matched and len(list(self.children())) > 0:
                return False
        return True
    return False


def exact_ordering_inclusion(self, other, substituents=True, tolerance=0):
    '''
    A generalization of :meth:`~pygly2.structure.monosaccharide.Monosaccharide.exact_ordering_equality` which
    allows for ``self`` to be matched to ``other``, but for ``other`` to include more. Consequently, this method is
    not commutative.
    '''
    obs, expect = monosaccharide_similarity(self, other, include_substituents=substituents)
    if (obs - expect) >= -tolerance:
        if substituents:
            for a_sub, b_sub in itertools.izip_longest(self.substituents(), other.substituents()):
                if b_sub is None:
                    return False
                if a_sub is None:
                    break
                if a_sub != b_sub:
                    return False
        for a_mod, b_mod in itertools.izip_longest(self.modifications.items(), other.modifications.items()):
            if b_mod is None:
                return False
            if a_mod is None:
                break
            if a_mod != b_mod:
                return False
        for a_child, b_child in itertools.izip_longest(self.children(), other.children()):
            if b_child is None:
                return False
            if a_child is None:
                break
            if a_child[0] == b_child[0]:
                if not exact_ordering_inclusion(a_child[1], b_child[1], substituents=substituents, tolerance=tolerance):
                    return False
            else:
                return False
        return True
    return False


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
    for node in tree:
        if comparator(subtree.root, node):
            return node.id
    return None


def edit_distance(seq_a, seq_b, exact=False):  # pragma: no cover
    comparator = exact_ordering_inclusion
    if exact:
        comparator = Monosaccharide.exact_ordering_equality
    previous = range(len(seq_b) + 1)
    for i, new_pos in enumerate(seq_a):
        current = [i + 1]
        for j, prev_pos in enumerate(seq_b):
            insertions = previous[j + 1] + 1
            deletions = current[j] + 1
            substitutions = previous[j] + (not comparator(new_pos, prev_pos))
            current.append(min(insertions, deletions, substitutions))
        previous = current
    return previous[-1]


def to_nested_sequence(node=None, p=None):  # pragma: no cover
    if isinstance(node, Glycan):
        node = node.root
    return [node, map(lambda x: to_nested_sequence(*x[::-1]), list(node.children()))]


def unfold(nested):  # pragma: no cover
    yield nested[0]
    for item in nested[1]:
        for k in unfold(item):
            yield k
    yield 1


def to_balanced_sequence(tree):  # pragma: no cover
    return list(unfold(to_nested_sequence(tree)))


def nested_sequence_contains(seq_a, seq_b):  # pragma: no cover
    if seq_a == seq_b:
        return True


#: Results Container for :func:`maximum_common_subgraph`
MaximumCommonSubtreeResults = make_struct(
    "MaximumCommonSubtreeResults", ("score", "tree", "similarity_matrix"))


def compare_nodes(node_a, node_b, include_substituents=True, exact=False):
    '''
    Score the pair of `node_a` and `node_b` and all pairings of their children.

    Returns
    -------
    float:
        The similarity score between the input nodes
    '''
    score = 0.
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
                               include_substituents=include_substituents, exact=exact)
    return score


def maximum_common_subgraph(seq_a, seq_b, exact=False):
    '''
    Find the maximum common subgraph between `seq_a` and `seq_b`.

    Returns
    -------
    MaximumCommonSubtreeResults
    '''
    solution_matrix = [
        [0. for i in range(len(seq_b))] for j in range(len(seq_a))]
    for i, a_node in enumerate(seq_a):
        for j, b_node in enumerate(seq_b):
            res = compare_nodes(a_node, b_node, exact=exact)
            solution_matrix[i][j] = res
    score, ix_a, ix_b = find_max(solution_matrix)
    node_a = seq_a[ix_a]
    node_b = seq_b[ix_b]
    return MaximumCommonSubtreeResults(
        score,
        extract_maximum_common_subgraph(node_a, node_b, exact=exact),
        solution_matrix)


def find_max(solution_matrix):
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


def depth(monosaccharide):
    '''
    Calculate  the distance from `monosaccharide` to its furthest grand-child node.
    '''
    d = 1
    if(len(monosaccharide.links) > 1):
        d += max(depth(ch) for p, ch in monosaccharide.children())
    return d


def extract_maximum_common_subgraph(node_a, node_b, exact=False):
    '''
    Given a pair of matched starting nodes from two separate glycan structures,
    traverse them together, copying the best matching branches into a new |Glycan|
    object.
    '''
    root = node_a.clone()
    node_stack = [(root, node_a, node_b)]
    b_taken = set()
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
            terminal = a_child.clone()
            link = [
                link for link in node_a.links[a_pos] if link.is_child(a_child)][0]
            link.clone(mcs_node, terminal)
            node_stack.append((terminal, a_child, matched_node))

    return Glycan(root)
