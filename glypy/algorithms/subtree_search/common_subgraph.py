import math
import operator
import itertools

from collections import deque, defaultdict

from glypy.algorithms.similarity import monosaccharide_similarity
from glypy.utils import make_struct, root, groupby, tree as treep
from glypy.structure import Glycan
from glypy.structure.monosaccharide import depth


#: Results Container for :func:`maximum_common_subgraph`
_MaximumCommonSubtreeResults = make_struct(
    "MaximumCommonSubtreeResults", ("score", "tree", "similarity_matrix"))

class MaximumCommonSubtreeResults(_MaximumCommonSubtreeResults):
    """Holds a maximum common subgraph solution

    Attributes
    ----------
    score: float
        The alignment score between the trees compared
    similarity_matrix: :class:`list` of :class:`list` of :class:`float`
        A simple dynamic programming solution matrix describing the alignment
        at each position.
    tree: :class:`~.Glycan`
        The maximum common subgraph, extracted as a separate :class:`~.Glycan`
        object.
    """
    __slots__ = ()


class MaximumCommonSubgraphSolver(object):
    '''
    Find the maximum common subgraph between :attr:`seq_a` and :attr:`seq_b`.

    Attributes
    ----------
    seq_a: :class:`~.Glycan`
    seq_b: :class:`~.Glycan`
    exact: bool
        Whether to use exact equality or fuzzy equality

    References
    ----------
    [1] K. F. Aoki, A. Yamaguchi, Y. Okuno, T. Akutsu, N. Ueda, M. Kanehisa, and H. Mamitsuka,
    "Efficient tree-matching methods for accurate carbohydrate database queries."
    Genome Inform. Jan. 2003.
    '''
    def __init__(self, seq_a, seq_b, exact=True):
        self.seq_a = seq_a
        self.seq_b = seq_b
        self.exact = exact
        self.solution_matrix = [
            [0. for i in range(len(seq_b))] for j in range(len(seq_a))]
        self.solution = None

        self.fit()

    def compare_nodes(self, node_a, node_b, visited=None):
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
            node_a, node_b, include_substituents=True, include_children=False)
        if self.exact:
            if observed == expected:
                score += 1
        else:
            score += observed / float(expected)
        for child_a, child_b in itertools.product((ch_a for p, ch_a in node_a.children()),
                                                  (ch_b for p, ch_b in node_b.children())):
            score += self.compare_nodes(child_a, child_b, visited=visited)
        return score

    def fit(self):
        for i, a_node in enumerate(self.seq_a):
            for j, b_node in enumerate(self.seq_b):
                res = self.compare_nodes(a_node, b_node)
                self.solution_matrix[i][j] = res
        score, ix_a, ix_b = self._find_max_of_matrix(self.solution_matrix)
        node_a = self.seq_a[ix_a]
        node_b = self.seq_b[ix_b]
        self.solution = MaximumCommonSubtreeResults(
            score,
            self._extract_maximum_common_subgraph(node_a, node_b),
            self.solution_matrix)

    def _find_max_of_matrix(self, solution_matrix):
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

    def _coordinates_for(self, solution_matrix, value):
        solutions = []
        for i in range(len(solution_matrix)):
            for j in range(len(solution_matrix[0])):
                if solution_matrix[i][j] == value:
                    solutions.append((i, j))
        return solutions

    def _extract_maximum_common_subgraph(self, node_a, node_b):
        '''
        Given a pair of matched starting nodes from two separate glycan structures,
        traverse them together, copying the best matching branches into a new |Glycan|
        object.

        Parameters
        ----------
        node_a: :class:`~.Monosaccharide`
        node_b: :class:`~.Monosaccharide`
        exact: bool
            Whether or not to take exact matches, or just take the best
            pairing. If `exact` = |True| and there is no exact match, the
            branch will terminate.

        Returns
        -------
        :class:`~.Glycan`
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
                    if self.exact and observed == expected:
                        matched_node = b_child
                        break
                    else:
                        score_pairs[b_child.id] = (expected - observed, b_child)
                if not self.exact and len(score_pairs) > 0:
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

    @classmethod
    def maximum_common_subgraph(cls, seq_a, seq_b, exact=True):
        '''
        Find the maximum common subgraph between `seq_a` and `seq_b`.

        Parameters
        ----------
        seq_a: :class:`~.Glycan`
        seq_b: :class:`~.Glycan`
        exact: bool
            Whether to use exact equality or fuzzy equality

        Returns
        -------
        :class:`MaximumCommonSubtreeResults`

        References
        ----------
        [1] K. F. Aoki, A. Yamaguchi, Y. Okuno, T. Akutsu, N. Ueda, M. Kanehisa, and H. Mamitsuka,
        "Efficient tree-matching methods for accurate carbohydrate database queries."
        Genome Inform. Jan. 2003.
        '''
        inst = cls(seq_a, seq_b, exact=exact)
        return inst.solution


maximum_common_subgraph = MaximumCommonSubgraphSolver.maximum_common_subgraph


def n_saccharide_similarity(self, other, n=2, exact=False):
    """Calculate n-saccharide similarity between two structures

    Parameters
    ----------
    self: :class:`~.Glycan`
    other: :class:`~.Glycan`
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

    References
    ----------
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
    self: :class:`~.Glycan` or list
    other : :class:`~.Glycan` or list
        :class:`~.Glycan` objects whose fragments will be compared or lists of Fragment
        to compare.
    fragmentation_parameters : :class:`dict`, optional
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


class Treelet(object):
    """Represents a subgraph of a larger :class:`~.Glycan`, with a frontier
    of node ids which are children of the current subgraph.

    Attributes
    ----------
    frontier_ids : :class:`set`
        The id values of the nodes from the parent :class:`~.Glycan` which
        are children of members of :attr:`subtree`
    subtree : :class:`~.Glycan`
        The subgraph defining the treelet
    """

    def __init__(self, subtree, frontier_ids):
        self.subtree = subtree
        self.frontier_ids = set(frontier_ids)

    @classmethod
    def from_monosaccharide(cls, monosaccharide):
        monosaccharide = root(monosaccharide)
        subtree = Glycan(
            monosaccharide.clone(prop_id=True), index_method=None)
        frontier_ids = set()
        for pos, child in monosaccharide.children():
            frontier_ids.add(child.id)
        return cls(subtree, frontier_ids)

    def __len__(self):
        return len(self.subtree)

    def __root__(self):
        return root(self.subtree)

    def __tree__(self):
        return treep(self.subtree)

    def __eq__(self, other):
        try:
            return self.subtree == other.subtree and self.frontier_ids == other.frontier_ids
        except AttributeError:
            return treep(self) == treep(other)

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash(self.subtree)

    def canonicalize(self):
        self.subtree.canonicalize()

    def expand(self, reference, frontier_id):
        node = reference.get(frontier_id)
        new_node = node.clone(prop_id=True)
        tree = self.subtree.clone(index_method=None)
        for pos, link in node.parents(True):
            new_parent = tree.get(link.parent.id)
            link.clone(new_parent, new_node)
        new_frontier = set(self.frontier_ids)
        new_frontier.remove(frontier_id)
        for pos, child in node.children():
            new_frontier.add(child.id)
        return self.__class__(tree, new_frontier)

    def expand_all(self, reference):
        extent = []
        for node_id in self.frontier_ids:
            extent.append(self.expand(reference, node_id))
        return extent


class TreeletIterator(object):
    """Iterator over all distinct :math:`k`-treelets of :attr:`tree`, all
    unique sub-trees of :attr:`tree` with :attr:`k` monosaccharides.

    Attributes
    ----------
    k : int
        The number monosaccharides per treelet
    tree : :class:`~.Glycan`
        The glycan to extract reelets from
    distinct : bool
        Whether or not to filter out duplicate treelets
    """

    def __init__(self, tree, k, distinct=True):
        self.tree = tree
        self.k = k
        self.distinct = distinct
        self.node_queue = None
        self.seen = None
        self._init_node_queue()
        self.iterator = self._make_iterator()

    def _init_node_queue(self):
        self.node_queue = deque([root(self.tree)])
        self.seen = set()

    def get_next_set(self):
        node = self.node_queue.popleft()
        tree = Treelet.from_monosaccharide(node)
        for i, child in node.children():
            self.node_queue.append(child)
        return self._extend_tree(tree, 2)

    def _make_iterator(self):
        while self.node_queue:
            treelet_set = self.get_next_set()
            for treelet in treelet_set:
                treelet.canonicalize()
                if self.distinct and treelet in self.seen:
                    continue
                self.seen.add(treelet)
                yield treelet

    def __iter__(self):
        return self.iterator

    def next(self):
        return next(self.iterator)

    def __next__(self):
        return next(self.iterator)

    def _get_next_nodes(self, start):
        child_sets = []
        for pos, child_link in self.tree.get(start.id).children(links=True):
            child_sets.append(child_link)
        return child_sets

    def _extend_tree(self, tree, n):
        extents = tree.expand_all(self.tree)
        if self.k == n:
            for t in extents:
                yield t
        else:
            for subtreelet in extents:
                for t in self._extend_tree(subtreelet, n + 1):
                    yield t


def treelets(glycan, k, distinct=True):
    """Iterator over all distinct :math:`k`-treelets of :attr:`tree`, all
    unique sub-trees of :attr:`tree` with :attr:`k` monosaccharides.

    Parameters
    ----------
    glycan : :class:`~.Glycan`
        The glycan to extract reelets from
    k : int
        The number monosaccharides per treelet
    distinct : bool
        Whether or not to filter out duplicate treelets
    """
    for treelet in TreeletIterator(glycan, k, distinct=distinct):
        yield treep(treelet)


class TreeletEnrichmentTest(object):
    '''A test to calculate the probability that the frequency of each treelet
    from glycans in :attr:`cond1` is equal to the frequency of that treelet in
    :attr:`cond2` using the Fisher's Exact Test.

    Attributes
    ----------
    cond1: list of :class:`~.Glycan`
        Glycans from the first condition, the condition to be tested for enrichment
    cond2: list of :class:`~.Glycan`
        Glycans from the second condition, the background condition
    k: int
        The size of the treelet to use
    distinct: bool
        Whether or not to count redundant treelets. Defaults to |True|
    enrichment_probabilities: dict
        Holds the enrichment p values for each treelet

    References
    ----------
    Pevzner, P., & Shamir, R. (2011). Bioinformatics for Biologists.
    New York, NY, USA: Cambridge University Press.
    '''
    def __init__(self, cond1, cond2, k, distinct=True):
        self.cond1 = list(cond1)
        self.cond2 = list(cond2)
        self.k = k
        self.distinct = distinct
        self.total = len(self.cond1) + len(self.cond2)
        self.cond1_treelets = defaultdict(int)
        self.cond2_treelets = defaultdict(int)
        self.cond1_treelet_to_glycan = defaultdict(set)
        self.cond2_treelet_to_glycan = defaultdict(set)
        self.enrichment_probabilities = dict()
        self.fit()

    def count_treelets(self):
        for g1 in self.cond1:
            for treelet in treelets(g1, self.k, distinct=self.distinct):
                self.cond1_treelets[treelet] += 1
                self.cond1_treelet_to_glycan[treelet].add(g1)
        for g2 in self.cond2:
            for treelet in treelets(g2, self.k, distinct=self.distinct):
                self.cond2_treelets[treelet] += 1
                self.cond2_treelet_to_glycan[treelet].add(g2)

    def fisher_exact_test(self, treelet):
        M = self.total
        N = len(self.cond1)
        n_pos = len(self.cond1_treelet_to_glycan[treelet])
        n_neg = len(self.cond2_treelet_to_glycan[treelet])

        fac = math.factorial
        numerator = fac(M) / float(fac(n_pos) * fac(n_neg) * fac(N - n_pos) * fac(M - N - fac(n_neg)))
        d1 = fac(M) / float(fac(n_pos + n_neg) * fac(M - n_pos - n_neg))
        d2 = fac(M) / float(fac(M - N) * fac(N) * N)
        p = numerator / (d1 * d2)
        return p

    def fit(self):
        self.count_treelets()
        for treelet in self.cond1_treelet_to_glycan:
            self.enrichment_probabilities[treelet] = self.fisher_exact_test(treelet)

    @classmethod
    def treelet_enrichment(cls, cond1, cond2, k, distinct=True):
        '''Perform a test for each treelet to determine if it is enriched in *cond1*
        compared to *cond2*.

        Parameters
        ----------
        cond1: list of :class:`~.Glycan`
            Glycans from the first condition, the condition to be tested for enrichment
        cond2: list of :class:`~.Glycan`
            Glycans from the second condition, the background condition
        k: int
            The size of the treelet to use
        distinct: bool
            Whether or not to count redundant treelets. Defaults to |True|

        Returns
        -------
        dict
            A mapping from treelet to p value from Fisher's Exact Test for
            that treelet being found with its observed frequency in *cond1*
            if it is as common in *cond1* as in *cond2*.

        References
        ----------
        Pevzner, P., & Shamir, R. (2011). Bioinformatics for Biologists.
        New York, NY, USA: Cambridge University Press.
        '''
        inst = cls(cond1, cond2, k, distinct=distinct)
        return inst.enrichment_probabilities


treelet_enrichment = TreeletEnrichmentTest.treelet_enrichment
