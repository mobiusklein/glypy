"""Implementation of explicit sub-graph inclusion and co-traversal algorithms.

Graph inclusion here means that there exists a path through graph A that can be
mapped onto equivalent nodes in graph B that are connected in the same way.

"""
from collections import deque, defaultdict

from glypy.structure import UnknownPosition
from glypy.algorithms.similarity import commutative_similarity, commutative_similarity_score_with_tolerance
from glypy.utils import root


class TopologicalInclusionMatcher(object):
    """A recursive topological traversal.

    Attributes
    ----------
    target: :class:`~.Monosaccharide`
        The target to test for inclusion with
    reference: :class:`~.Monosaccharide`
        The reference to test for inclusion in
    substituents: bool
        Whether to examine substituents when computing similarity and inclusion.
    tolerance: float
        The magnitude of the similarity error to tolerate.
    visited: set
        The set of target and reference node id pairs that have already been considered.

    """
    def __init__(self, target, reference, substituents=True, tolerance=0, visited=None):
        self.target = target
        self.reference = reference
        self.substituents = substituents
        self.tolerance = tolerance
        self.visited = visited or set()

    @classmethod
    def compare(cls, target, reference, substituents=True, tolerance=0, visited=None):
        '''
        A generalization of :meth:`~Monosaccharide.topological_equality` which allows for ``target``
        to be matched to ``reference``, but for ``reference`` to include more. Consequently,
        this method is not commutative.

        Parameters
        ----------
        target: :class:`~.Monosaccharide`
            The monosaccharide to compare
        reference: :class:`~.Monosaccharide`
            The monosaccharide to compare against
        substituents: :class:`bool`, optional
            Whether or not to compare with the substituents of each node. Defaults
            to :const:`True`.
        tolerance: :class:`float`, optional
            The maximum difference to permit between nodes for inclusion. Defaults
            to 0

        Returns
        -------
        :class:`float`
            The inclusion score. Greater than 0 indicates topological inclusion,
            though larger corresponds to better alignment.

        See Also
        --------
        :func:`~.commutative_similarity_score_with_tolerance`
        exact_ordering_inclusion
        '''
        inst = cls(
            target, reference, substituents=substituents,
            tolerance=tolerance, visited=visited)
        score = inst.test()
        return score

    def test_similarity(self):
        """Calculate the commutative similarity between :attr:`target` and
        :attr:`reference`.

        Returns
        -------
        float

        See Also
        --------
        :func:`~.commutative_similarity_score_with_tolerance`
        """
        score, similar = commutative_similarity_score_with_tolerance(
            self.target, self.reference, self.tolerance,
            include_substituents=self.substituents)
        return score, similar

    def test_children(self):
        """Find the optimal similarity pairing between descendents of
        :attr:`target` and :attr:`reference`.

        """
        match_index = self._build_child_pair_score_map()
        required_nodes = {node.id for p, node in self.target.children()}
        optimal_assignment, score = self.optimal_assignment(match_index, required_nodes)
        return optimal_assignment, score, bool(required_nodes)

    def test(self):
        """Compute the similarity metric for the current node pair.

        Returns
        -------
        float:
            The similarity score for the current pair.
        """
        key = (self.target.id, self.reference.id)
        if key in self.visited:
            return True
        self.visited.add(key)
        pair_score, similar = self.test_similarity()
        if similar:
            child_pairs, score, had_children = self.test_children()
            if child_pairs:
                return pair_score + score
            elif not had_children:
                return pair_score
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


def exact_ordering_inclusion(target, reference, substituents=True, tolerance=0, visited=None):
    '''
    A generalization of :meth:`~Monosaccharide.exact_ordering_equality` which allows for ``target``
    to be matched to ``reference``, but for ``reference`` to include more. Consequently,
    this method is not commutative.

    Parameters
    ----------
    target: :class:`~.Monosaccharide`
        The monosaccharide to compare
    reference: :class:`~.Monosaccharide`
        The monosaccharide to compare against
    substituents: :class:`bool`, optional
        Whether or not to compare with the substituents of each node. Defaults
        to :const:`True`.
    tolerance: :class:`float`, optional
        The maximum difference to permit between nodes for inclusion. Defaults
        to 0

    Returns
    -------
    :class:`float`:
        The similarity score between the two structures. A score of 0 means no
        inclusion, and a non-zero score means inclusion, with larger scores indicating
        more node pairs matching.

    See Also
    --------
    commutative_similarity_score_with_tolerance
    topological_inclusion
    '''
    if visited is None:
        visited = set()
    if (target.id, reference.id) in visited:
        return True
    node_score, similar = commutative_similarity_score_with_tolerance(
        target, reference, tolerance, include_substituents=substituents)
    if similar:
        if substituents:
            reference_substituents = dict(reference.substituents())
            for a_pos, a_sub in target.substituents():
                b_sub = reference_substituents.get(a_pos)
                if b_sub is None:  # pragma: no cover
                    return 0
                if a_sub != b_sub:  # pragma: no cover
                    return 0
        reference_mods = dict(reference.modifications.items())
        for a_pos, a_mod in target.modifications.items():
            b_mod = reference_mods.get(a_pos)
            if b_mod is None:  # pragma: no cover
                return 0
            if a_mod != b_mod:  # pragma: no cover
                return 0
        reference_children = dict(reference.children())
        for pos, a_child in target.children():
            b_child = reference_children.get(pos)
            if b_child is None:  # pragma: no cover
                return 0
            if a_child[0] == b_child[0]:
                match_score = exact_ordering_inclusion(a_child, b_child, substituents=substituents,
                                                       tolerance=tolerance, visited=visited)
                if not match_score:
                    return 0
                node_score += match_score
            else:
                return 0
        return node_score


def subtree_of(subtree, tree, exact=False, include_substituents=True, tolerance=0):
    '''
    Test to see if `subtree` is included in `tree` anywhere. Returns the
    node id number of the first occurence of `subtree` included in `tree` or |None|
    if it is not found.

    Parameters
    ----------
    subtree: :class:`~.Glycan`
        The structure to search for. The search attempts to match the complete structure of subtree.
    tree: :class:`~.Glycan`
        The sturcture to search in. The search iterates over each residue in `tree` and calls a comparator
        function, comparing the `subtree` to the substructure rooted at that residue.
    exact: :class:`bool`
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
        if comparator(tree_root, node, substituents=include_substituents, tolerance=tolerance):
            return node.id
    return None


def find_matching_subtree_roots(subtree, tree, exact=False, include_substituents=True, tolerance=0):
    '''
    Find the list of nodes where occurences of `subtree` included in `tree` are rooted.

    Parameters
    ----------
    subtree: :class:`~.Glycan`
        The structure to search for. The search attempts to match the complete structure of subtree.
    tree: :class:`~.Glycan`
        The sturcture to search in. The search iterates over each residue in `tree` and calls a comparator
        function, comparing the `subtree` to the substructure rooted at that residue.
    exact: :class:`bool`
        If |True|, use :func:`exact_ordering_inclusion` to compare nodes. Otherwise use :func:`topological_inclusion`.
        Defaults to |False|.

    Returns
    -------
    :class:`list` of :class:`~.Monosaccharide`
    '''
    if exact:
        comparator = exact_ordering_inclusion
    else:
        comparator = topological_inclusion
    tree_root = root(subtree)
    matched_nodes = []
    for node in tree:
        if comparator(tree_root, node, substituents=include_substituents, tolerance=tolerance):
            matched_nodes.append(node)
    return matched_nodes


def walk_with(query, reference, visited=None, comparator=commutative_similarity, include_substituents=True):
    """Walk the `query` along `reference`, yielding successive matched nodes along
    a subtree of `reference` using a :term:`Exact Matching` traversal.

    This function provides slightly more detail than :func:`subtree_of`
    as it exposes every step along the subgraph traversal, rather than just the root.


    Parameters
    ----------
    query : :class:`Glycan`
        The query structure to search with
    reference : :class:`Glycan`
        The reference structure to search in
    visited : :class:`set`, optional
        The set of node id pairs to ignore
    comparator : :class:`Callable`, optional
        The :class:`Callable` object which can compare two :class:`~.Monosaccharide`

    Yields
    ------
    (:class:`~.Monosaccharide`, :class:`~.Monosaccharide`)
        Pairs of matched nodes along a path
    """
    if visited is None:
        visited = set()
    query_root = root(query)
    reference_root = root(reference)
    node_stack = deque()
    if comparator(query_root, reference_root, include_substituents=include_substituents):
        node_stack.append((reference_root, query_root))
    while len(node_stack) != 0:
        rnode, qnode = node_stack.pop()
        key = (rnode.id, qnode.id)
        if key in visited:
            continue
        visited.add(key)
        # yield the next step along the matched path
        yield rnode, qnode

        for p, rlink in rnode.links.items():
            rparent = rlink.parent
            qparent = None
            rchild = rlink.child
            qchild = None
            # query link position is known
            if p != UnknownPosition:
                # check the parent of the link for matches, potentially flowing
                # through a cycle if one is present or flowing down the tree if
                # not starting at the root
                for qlink in qnode.links[p]:
                    if comparator(qlink.parent, rparent):
                        qparent = qlink.parent
                        break
                if qparent is not None:
                    key = (rparent.id, qparent.id)
                    if key not in visited:
                        node_stack.append((rparent, qparent))
                # check the child of the link for matches, proceding down the tree
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
