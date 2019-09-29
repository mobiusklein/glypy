import operator
from collections import defaultdict
import functools

from six import string_types as basestring

from glypy import Substituent, monosaccharides
from glypy.structure.constants import Modification, Stem, UnknownPosition


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
    ignore_reduction: bool
        Whether or not to include differences in reduction state as a
        mismatch
    ignore_ring: bool
        Whether or not to include differences in ring coordinates as
        a mismatch
    treat_null_as_wild: bool
        Whether or not to treat traits with a value of :const:`None` or
        :const:`~.UnknownPosition` as always matching when the null
        value is on the *target* residue (the residue that traits are being
        matched to).
    short_circuit_after: None or Number
        Controls whether to quit comparing nodes if the difference
        becomes too large, useful for speeding up pessimistic
        comparisons
    visited: set
        Tracks which node pairs have already been compared to break
        cycles. This carries state across multiple calls to :meth:`compare`
        and must be reset by calling :meth:`reset` before reusing an
        instance on new structures.
    '''
    def __init__(self, include_substituents=True, include_modifications=True,
                 include_children=False, exact=True, ignore_reduction=False,
                 ignore_ring=False, treat_null_as_wild=True,
                 match_attachement_positions=False, short_circuit_after=None,
                 visited=None):
        if visited is None:
            visited = set()
        self.include_substituents = include_substituents
        self.include_modifications = include_modifications
        self.include_children = include_children
        self.exact = exact
        self.ignore_ring = ignore_ring
        self.ignore_reduction = ignore_reduction
        self.treat_null_as_wild = treat_null_as_wild
        self.match_attachement_positions = match_attachement_positions
        self.visited = visited
        self.short_circuit_after = short_circuit_after

    def reset(self):
        self.visited.clear()

    @classmethod
    def similarity(cls, node, target, include_substituents=True,
                   include_modifications=True, include_children=False,
                   exact=True, ignore_reduction=False, ignore_ring=False,
                   treat_null_as_wild=True, match_attachement_positions=False,
                   short_circuit_after=None, visited=None):
        """A heuristic comparison for measuring similarity between monosaccharides.

        Compares:
            1. ring_start and ring_end
            2. superclass
            3. configuration
            4. stem
            5. anomer
            6. If `include_modifications`, each modification
            7. If `include_substituents`, each substituent
            8. If `include_children`, each child |Monosaccharide|

        The result is two numbers, the observed similarity between `node` and `target`,
        and the similarity between `target` and itself. `expected - observed` is there
        number of differences observed between the two monosaccharides, which can be useful
        for expressing how far apart two monosaccharides are in feature space. For more distant
        similarity testing, especially when considering children, the ratio `observed / expected`
        might be used instead.

        Similarity is not symmetric, e.g. ``a -> b != b -> a``. A commutative version of similarity
        can be used by calculating both directions, and taking the result with the smallest error.

        Parameters
        ----------
        node : :class:`~.Monosaccharide`
            The reference monosaccharide
        target : :class:`~.Monosaccharide`
            The monosaccharide to compare against
        include_substituents: bool
            Include substituents in comparison (Defaults |True|)
        include_modifications: bool
            Include modifications in comparison (Defaults |True|)
        include_children: bool
            Include children in comparison (Defaults |False|)
        exact: bool
            Penalize for having unmatched attachments (Defaults |True|)
        ignore_reduction: bool
            Whether or not to include differences in reduction state as a
            mismatch
        ignore_ring: bool
            Whether or not to include differences in ring coordinates as
            a mismatch
        treat_null_as_wild: bool
            Whether or not to treat traits with a value of :const:`None` or
            :const:`~.UnknownPosition` as always matching when the null
            value is on the *target* residue (the residue that traits are being
            matched to).
        short_circuit_after: None or Number
            Controls whether to quit comparing nodes if the difference
            becomes too large, useful for speeding up pessimistic
            comparisons
        visited: set
            Tracks which node pairs have already been compared to break
            cycles. This carries state across multiple calls to :meth:`compare`
            and must be reset by calling :meth:`reset` before reusing an
            instance on new structures.

        Returns
        -------
        :class:`int`: observed
            The number of observed features that matched
        :class:`int`: expected
            The number of features that could have been matched
        """
        inst = cls(
            include_substituents=include_substituents,
            include_modifications=include_modifications,
            include_children=include_children,
            exact=exact, ignore_reduction=ignore_reduction,
            ignore_ring=ignore_ring, treat_null_as_wild=treat_null_as_wild,
            match_attachement_positions=match_attachement_positions,
            short_circuit_after=short_circuit_after,
            visited=visited)
        return inst.compare(node, target)

    def compare_anomer(self, node, target):
        """Compare :attr:`~.Monosaccharide.anomer` of `node` and `target`

        Parameters
        ----------
        node : :class:`~.Monosaccharide`
            The reference monosaccharide
        target : :class:`~.Monosaccharide`
            The monosaccharide to compare against

        Returns
        -------
        test: int
            The similarity between node and target
        reference: int
            The expected similarity between target and itself
        """
        test = (node.anomer == target.anomer) or ((target.anomer.value is None) and self.treat_null_as_wild)
        reference = 1
        return int(test), reference

    def compare_compositions(self, node, target):
        """Compre :meth:`~.Monosaccharide.total_composition` of `node` and `target`.

        Parameters
        ----------
        node : :class:`~.Monosaccharide`
            The reference monosaccharide
        target : :class:`~.Monosaccharide`
            The monosaccharide to compare against

        Returns
        -------
        test: int
            The similarity between node and target
        reference: int
            The expected similarity between target and itself
        """
        test = int(node.total_composition() == target.total_composition())
        reference = 1
        return test, reference

    def compare_ring_structure(self, node, target):
        """Compare the :attr:`~.Monosaccharide.superclass`, :attr:`~.Monosaccharide.stem`,
        :attr:`~.Monosaccharide.configuration`, and possibly :attr:`~.Monosaccharide.ring_start`
        and :attr:`~.Monosaccharide.ring_end` of `node` and `target`.

        If :attr:`ignore_ring` is :const:`True`, ring positions will be assumed to match.

        Parameters
        ----------
        node : :class:`~.Monosaccharide`
            The reference monosaccharide
        target : :class:`~.Monosaccharide`
            The monosaccharide to compare against

        Returns
        -------
        test: int
            The similarity between node and target
        reference: int
            The expected similarity between target and itself
        """
        test = reference = 0
        test += (node.superclass == target.superclass) or ((target.superclass.value is None) and
                                                           self.treat_null_as_wild)
        reference += 1
        test += (node.stem == target.stem) or ((target.stem[0].value is None) and self.treat_null_as_wild)
        reference += 1
        test += (node.configuration == target.configuration) or ((target.configuration[0].value is None) and
                                                                 self.treat_null_as_wild)
        reference += 1
        if not self.ignore_ring:
            test += (node.ring_start == target.ring_start) or ((target.ring_start == UnknownPosition) and
                                                               self.treat_null_as_wild)
            reference += 1
            test += (node.ring_end == target.ring_end) or ((target.ring_end == UnknownPosition) and
                                                           self.treat_null_as_wild)
            reference += 1
        else:
            test += 2
            reference += 2
        return test, reference

    def compare_modifications(self, node, target):
        """Compare the modifications of `node` and `target`

        Parameters
        ----------
        node : :class:`~.Monosaccharide`
            The reference monosaccharide
        target : :class:`~.Monosaccharide`
            The monosaccharide to compare against

        Returns
        -------
        test: int
            The similarity between node and target
        reference: int
            The expected similarity between target and itself
        """
        test = reference = 0
        node_reduced = False
        target_reduced = False
        n_mods = 0
        if self.match_attachement_positions:
            node_mods = node.modifications
            n_mods = len(node_mods)
            for pos, mod in target.modifications.items():
                if mod == 'aldi':
                    target_reduced = True
                check = (mod in node_mods[pos])
                if check:
                    if mod == 'aldi':
                        node_reduced = True
                    test += 1
                    n_mods -= 1
                reference += 1
        else:
            node_mods = list(node.modifications.values())
            n_mods = len(node_mods)
            for mod in target.modifications.values():
                if mod == 'aldi':
                    target_reduced = True
                check = (mod in node_mods)
                if check:
                    if mod == 'aldi':
                        node_reduced = True
                    test += 1
                    node_mods.pop(node_mods.index(mod))
                    n_mods -= 1
                reference += 1

        if self.ignore_reduction:
            if target_reduced:
                reference -= 1
            if node_reduced:
                test -= 1
        reference += n_mods if self.exact else 0
        return test, reference

    def compare_substituents(self, node, target):
        """Compare the substituents of `node` and `target`.

        If :attr:`match_attachement_positions` is :const:`True`,
        the positions must match exactly, otherwise only the substituent
        types must match.

        Parameters
        ----------
        node : :class:`~.Monosaccharide`
            The reference monosaccharide
        target : :class:`~.Monosaccharide`
            The monosaccharide to compare against

        Returns
        -------
        test: int
            The similarity between node and target
        reference: int
            The expected similarity between target and itself
        """
        test = reference = 0
        if self.match_attachement_positions:
            node_subs = defaultdict(list)
            n_subs = 0
            for p, sub in node.substituents():
                n_subs += 1
                node_subs[p].append(sub)
            for pos, sub in target.substituents():
                if (sub in node_subs[pos]):
                    test += 1
                    n_subs -= 1
                reference += 1
        else:
            node_subs = [node_sub for p, node_sub in node.substituents()]
            n_subs = len(node_subs)
            for pos, sub in target.substituents():
                if (sub in node_subs):
                    test += 1
                    node_subs.pop(node_subs.index(sub))
                    n_subs -= 1
                reference += 1
        reference += n_subs if self.exact else 0
        return test, reference

    def _build_child_pair_score_map(self, node, target):
        '''Compute all pair-wise similarities between children
        of ``node`` and ``target``

        Parameters
        ----------
        node : :class:`~.Monosaccharide`
            The reference monosaccharide
        target : :class:`~.Monosaccharide`
            The monosaccharide to compare against

        Returns
        -------
        dict
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
        """Compute the similarity of the set of children between `node` and `target`

        If :attr:`match_attachement_positions` is :const:`True`, this will require
        the positions of child nodes to match exactly, otherwise, all pair-wise combinations
        will be considered and the optimal solution will be selected using :meth:`optimal_assignment`

        Parameters
        ----------
        node : :class:`~.Monosaccharide`
            The reference monosaccharide
        target : :class:`~.Monosaccharide`
            The monosaccharide to compare against

        Returns
        -------
        test: int
            The similarity between node and target
        reference: int
            The expected similarity between target and itself
        """
        test = reference = 0
        if self.match_attachement_positions:
            node_children = defaultdict(list)
            for pos, child in node.children():
                node_children[pos].append(child)
            target_children = defaultdict(list)
            for pos, child in target.children():
                target_children[pos].append(child)
            for pos, children in target_children.items():
                for t_child, n_child in zip(children, node_children[pos]):
                    test_child, reference_child = self.compare(t_child, n_child)
                    test += test_child
                    reference += reference_child
        else:
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
        """Calculate the similarity between `node` and `target`.

        This method does most of the organizational work, calling the appropriate
        methods and checking for short-circuiting.

        Parameters
        ----------
        node : :class:`~.Monosaccharide`
            The reference monosaccharide
        target : :class:`~.Monosaccharide`
            The monosaccharide to compare against

        Returns
        -------
        test: int
            The similarity between node and target
        reference: int
            The expected similarity between target and itself
        """
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
    """Apply :func:`monosaccharide_similarity` to ``node`` and ``target`` for both
    ``node --> target`` and ``target --> node``, returning whether either comparison
    passes the tolerance threshold.

    Parameters
    ----------
    node: :class:`~.Monosaccharide`
        The reference monosaccharide
    target: :class:`~.Monosaccharide`
        The monosaccharide to compare against
    tolerance: :class:`int`, optional
        The minimum number of errors to tolerate
    *args:
        Forwarded to :func:`monosaccharide_similarity`
    **kwargs:
        Forwarded to :func:`monosaccharide_similarity`

    Returns
    -------
    :class:`bool`
    """
    obs, expect = monosaccharide_similarity(node, target, *args, **kwargs)
    if (expect - obs) <= tolerance:
        return True
    else:
        obs, expect = monosaccharide_similarity(target, node, *args, **kwargs)
        return (expect - obs) <= tolerance


def commutative_similarity_score(node, target, *args, **kwargs):
    """Apply :func:`monosaccharide_similarity` to ``node`` and ``target`` for both
    ``node --> target`` and ``target --> node``, returning the maximally normalized
    ratio of `observed / expected`.

    Parameters
    ----------
    node: :class:`~.Monosaccharide`
        The reference monosaccharide
    target: :class:`~.Monosaccharide`
        The monosaccharide to compare against
    *args:
        Forwarded to :func:`monosaccharide_similarity`
    **kwargs:
        Forwarded to :func:`monosaccharide_similarity`

    Returns
    -------
    :class:`float`:
        The maximal similarity score ratio
    """
    a_b, b_b = monosaccharide_similarity(node, target, *args, **kwargs)
    b_a, a_a = monosaccharide_similarity(target, node, *args, **kwargs)
    return max(a_b / (1. * b_b), b_a / (1. * a_a))


def commutative_similarity_score_with_tolerance(node, target, tolerance, *args, **kwargs):
    """Apply :func:`monosaccharide_similarity` to ``node`` and ``target`` for both
    ``node --> target`` and ``target --> node``, returning the maximally normalized
    ratio score, and whether there was a pair error less than `tolerance`.

    This can be viewed as a combination of :func:`commutative_similarity` and
    :func:`commutative_similarity_score` while making fewer calls to
    :func:`monosaccharide_similarity`.

    Parameters
    ----------
    node: :class:`~.Monosaccharide`
        The reference monosaccharide
    target: :class:`~.Monosaccharide`
        The monosaccharide to compare against
    tolerance: :class:`int`
        The minimum number of errors to tolerate
    *args:
        Forwarded to :func:`monosaccharide_similarity`
    **kwargs:
        Forwarded to :func:`monosaccharide_similarity`

    Returns
    -------
    :class:`float`:
        The maximal similarity score ratio
    :class:`bool`:
        Whether the difference passes error tolerance
    """
    a_b, b_b = monosaccharide_similarity(node, target, *args, **kwargs)
    b_a, a_a = monosaccharide_similarity(target, node, *args, **kwargs)
    pass_threshold = False
    if (b_b - a_b) <= tolerance or (a_a - b_a) <= tolerance:
        pass_threshold = True
    return max(a_b / (1. * b_b), b_a / (1. * a_a)), pass_threshold



def has_substituent(monosaccharide, substituent):
    """Checks whether ``monosaccharide`` has any substituent groups
    matching ``substituent``.

    Parameters
    ----------
    monosaccharide : :class:`~.Monosaccharide`
        The monosaccharide to check
    substituent : :class:`~.Substituent` or :class:`str`
        The substituent to check for

    Returns
    -------
    :class:`bool`
    """
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
    """Checks whether ``monosaccharide`` has any modification sites
    matching ``modification``.

    Parameters
    ----------
    monosaccharide : :class:`~.Monosaccharide`
        The monosaccharide to check
    modification : :class:`~.constants.Modification` or :class:`str`
        The modification to check for

    Returns
    -------
    :class:`bool`
    """
    for position, mod in monosaccharide.modifications.items():
        if mod == modification:
            return True
    return False


def has_monosaccharide(glycan, monosaccharide, tolerance=0, *args, **kwargs):
    """Checks whether ``glycan`` has any monosaccharide nodes
    matching ``monosaccharide`` within ``tolerance`` using
    :func:`commutative_similarity`

    Parameters
    ----------
    glycan: :class:`~.SaccharideCollection`
        The glycan structure or composition to search
    monosaccharide: :class:`~.Monosaccharide`
        The monosaccharide to search for
    tolerance : int, optional
        The error tolerance to use
    *args:
        Forwarded to :func:`monosaccharide_similarity`
    **kwargs:
        Forwarded to :func:`monosaccharide_similarity`

    Returns
    -------
    :class:`bool`
    """
    if isinstance(monosaccharide, basestring):
        monosaccharide = monosaccharides[monosaccharide]
    visited = set()
    for node in glycan:
        if commutative_similarity(
                node, monosaccharide, tolerance=tolerance, visited=visited, *args, **kwargs):
            return node
    return False


def is_reduced(obj):
    """A simple predicate to test whether an object has a reduced structure.

     If `obj` does not have a `reducing_end` attribute, this will return :const:`False`

    Parameters
    ----------
    obj : object
        The object to check

    Returns
    -------
    bool
    """
    try:
        return obj.reducing_end is not None
    except AttributeError:
        return False


def is_amine(substituent):
    """A simple predicate to test whether a substituent has an amine group adjacent
    attached to the carbon backbone by naming convention.

    This predicate checks to see if the name of the substituent is "amino" or if it
    starts with the phrase "n_" only.

    Parameters
    ----------
    substituent : Substituent  or str
        The object to test

    Returns
    -------
    bool
    """
    if isinstance(substituent, Substituent):
        name = substituent.name
    else:
        name = substituent
    return name.startswith("n_") or name == "amino"


def is_aminated(monosaccharide):
    """Tests to see if any substituents of `monosaccharide` are amines.

    Each substituent is tested using :func:`is_amine`, with all the caveats that
    entails.

    Parameters
    ----------
    monosaccharide : Monosaccharide
        The monosaccharide to test

    Returns
    -------
    bool

    See Also
    --------
    is_amine
    """
    for p, substituent in monosaccharide.substituents():
        if is_amine(substituent):
            return True
    return False


has_fucose = functools.partial(has_monosaccharide, monosaccharide=monosaccharides["Fucose"])
has_n_acetyl = functools.partial(has_substituent, substituent=Substituent("n-acetyl"))
is_acidic = functools.partial(has_modification, modification=Modification.Acidic)
is_sulfated = functools.partial(has_substituent, substituent=Substituent("sulfate"))


def is_generic_monosaccharide(monosaccharide):
    """Tests if the :attr:`~.Monosaccharide.stem` is unknown.

    Parameters
    ----------
    monosaccharide : Monosaccharide
        The object to test

    Returns
    -------
    bool
    """
    return monosaccharide.stem[0] is Stem.x


def is_derivatized(monosaccharide):
    """Tests whether any of the substituents attached to `monosaccharide` were
    added by derivatization.

    Parameters
    ----------
    monosaccharide : Monosaccharide
        The object to test

    Returns
    -------
    bool
    """
    for pos, sub in monosaccharide.substituents():
        if sub._derivatize:
            return True
    return False
