from collections import defaultdict
import operator

from six import string_types as basestring

from ...utils import groupby
from ...structure import (
    named_structures, Monosaccharide, Substituent,
    Anomer, Stem, RingType,
    SuperClass, Configuration)
from ...algorithms.similarity import (monosaccharide_similarity, has_substituent,
                                      has_modification, has_monosaccharide,
                                      is_generic_monosaccharide)
from ...composition.composition_transform import strip_derivatization
from .synonyms import monosaccharides as monosaccharide_synonyms


def has_ambiguity(node):
    ambiguous = node.stem[0] is Stem.x or node.anomer is Anomer.x or\
        node.superclass is SuperClass.x or node.ring_type is RingType.x
    return ambiguous


# A static copy of monosaccharide names to structures for copy-free comparison
monosaccharides = dict(named_structures.monosaccharides)
monosaccharides_ordered = sorted(list(monosaccharides.items()), key=lambda x: has_ambiguity(x[1]))


def get_preferred_name(name, selector=min, key=len):
    '''
    Given a name, of its synonyms, find the name that satisfies the `selector`
    criterion function (:func:`min`) based on some `key` function of the name (:func:`len`)

    Parameters
    ----------
    name: str
        Given name to compare to synonyms
    selector: function
        Function to use to select the preferred name by some statistic
    key: function
        Function to use to convert names into statistics

    Returns
    -------
    str
    '''
    preferred_name = selector(monosaccharide_synonyms.get(name, [name]) + [name], key=key)
    return preferred_name


def is_a(node, target, tolerance=0, include_modifications=True, include_substituents=True, exact=True,
         short_circuit=False, ignore_ring=True, **kwargs):
    '''
    Perform a semi-fuzzy match between `node` and `target` where node is the unqualified
    residue queried and target is the known residue to be matched against.

    Forwards all unmatched arguments to :func:`~.monosaccharide_similarity`

    Parameters
    ----------
    node: Monosaccharide or Substituent
        Object to be identified
    target: Monosaccharide, Substituent or str
        The reference type. May be a |str| object which is used to look up a |Monosaccharide| by name in
        :obj:`glypy.monosaccharides`
    tolerance: int
        The error tolerance for the search
    include_modifications: bool
        Whether or not to include modifications in comparison. Defaults to |True|
    include_substituents: bool
        Whether or not to include substituents in comparison. Defaults to |True|
    exact: bool
        Whether or not to penalize for unmatched attachments. Defaults to |True|
    Returns
    -------
    bool

    '''
    res = 0
    qs = 0
    if isinstance(target, basestring):
        target = monosaccharides[target]

    if isinstance(node, Substituent):
        if not isinstance(target, Substituent):
            return False
        else:
            res += node.name == target.name
            qs += 1
    else:
        if not isinstance(target, Monosaccharide):
            return False
        res, qs = monosaccharide_similarity(node, target, include_modifications=include_modifications,
                                            include_substituents=include_substituents,
                                            include_children=False, exact=exact, ignore_ring=ignore_ring,
                                            short_circuit_after=tolerance if short_circuit else None, **kwargs)
    threshold = (qs - res) <= tolerance
    return threshold


def identify(node, blacklist=None, tolerance=0, include_modifications=True, include_substituents=True,
             ignore_ring=True, **kwargs):
    '''
    Attempt to find a common usage name for the given |Monosaccharide|, `node`. The name is determined by
    performing an incremental comparison of the traits of `node` with each named residue in the database
    accessed at :obj:`glypy.monosaccharides`.

    Forwards all unmatched arguments to :func:`~.monosaccharide_similarity`

    Parameters
    ----------
    node: Monosaccharide
        Object to be identified
    blacklist: list
        The set of all monosaccharides to not attempt matching against, because they are too general.
    tolerance: int
        The error tolerance for the search
    include_modifications: bool
        Whether or not to include modifications in comparison. Defaults to |True|
    include_substituents: bool
        Whether or not to include substituents in comparison. Defaults to |True|

    Returns
    -------
    str

    Raises
    ------
    IdentifyException:
        When a suitable name cannot be found.

    See Also
    --------
    is_a
    preferred_name
    monosaccharide_similarity
    '''
    if blacklist is None:
        blacklist = {"Pen", "Hex", "Hep", "Oct", "Non"}
    for name, structure in monosaccharides_ordered:
        if name in blacklist:
            continue
        if is_a(node, structure, tolerance, include_modifications, include_substituents, ignore_ring=ignore_ring,
                **kwargs):
            return get_preferred_name(name)
    raise IdentifyException("Could not identify {}".format(node))


class IdentifyException(KeyError):
    pass


def naive_name_monosaccharide(monosaccharide):
    '''
    Generate a generic name for `monosaccharide`, based loosely on IUPAC
    naming schema without including information about linkage.

    The tendency for monosaccharides of superclass > 7 to have special names,
    which will be used preferentially if possible.

    Parameters
    ----------
    monosaccharide: Monosaccharide

    Returns
    -------
    str:
        A simple name based on `SuperClass`, modifications, and substituents.

    See Also
    --------
    :func:`glypy.io.nomenclature.identity.identify`

    '''
    try:
        c = monosaccharide.clone()
        if not isinstance(c, Monosaccharide):
            return None
        strip_derivatization(c)
        try:
            if monosaccharide.superclass.value > 6:
                return identify(c, tolerance=0)
        except Exception:
            pass
        c.anomer = None
        return identify(c)
    except IdentifyException:
        try:
            c.stem = None
            c.configuration = None
            return identify(c)
        except IdentifyException:
            return "".join(mod.name for mod in list(c.modifications.values()) if mod.name != 'aldi') +\
                   c.superclass.name.title() + ''.join([''.join(map(str.title, subst.name.split("_")))[:3]
                                                        for p, subst in c.substituents()])


def split_along_axis(monosaccharides, axis):
    getter = operator.attrgetter(axis)
    groups = groupby(monosaccharides, getter)
    return groups


def residue_list_to_tree(monosaccharides, axes=('anomer', 'superclass', 'stem', 'configuration')):
    root = split_along_axis(monosaccharides, axes[0])
    if len(axes) > 1:
        for level, group in list(root.items()):
            root[level] = residue_list_to_tree(group, axes[1:])
    return root


class MonosaccharideIdentifier(object):
    def __init__(self, reference_index=None, **kwargs):
        if reference_index is None:
            reference_index = dict(named_structures.monosaccharides)
        self.reference_index = dict(reference_index)
        self.trait_tree = residue_list_to_tree(set(self.reference_index.values()))
        self.name_map = self._build_name_map()

    def _build_name_map(self):
        by_monosaccharide = groupby(self.reference_index.items(), lambda x: x[1])
        monosaccharide_to_name = {
            k: min([vi[0] for vi in v], key=len)
            for k, v in by_monosaccharide.items()
        }
        return monosaccharide_to_name

    def _find_potential_matches(self, monosaccharide, exact_candidates=False, **kwargs):
        anomer = monosaccharide.anomer
        candidates = []
        members = self.trait_tree[anomer]
        candidates.append(members)
        if anomer != Anomer.x:
            candidates.append(self.trait_tree[Anomer.x])
        superclass = monosaccharide.superclass
        next_candidates = []
        for candidate in candidates:
            if not candidate:
                continue
            next_candidates.append(candidate[superclass])
            if superclass != SuperClass.x:
                next_candidates.append(candidate[SuperClass.x])
        candidates = next_candidates
        next_candidates = []
        stem = monosaccharide.stem
        for candidate in candidates:
            if not candidate:
                continue
            next_candidates.append(candidate[stem])
            if stem != (Stem.x, ):
                next_candidates.append(candidate[(Stem.x, )])
        candidates = next_candidates
        next_candidates = []
        configuration = monosaccharide.configuration
        for candidate in candidates:
            if not candidate:
                continue
            next_candidates.append(candidate[configuration])
            if configuration != (Configuration.x, ):
                next_candidates.append(candidate[(Configuration.x, )])
        candidates = []
        for c in next_candidates:
            candidates.extend(c)
        is_a_potential = {}
        kwargs.setdefault('exact', True)
        kwargs.setdefault('treat_null_as_wild', False)
        kwargs.setdefault('match_attachement_positions', True)
        for c in candidates:
            if is_a(monosaccharide, c, exact=exact_candidates):
                a, b = monosaccharide_similarity(
                    monosaccharide, c, **kwargs)
                is_a_potential[c] = a / float(b)
        return is_a_potential

    def query(self, monosaccharide, **kwargs):
        is_a_potential = self._find_potential_matches(monosaccharide, **kwargs)
        if not is_a_potential:
            return None
        match = max(is_a_potential.items(), key=lambda x: x[1])[0]
        return match

    def identify(self, monosaccharide, **kwargs):
        template = self.query(monosaccharide, **kwargs)
        if template is not None:
            return self.name_map[template]
        else:
            raise IdentifyException(monosaccharide)
