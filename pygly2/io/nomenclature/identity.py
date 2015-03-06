import functools
from ...structure import named_structures, Monosaccharide, Substituent
from ...algorithms.similarity import monosaccharide_similarity
from .synonyms import monosaccharides as monosaccharide_synonyms


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


def is_a(node, target, tolerance=0, include_modifications=True, include_substituents=True, exact=True):
    '''
    Perform a semi-fuzzy match between `node` and `target` where node is the unqualified
    residue queried and target is the known residue to be matched against

    Parameters
    ----------
    node: Monosaccharide or Substituent
        Object to be identified
    target: Monosaccharide, Substituent or str
        The reference type. May be a |str| object which is used to look up a |Monosaccharide| by name in
        :obj:`pygly2.monosaccharides`
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
        target = named_structures.monosaccharides[target]
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
                                            include_children=False, exact=exact)
    return (qs - res) <= tolerance


def identify(node, blacklist=None, tolerance=0, include_modifications=True, include_substituents=True):
    '''
    Attempt to find a common usage name for the given |Monosaccharide|, `node`. The name is determined by
    performing an incremental comparison of the traits of `node` with each named residue in the database
    accessed at :obj:`pygly2.monosaccharides`.

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
        blacklist = ["Hex"]
    for name, structure in named_structures.monosaccharides.items():
        if name in blacklist:
            continue
        if is_a(node, structure, tolerance, include_modifications, include_substituents):
            return get_preferred_name(name)
    raise IdentifyException("Could not identify {}".format(node))


class IdentifyException(Exception):
    pass
