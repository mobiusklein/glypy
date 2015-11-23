from collections import defaultdict
import operator
from ...utils import groupby
from ...structure import named_structures, Monosaccharide, Substituent, Anomer, Stem, RingType, SuperClass
from ...algorithms.similarity import (monosaccharide_similarity, has_substituent,
                                      has_modification, has_monosaccharide)
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
                                            include_children=False, exact=exact)
    threshold = (qs - res) <= tolerance
    return threshold


def identify(node, blacklist=None, tolerance=0, include_modifications=True, include_substituents=True):
    '''
    Attempt to find a common usage name for the given |Monosaccharide|, `node`. The name is determined by
    performing an incremental comparison of the traits of `node` with each named residue in the database
    accessed at :obj:`glypy.monosaccharides`.

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
        if is_a(node, structure, tolerance, include_modifications, include_substituents):
            return get_preferred_name(name)
    raise IdentifyException("Could not identify {}".format(node))


class IdentifyException(Exception):
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
        except:
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
