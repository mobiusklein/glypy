import functools
from ...structure import named_structures, Monosaccharide, Substituent
from ...algorithms.similarity import monosaccharide_similarity
from .synonyms import monosaccharides as monosaccharide_synonyms


def get_preferred_name(name, selector=min, key=len):
    preferred_name = selector(monosaccharide_synonyms.get(name, [name]) + [name], key=key)
    return preferred_name


def is_a(node, target, tolerance=0, include_modifications=True, include_substituents=True):
    '''Perform a semi-fuzzy match between `node` and `target` where node is the unqualified
    residue queried and target is the known residue to be matched against'''
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
                                            include_substituents=include_substituents, include_children=False)
    return (qs - res) <= tolerance


def identify(node, blacklist=None, tolerance=0, include_modifications=True, include_substituents=True):
    if blacklist is None:
        blacklist = ["Hex"]
    for name, structure in named_structures.monosaccharides.items():
        if name in blacklist:
            continue
        if is_a(node, structure, tolerance, include_modifications, include_substituents):
            return get_preferred_name(name)
    # if tolerance < 4:
    #     return identify(node, blacklist, tolerance=tolerance + 1)
    raise IdentifyException("Could not identify {}".format(node))


class IdentifyException(Exception):
    pass


class WrappedQuery(object):
    '''
    A little python sorcery to conveniently dynamically generate identity methods
    at run-time. Once a predicate is made, it is cached to save time. This object
    is substituted for the module in `sys.modules`, allowing us to override normal behavior
    like `:func:getattr`. Also can fall back to the real module directly by accessing :attr:`module`
    '''
    def __init__(self, module):
        self._module = module

    def __getattr__(self, name):
        try:
            return object.__getattr__(self, name)
        except AttributeError:
            try:
                partial_fn = functools.partial(is_a, target=named_structures.monosaccharides[name.replace('is_', "")])
                setattr(self, name, partial_fn)
                return partial_fn
            except:
                return getattr(self._module, name)

    def is_a(self, node, target, tolerance=0):
        return is_a(node, target, tolerance)

    def identify(self, node, blacklist=None):
        return identify(node, blacklist)

if __name__ != '__main__':
    import sys
    self = sys.modules[__name__]
    sys.modules[__name__] = WrappedQuery(self)
