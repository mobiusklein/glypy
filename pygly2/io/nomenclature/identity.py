import functools
from ...structure import named_structures
from .synonyms import monosaccharides as monosaccharide_synonyms


def get_preferred_name(name, selector=min, key=len):
    preferred_name = selector(monosaccharide_synonyms[name], key=key)
    return preferred_name


def is_a(node, target, tolerance=0):
    '''Perform a semi-fuzzy match between `node` and `target` where node'''
    if isinstance(target, basestring):
        target = named_structures.monosaccharides[target]
    res = 0
    qs = 0
    res += (node.superclass == target.superclass) or (target.superclass.value is None)
    qs += 1
    res += (node.stem == target.stem) or (target.stem[0].value is None)
    qs += 1
    res += (node.configuration == target.configuration) or (target.configuration[0].value is None)
    qs += 1
    node_mods = list(node.modifications.values())
    for pos, mod in target.modifications.items():
        check = (mod in node_mods)
        res += check
        if check:
            node_mods.pop(node_mods.index(mod))
        qs += 1
    qs += len(node_mods)
    node_subs = list(node for p, node in node.substituents())
    for pos, sub in target.substituents():
        check = (sub in node_subs)
        res += check
        if check:
            node_subs.pop(node_subs.index(sub))
        qs += 1
    qs += len(node_subs)
    return (res - qs) >= tolerance


def identify(node, blacklist=None):
    if blacklist is None:
        blacklist = ["Hex", "Hexose"]
    for name, structure in named_structures.monosaccharides.items():
        if name in blacklist:
            continue
        if is_a(node, structure):
            return get_preferred_name(name)
    for name in blacklist:
        if is_a(node, named_structures.monosaccharides[name]):
            return get_preferred_name(name)
    raise IdentifyException("Could not identify {}" % node)


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
