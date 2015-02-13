import functools
from ...structure import named_structures


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
    for pos, mod in target.modifications.items():
        res += (mod in node.modifications.values())
        qs += 1
    node_subs = list(node for p, node in node.substituents())
    for pos, sub in target.substituents():
        res += (sub in node_subs)
    return (res - qs) >= tolerance


class WrappedQuery(object):
    '''
    A little python sorcery to conveniently dynamically generate identity methods
    at run-time. Once a predicate is made, it is cached to save time. This object
    is substituted for the module in `sys.modules`, allowing us to override normal behavior
    like `:func:getattr`. Also can fall back to the real module directly by accessing :attr:`module`
    '''
    def __init__(self, module):
        self.module = module

    def __getattr__(self, name):
        try:
            return object.__getattr__(self, name)
        except AttributeError:
            try:
                partial_fn = functools.partial(is_a, target=named_structures.monosaccharides[name.replace('is_', "")])
                setattr(self, name, partial_fn)
                return partial_fn
            except:
                return getattr(self.module, name)


if __name__ != '__main__':
    import sys
    self = sys.modules[__name__]
    sys.modules[__name__] = WrappedQuery(self)
