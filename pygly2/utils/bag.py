from collections import defaultdict

class Bag(object):
    '''Implements a simple MultiSet-like data structure on top of a dictionary of lists'''
    def __init__(self, **kwargs):
        self.contents = defaultdict(list)
        for k, v in kwargs.items():
            self.contents[k].append(v)

    def __getitem__(self, key):
        return self.contents[key]

    def __setitem__(self, key, value):
        self.contents[key].append(value)

    def pop(self, key, value):
        objs = self.contents[key]
        objs.pop(objs.index(value))
        if len(objs) == 0:
            self.contents.pop(key)
        return len(objs)

    def __iter__(self):
        return iter(self.contents)

    def __contains__(self, key):
        return key in self.contents

    def keys(self):
        return iter(self)

    def values(self):
        for k in self:
            for v in self[k]:
                yield v

    def items(self):
        for k in self:
            for v in self[k]:
                yield (k, v)

    def __len__(self):
        return sum(len(self[k]) for k in self)

    def __repr__(self):
        return repr(self.contents)

    def __eq__(self, other):
        return self.contents == other.contents

    def __ne__(self, other):
        return self.contents != other.contents

class OrderedBag(Bag):
    '''
    Implements a simple MultiSet-like data structure on top of a dictionary of lists
    that remembers the order keys were first inserted in.
    '''
    def __init__(self, **kwargs):
        self.contents = defaultdict(list)
        self.key_order = []
        for k, v in kwargs.items():
            if k not in self.key_order:
                self.key_order.append(k)
            self.contents[k].append(v)

    def __iter__(self):
        for key in self.key_order:
            yield key

    def items(self):
        for key in self.key_order:
            for v in self[key]:
                yield (key, v)

    def __setitem__(self, key, value):
        if key not in self.key_order:
            self.key_order.append(key)
        self.contents[key].append(value)

    def pop(self, key, value):
        rv = super(OrderedBag, self).pop(key, value)
        if rv == 0:
            self.key_order.pop(self.key_order.index(key))

    def __repr__(self):
        return ''.join((repr(self.key_order), '\n', repr(self.contents)))


