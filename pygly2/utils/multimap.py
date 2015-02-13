from collections import defaultdict


class MultiMap(object):
    '''Implements a simple MultiMap data structure on top of a dictionary of lists'''
    def __init__(self, **kwargs):
        self.contents = defaultdict(list)
        for k, v in kwargs.items():
            self.contents[k].append(v)

    def __getitem__(self, key):
        return self.contents[key]

    def __setitem__(self, key, value):
        self.contents[key].append(value)

    def pop(self, key, value):
        '''
        Removes `value` from the collection of values stored at `key`
        and returns the `tuple` `(key, value)`

        Raises
        ------
        IndexError
        KeyError
        '''
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
        '''
        Returns an iterator over the keys of :attr:`contents`
        An alias of :meth:`__iter__`
        '''
        return iter(self)

    def values(self):
        '''
        Returns an iterator over the values of :attr:`contents`
        '''
        for k in self:
            for v in self[k]:
                yield v

    def items(self):
        '''
        Returns an iterator over the items of :attr:`contents`. Each item
        takes the form of `(key, value)`.
        '''
        for k in self:
            for v in self[k]:
                yield (k, v)

    def __len__(self):
        '''
        Returns the number of items in :attr:`contents`
        '''
        return sum(len(self[k]) for k in self)

    def __repr__(self):
        return repr(self.contents)

    def __eq__(self, other):
        for a, b in zip(self.items(), other.items()):
            if a != b:
                return False
        return True

    def __ne__(self, other):
        return self.contents != other.contents


class OrderedMultiMap(MultiMap):
    '''
    Implements a simple MultiMap data structure on top of a dictionary of lists
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
        '''
        Returns an iterator over the keys of :attr:`contents` in the order
        they were added.
        '''
        for key in self.key_order:
            yield key

    #: Alias of :meth:`__iter__`
    keys = __iter__

    def items(self):
        '''
        As in :class:`MultiMap`, but items are yielded in the order their keys were
        first added.
        '''
        for key in self.key_order:
            for v in self[key]:
                yield (key, v)

    def __setitem__(self, key, value):
        if key not in self.key_order:
            self.key_order.append(key)
        self.contents[key].append(value)

    def pop(self, key, value):
        rv = super(OrderedMultiMap, self).pop(key, value)
        if rv == 0:
            self.key_order.pop(self.key_order.index(key))

    def __repr__(self):
        return ''.join((repr(self.key_order), '\n', repr(self.contents)))
