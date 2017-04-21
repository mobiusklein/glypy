import logging
import json
try:
    from itertools import izip_longest
except:
    from itertools import zip_longest as izip_longest
from collections import defaultdict

logger = logging.getLogger(__name__)


def _str_dump_multimap(mm):  # pragma: no cover
    d = defaultdict(list)
    for k, v in mm.items():
        d[k].append(repr(v))
    return json.dumps(d, sort_keys=True, indent=4)


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

    def popv(self, value):
        for k in self:
            if value in self[k]:
                self.pop(k, value)
                return len(self[k])
        return None

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

    def lists(self):
        return self.contents.items()

    def __len__(self):
        '''
        Returns the number of items in :attr:`contents`
        '''
        return sum(len(self[k]) for k in self)

    def __repr__(self):  # pragma: no cover
        return _str_dump_multimap(self)

    def __eq__(self, other):
        if other is None:  # pragma: no cover
            return False
        for a in self.keys():
            if a in other:
                if self[a] != other[a]:
                    return False
            else:
                if self[a] != []:
                    return False
        for b in other.keys():
            if b in self:
                continue
            else:
                if other[b] != []:
                    return False
        # for a, b in izip_longest(self.items(), other.items()):
        #     if a != b:
        #         return False
        return True

    def __ne__(self, other):
        return not self == other

    def update(self, mapping):
        for k, v in mapping.items():
            self[k] = v

    def has_value(self, value):
        for v in self.values():
            if v == value:
                return True
        return False


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

    def values(self):
        for key in self.key_order:
            for v in self[key]:
                yield v

    def items(self):
        '''
        As in :class:`MultiMap`, but items are yielded in the order their keys were
        first added.
        '''
        for key in self.key_order:
            for v in self[key]:
                yield (key, v)

    def lists(self):
        for key in self.key_order:
            yield key, self[key]

    def reorder(self, new_order):
        assert set(new_order) == set(self.key_order)
        self.key_order = list(new_order)

    def __setitem__(self, key, value):
        if key not in self.key_order:
            self.key_order.append(key)
        self.contents[key].append(value)

    def __repr__(self):  # pragma: no cover
        return ''.join((repr(self.key_order), '\n', _str_dump_multimap(self)))
