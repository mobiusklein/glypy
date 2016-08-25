import warnings
import os
import random
import sys
import gzip

from collections import defaultdict
try:  # pragma: no cover
    import cPickle as pickle
except:  # pragma: no cover
    import pickle
try:  # pragma: no cover
    from lxml import etree as ET
except ImportError:  # pragma: no cover
    try:
        from xml.etree import cElementTree as ET
    except:
        from xml.etree import ElementTree as ET

try:   # pragma: no cover
    from cStringIO import StringIO
except:  # pragma: no cover
    try:
        from StringIO import StringIO
    except:
        from io import StringIO

try:  # pragma: no cover
    i128 = long
    basestring = basestring
except:  # pragma: no cover
    i128 = int
    basestring = (bytes, str)


def cyclewarning():
    '''
    Used to warn users about the presence of cylcical glycans, which are harder to
    reason about for crossring_cleavages.
    '''
    warnings.warn("A cyclic glycan may be present. They may not cross-ring fragment correctly.", stacklevel=3)


def opener(obj, mode='r'):
    '''
    Try to use `obj` to access a file-object. If `obj` is a string, assume
    it denotes a path to a file, and open that file in the specified mode.
    If `obj` has an attribute `read`, assume it
    itself is a file-like object and return it.

    Parameters
    ----------
    obj: basestring or file-like object
        If `obj` is a base string it is treated like a file path, else if it supports
        the file-like operation `read`, return the object unchanged.
    mode: str, optional
        The mode, if any, to open `obj` with if it is a file path. Defaults to 'r', `read`

    '''
    if isinstance(obj, basestring):
        if obj[-2:] == 'gz':  # pragma: no cover
            return gzip.open(obj, mode)
        return open(obj, mode)
    elif hasattr(obj, "read"):
        return obj
    else:  # pragma: no cover
        raise IOError("Can't find a way to open {}".format(obj))


def invert_dict(d):
    return {v: k for k, v in d.items()}


def make_counter(start=1):
    '''
    Create a functor whose only internal piece of data is a mutable container
    with a reference to an integer, `start`. When the functor is called, it returns
    current `int` value of `start` and increments the mutable value by one.

    Parameters
    ----------
    start: int, optional
        The number to start counting from. Defaults to `1`.

    Returns
    -------
    int:
        The next number in the count progression.
    '''
    start = [start]

    def count_up():
        ret_val = start[0]
        start[0] += 1
        return ret_val
    return count_up


def identity(x):   # pragma: no cover
    return x


def nullop(*args, **kwargs):   # pragma: no cover
    pass


def chrinc(a='a', i=1):
    return chr(ord(a) + i)


def make_struct(name, fields, debug=False):
    '''
    A convenience function for defining plain-old-data (POD) objects that are optimized
    for named accessor lookup, unlike `namedtuple`. If the named container does not
    require any special logic and won't be extended, the resulting structure is best for
    storing and accessing the data.

    Parameters
    ----------
    name: str
        The name of the new class structure
    fields: iterable of str

    '''
    template = ('''
class {name}(object):
    __slots__ = {fields!r}
    def __init__(self, {args}):
        {self_fields} = {args}
    def __getitem__(self, idx):
        return getattr(self, fields[idx])
    def __getstate__(self):
        return ({self_fields})
    def __setstate__(self, state):
        {self_fields} = state
    def __repr__(self):
        rep = "<{name}"
        for f in {fields!r}:
            rep += " " + f + "=" + str(getattr(self, f))
        rep += ">"
        return rep
    def __eq__(self, other):
        for f in {fields!r}:
            if getattr(self, f) != getattr(other, f):
                return False
        return True
    def __ne__(self, other):
        return not self == other

    @property
    def __dict__(self):
        d = dict()
        for f in {fields!r}:
            d[f] = getattr(self, f)
        return d
    ''').format(
        name=name,
        fields=fields,
        args=','.join(fields),
        self_fields=','.join('self.' + f for f in fields))
    d = {'fields': fields}
    if debug:  # pragma: no cover
        print(template)
    exec(template, d)
    result = d[name]
    # Patch the class to support pickling, as is done for namedtuple
    try:
        result.__module__ = sys._getframe(1).f_globals.get('__name__', '__main__')
    except (AttributeError, ValueError):  # pragma: no cover
        pass
    return result


class ClassPropertyDescriptor(object):
    '''
    Standard Class Property Descriptor Implementation
    '''
    def __init__(self, fget, fset=None):
        self.fget = fget
        self.fset = fset

    def __get__(self, obj, klass=None):
        if klass is None:  # pragma: no cover
            klass = type(obj)
        return self.fget.__get__(obj, klass)()

    def __set__(self, obj, value):
        if not self.fset:  # pragma: no cover
            raise AttributeError("can't set attribute")
        type_ = type(obj)
        return self.fset.__get__(obj, type_)(value)

    def setter(self, func):
        if not isinstance(func, (classmethod, staticmethod)):
            func = classmethod(func)
        self.fset = func
        return self


def classproperty(func):
    '''
    Applies ClassPropertyDescriptor as you would a normal
    @property descriptor
    '''
    if not isinstance(func, (classmethod, staticmethod)):
        func = classmethod(func)

    return ClassPropertyDescriptor(func)


def root(structure):
    """The `root` protocol function.

    Creates a generic method for obtaining the root of a structure representing
    or containing a glycan graph with a single distinct root.

    Parameters
    ----------
    structure : any
        An object that implements the `root` protocol, containing
        a tree structure somewhere inside it.

    Returns
    -------
    Monosaccharide : The root of the |Glycan| tree
    """
    try:
        root = structure.__root__()
        return root
    except AttributeError:
        raise TypeError("{} does not support the `root` protocol".format(structure.__class__.__name__))


def tree(structure):
    """The `tree` protocol function.

    Creates a generic method for obtaining the :class:`Glycan` of a structure representing
    or containing a glycan graph.

    Parameters
    ----------
    structure : any
        An object that implements the `root` protocol, containing
        a tree structure somewhere inside it.

    Returns
    -------
    Glycan
    """
    try:
        tree = structure.__tree__()
    except AttributeError:  # pragma: no cover
        raise TypeError("{} does not support the `tree` protocol.".format(structure.__class__.__name__))
    return tree


def groupby(ungrouped_list, key_fn=identity):
    groups = defaultdict(list)
    for item in ungrouped_list:
        key_value = key_fn(item)
        groups[key_value].append(item)
    return groups


def where(iterable, fn):
    return [i for i, k in enumerate(iterable) if fn(k)]


def uid(n=128):
    int_ = random.getrandbits(n)
    return int_
