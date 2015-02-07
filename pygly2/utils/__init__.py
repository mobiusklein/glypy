__all__ = ["enum", 'opener', 'make_counter', 'invert_dict', 'identity', 'nullop']


try:
    from lxml import etree as ET
except ImportError:
    try:
        from xml.etree import cElementTree as ET
    except:
        from xml.etree import ElementTree as ET

try:
    from cStringIO import StringIO
except:
    try:
        from StringIO import StringIO
    except:
        from io import StringIO


import enum


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
        return open(obj, mode)
    elif hasattr(obj, "read"):
        return obj
    else:
        raise IOError("Can't find a way to open {}".format(obj))


def invert_dict(d):
    return {v: k for k, v in d.items()}


def make_counter(start=1, label=""):
    '''
    Create a functor whose only internal piece of data is a mutable container
    with a reference to an integer, `start`. When the functor is called, it returns
    current `int` value of `start` and increments the mutable value by one. 

    Parameters
    ----------
    start: int, optional
        The number to start counting from. Defaults to `1`.
    return: int
        The next number in the count progression.
    '''
    start = [start]
    def count_up():
        ret_val = start[0]
        start[0] += 1
        return ret_val
    return count_up


def identity(x):
    return x


def nullop(*args, **kwargs):
    pass