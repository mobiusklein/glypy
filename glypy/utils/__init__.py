import pkg_resources

from .base import (opener, make_counter, invert_dict, identity,
                   nullop, chrinc, make_struct, classproperty, cyclewarning,
                   root, tree, groupby, pickle, ET, StringIO)

__all__ = ['opener', 'make_counter', 'invert_dict', 'identity', 'nullop',
           "chrinc", "make_struct", "classproperty", "cyclewarning",
           "root", "tree", "groupby"]

pkg_resources.declare_namespace('glypy.utils')
