import pkg_resources
try: # pragma: no cover
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
from .base import opener, make_counter, invert_dict, identity, nullop, chrinc, make_struct, classproperty

__all__ = ['opener', 'make_counter', 'invert_dict', 'identity', 'nullop']

pkg_resources.declare_namespace('pygly2.utils')
