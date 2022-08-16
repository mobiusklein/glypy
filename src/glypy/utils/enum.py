from functools import total_ordering

from six import add_metaclass
try:
    intern
except NameError:
    from sys import intern # pylint: disable=no-name-in-module


class EnumValue(object):
    '''Represents a wrapper around an value with a name to identify it and
    more rich comparison logic. A value of an enumerated type'''

    __slots__ = ('group', 'name', 'value', 'names', '_hash')

    def __init__(self, group, name, value, other_names=None):
        self.name = intern(name)
        self.value = value
        self.names = {name} | (other_names or set())
        self.group = group
        self._hash = None

    def __hash__(self):
        if self._hash is None:
            self._hash = hash(self.name)
        return self._hash

    def __index__(self):  # pragma: no cover
        return self.value

    def __int__(self):
        return int(self.value)

    def __eq__(self, other):
        try:
            if self.group is not other.group:
                return False
            # if self is other:
            #     return True
            # return self.value == other.value or self.names == other.names
            return self is other
        except AttributeError:
            return self.value == other or other in self.names

    def __and__(self, other):
        return self.value & other

    def __or__(self, other):
        return self.value | other

    def __rand__(self, other):
        return self.value & other

    def __ror__(self, other):
        return self.value | other

    def __xor__(self, other):
        return self.value ^ other

    def __rxor__(self, other):
        return self.value ^ other

    def __ne__(self, other):
        return not self == other

    def __repr__(self):  # pragma: no cover
        return "<{group_name} {name}:{value}>".format(name=self.name,
                                                      group_name=self.group.__name__,
                                                      value=self.value)

    def __str__(self):
        return self.name

    def __reduce__(self):
        return self.group, (self.name,)

    def add_name(self, name, force=False):
        if name not in self.group or force:
            self.names.add(name)
            self.group[name] = self
        else:
            raise KeyError("{} already exists in {}".format(name, self.group))

    def resolve(self, mapping):
        for name in self.names:
            try:
                if name in mapping:
                    return mapping[name]
            except KeyError:  # pragma: no cover
                pass
        raise KeyError("Could not resolve {} against {}".format(self, mapping))

    def __lt__(self, other):
        return self.value < other.value

    def __gt__(self, other):
        return self.value > other.value

    def __le__(self, other):
        return self.value <= other.value

    def __ge__(self, other):
        return self.value >= other.value

    def int_value(self):
        return (int(self.value))

try:
    _EnumValue = EnumValue
    _has_c = True
    from glypy.utils.cenum import EnumValue
except ImportError:
    _has_c = False


debug = False
QMARK = "?"


class EnumMeta(type):
    '''
    A metaclass for types hosting enumerated members. Class attributes are
    automatically translated into ``EnumValue`` objects with the specified value
    and a name string matching the specified attribute name. Newly assigned attributes
    are wrapped dynamically.

    The class itself can be treated like a dictionary for looking up these attributes
    by name or value, and these values can be iterated over.

    ..note::
        Why is this necessary? It's probably not. I chose to write it initially for
        compatibility with a previous iteration of the library, but later decided it
        was worth keeping for two reasons:

        1. Avoids ``stringly-typing`` the library. Comparison of strings is a slippery slope
           to raw magic strings littering the codebase.
        2. Richer comparison behavior allows these same names to be used in different modes with the
           same symbol.
        3. Namespacing of EnumValue objects makes it easier to avoid accidental name collisions
           when comparing EnumValues instead of magic strings.

    '''

    def __new__(cls, name, parents, attrs):
        if attrs.get('__doc__') is None:
            attrs['__doc__'] = "EnumType"
        enum_type = type.__new__(cls, name, parents, attrs)
        mapped = {}
        attr_pairs = list(attrs.items())
        EnumType = attrs.get("__enum_type__", EnumValue)
        for label, value in attr_pairs:
            if not label.startswith("__") or label == "mro":
                attrs.pop(label)
                delattr(enum_type, label)
                enum_value = EnumType(enum_type, label, value)
                if value in mapped:
                    try:
                        mapped[value].add_name(label)
                        setattr(enum_type, label, mapped[value])
                    except KeyError as e:
                        print(e)
                else:
                    mapped[value] = enum_value
                    setattr(enum_type, label, enum_value)
        enum_type[QMARK] = None
        for tp in parents:
            if isinstance(tp, EnumMeta):
                for label, value in tp:
                    if label == QMARK:
                        continue
                    super(EnumMeta, enum_type).__setattr__(label, value)
        return enum_type

    def __iter__(self):
        for attr, val in self.__dict__.items():
            if not attr.startswith("__") or attr == "mro":
                yield (attr, val)

    def __contains__(self, k):
        val = (k in self.__dict__) or (k in self.__dict__.values())
        return val

    def __getitem__(self, k):
        return self.translate(k)  # pylint: disable=no-value-for-parameter

    def __setattr__(self, k, v):
        """Intercept attribute assignment, wrapping values in
        :class:`EnumValue`

        Parameters
        ----------
        k : str
            Name to be set
        v : object
            The value to assign. If it is not of type :class:`EnumValue` it
            will be wrapped like `EnumValue(self, name=k, value=v)`
        """
        if isinstance(v, EnumValue):
            v.names.add(k)
            super(EnumMeta, self).__setattr__(k, v)
        else:
            name = self.name(v)  # pylint: disable=no-value-for-parameter
            if name is not None:
                self[name].add_name(k)  # pylint: disable=unsubscriptable-object
            else:
                super(EnumMeta, self).__setattr__(k, EnumValue(self, k, v))

    def __setitem__(self, k, v):
        setattr(self, k, v)

    def translate(self, k):
        '''
        Attempt to translate the input object ``k`` into a data member of the Enum.

        First try to find an element of ``self`` by hashing it against member names.

        Then try to find an element of ``self`` by searching for a member in self that
        is value-equal to ``k``

        Otherwise throw a :exc:`KeyError`

        Parameters
        ----------
        k: object
            The value to be translated.

        Returns
        -------
        :class:`EnumValue`

        '''
        # global debug
        # if debug:
        #     print "Translating", k

        if k in self.__dict__:
            return self.__dict__[k]
        elif k in self.__dict__.values():
            return self[self.name(k)]  # pylint: disable=unsubscriptable-object,no-value-for-parameter
        else:
            raise KeyError("Could not translate {0} through {1}".format(k, self))

    def name(self, v):
        for k, val in reversed(list(self)):
            if v == val:
                return k

    def __repr__(self):
        return "<Enum {0}>".format(self.__name__)

    __call__ = translate


try:
    _EnumMeta = EnumMeta
    _has_c = True
    from glypy.utils.cenum import EnumMeta
except ImportError:
    _has_c = False


@add_metaclass(EnumMeta)
class Enum(object):
    '''
    A simple class implementing :class:`EnumMeta`. Useful base type for other
    enumerated types.
    '''
    # __metaclass__ = EnumMeta

    def __init__(self):
        raise Exception("This class is not meant to be instantiated. Reference its attribute members directly")
