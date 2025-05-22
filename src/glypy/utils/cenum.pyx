cimport cython
from cpython cimport PyErr_SetString, PyErr_SetObject
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem
from cpython.object cimport PyObject
from glypy._c.compat cimport PyInt_AsLong


cdef str QMARK = "?"


cdef class EnumMeta(type):
    def __cinit__(self, name, parents, attrs):
        self.name_map = {}
        self.value_map = {}
        EnumType = attrs.get("__enum_type__", EnumValue)
        mapped = dict()
        for label, value in list(attrs.items()):
            if not (label.startswith("__") or label == "mro"):
                attrs.pop(label)
                delattr(self, label)
                enum_value = EnumType(self, label, value)
                if value in mapped:
                    mapped[value].add_name(label)
                    self.put(label, mapped[value])
                else:
                    mapped[value] = enum_value
                    self.put(label, enum_value)

        self.put(QMARK, None)
        for tp in parents:
            if isinstance(tp, EnumMeta):
                for label, value in tp:
                    if label == QMARK:
                        continue
                    self.put(label, value)

    def __dir__(self):
        return list(self.name_map.keys())

    def __iter__(self):
        for attr, val in self.name_map.items():
            if not attr.startswith("__") or attr == "mro":
                yield (attr, val)

    def __contains__(self, k):
        val = (k in self.name_map) or (k in self.value_map)
        return val

    def __getattr__(self, name):
        try:
            return self.translate(name)
        except KeyError:
            raise AttributeError(name)

    def __setattr__(self, name, value):
        self.put(name, value)

    def __getitem__(self, name):
        return self.translate(name)

    def __setitem__(self, name, value):
        self.put(name, value)

    cdef int put(self, name, value) except -1:
        if not isinstance(value, EnumValue):
            existing_value = self.get(value)
            if existing_value is not None:
                value = existing_value
                try:
                    value.add_name(name)
                except KeyError as e:
                    PyErr_SetObject(KeyError, str(e))
                    return -1
            else:
                value = EnumValue(self, name, value)
        PyDict_SetItem(self.name_map, name, value)
        PyDict_SetItem(self.value_map, value.value, name)
        return 0

    cdef EnumValue get(self, k):
        cdef:
            PyObject* presult
        presult = PyDict_GetItem(self.name_map, k)
        if presult == NULL:
            presult = PyDict_GetItem(self.value_map, k)
            if presult == NULL:
                return None
            else:
                presult = PyDict_GetItem(self.name_map, <object>presult)
                if presult == NULL:
                    return None
        return <EnumValue>presult

    cpdef translate(self, k):
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
        cdef:
            EnumValue result
        result = self.get(k)
        if result is None:
            raise KeyError("Could not translate {0} through {1}".format(k, self))
        return result

    cpdef EnumValue name(self, v):
        cdef:
            PyObject* presult
            object result
        presult = PyDict_GetItem(self.value_map, v)
        if presult == NULL:
            return None
        else:
            result = <object>presult
            return result

    def __repr__(self):
        return "<Enum {0}>".format(self.__name__)

    def __call__(self, v):
        return self.translate(v)


@cython.c_api_binop_methods(True)
cdef class EnumValue(object):
    '''Represents a wrapper around an value with a name to identify it and
    more rich comparison logic. A value of an enumerated type'''

    def __init__(self, group, name, value, other_names=None):
        self.name = name
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
        cdef:
            EnumValue other_ev
            object temp
        if not isinstance(self, EnumValue):
            temp = self
            self = <EnumValue>other
            other = temp
        if isinstance(other, EnumValue):
            other_ev = <EnumValue>other
            if self.group is not other_ev.group:
                return False
            return self is other_ev
        else:
            return self.value == other or other in self.names

    def __and__(self, other):
        cdef:
            object temp
            EnumValue self_t
        if not isinstance(self, EnumValue):
            temp = self
            self_t = <EnumValue>other
            other = temp
        else:
            self_t = <EnumValue>self
        return self_t.value & other

    def __or__(self, other):
        cdef:
            object temp
            EnumValue self_t
        if not isinstance(self, EnumValue):
            temp = self
            self_t = <EnumValue>other
            other = temp
        else:
            self_t = <EnumValue>self
        return self_t.value | other

    def __xor__(self, other):
        cdef:
            object temp
            EnumValue self_t
        if not isinstance(self, EnumValue):
            temp = self
            self_t = <EnumValue>other
            other = temp
        else:
            self_t = <EnumValue>self
        return self_t.value ^ other

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

    cpdef add_name(self, basestring name, bint force=False):
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

    cpdef int int_value(self) except *:
        return PyInt_AsLong(int(self.value))



cdef class IntEnumValue(EnumValue):
    def __init__(self, group, name, value, other_names=None):
        super(IntEnumValue, self).__init__(group, name, value, other_names)
        self._int_value = PyInt_AsLong(self.value)

    cpdef int int_value(self) except *:
        return self._int_value