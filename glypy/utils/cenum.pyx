# from glypy.composition.compat cimport PyStr_InternInPlace

cdef class EnumValue(object):
    '''Represents a wrapper around an value with a name to identify it and
    more rich comparison logic. A value of an enumerated type'''


    def __init__(self, group, name, value, other_names=None):
        # PyStr_InternInPlace(name)
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
            self = <EnumValue>other
            other = temp
        else:
            self_t = <EnumValue>self
        return self_t.value | other


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
