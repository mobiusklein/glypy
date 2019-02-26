
cdef class EnumValue(object):
    cdef:
        public object group
        public str name
        public object value
        public set names
        public object _hash

    cpdef add_name(self, basestring name, bint force=*)


cdef class EnumMeta(type):
    cdef:
        public dict name_map
        public dict value_map

    cpdef EnumValue name(self, v)
    cpdef translate(self, k)

    cdef EnumValue get(self, k)
    cdef int put(self, name, value) except -1


cdef str QMARK