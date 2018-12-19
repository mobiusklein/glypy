cdef class EnumValue(object):
    cdef:
        public object group
        public str name
        public object value
        public set names
        public object _hash
