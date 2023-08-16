from glypy.composition.ccomposition cimport CComposition

cdef class _CompositionBase(dict):
    cdef:
        public object _mass
        public object _reducing_end
        public CComposition _composition_offset

    cpdef object _getitem_fast(self, object key)
    cpdef object _setitem_fast(self, object key, object value)

    cpdef object _update_from_typed_map(self, _CompositionBase template, bint copy_nodes=*)
    cpdef str serialize(self)

    cdef void _add_from(self, _CompositionBase other)
    cdef void _subtract_from(self, _CompositionBase other)

    cpdef set_composition_offset(self, CComposition composition)
    cpdef CComposition get_composition_offset(self)
    cpdef _invalidate(self)
