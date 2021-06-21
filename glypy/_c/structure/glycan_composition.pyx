from cpython.object cimport PyObject
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem, PyDict_Next

from glypy.composition.ccomposition cimport CComposition

cdef CComposition water = CComposition("H2O")


cdef object ZERO = 0


cdef class _CompositionBase(dict):

    @classmethod
    def _empty(cls):
        cdef _CompositionBase inst = cls.__new__(cls)
        inst._mass = None
        inst._reducing_end = None
        inst._composition_offset = CComposition._create(water)
        return inst

    cpdef object _getitem_fast(self, object key):
        cdef:
            PyObject* ptmp

        ptmp = PyDict_GetItem(self, key)
        if ptmp == NULL:
            return ZERO
        return <object>ptmp

    cpdef object _setitem_fast(self, object key, object value):
        PyDict_SetItem(self, key, value)

    def __reduce__(self):
        return self.__class__, (), self.__getstate__()

    def __getstate__(self):
        d = {
            'mapping': dict(self),
            'reducing_end': self._reducing_end,
            'composition_offset': self._composition_offset
        }
        return d

    def __setstate__(self, state):
        self.update(state['mapping'])
        self._reducing_end = state['reducing_end']
        self._composition_offset = state['composition_offset']

    cpdef _update_from_typed_map(self, _CompositionBase template, bint copy_nodes=False):
        cdef:
            PyObject* pkey
            PyObject* pval
            Py_ssize_t pos
        pos = 0
        if copy_nodes:
            while PyDict_Next(<dict>template, &pos, &pkey, &pval):
                self._setitem_fast((<object>pkey).clone(), <object>pval)
        else:
            while PyDict_Next(<dict>template, &pos, &pkey, &pval):
                self._setitem_fast(<object>pkey, <object>pval)
        reduced = template._reducing_end
        if reduced is not None:
            self._reducing_end = reduced.clone()
        self._mass = None
