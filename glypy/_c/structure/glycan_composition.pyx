from cpython.object cimport PyObject
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem, PyDict_Next
from cpython.list cimport PyList_Append, PyList_GET_ITEM, PyList_GET_SIZE, PyList_Sort
from cpython.tuple cimport PyTuple_GET_ITEM
from cpython.int cimport PyInt_AsLong, PyInt_FromLong
from glypy.composition.ccomposition cimport CComposition
from glypy.composition import formula

from glypy._c.utils cimport _prepare_glycan_composition_string

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

    cdef void _add_from(self, _CompositionBase other):
        cdef:
            Py_ssize_t i
            PyObject* pkey
            PyObject* pval
            PyObject* ptmp
            long value

        i = 0
        while PyDict_Next(other, &i, &pkey, &pval):
            ptmp = PyDict_GetItem(self, <object>pkey)
            if ptmp == NULL:
                value = 0
            else:
                value = PyInt_AsLong(<object>ptmp)
            value += PyInt_AsLong(<object>pval)
            PyDict_SetItem(self, <object>pkey, PyInt_FromLong(value))

    cdef void _subtract_from(self, _CompositionBase other):
        cdef:
            Py_ssize_t i
            PyObject* pkey
            PyObject* pval
            PyObject* ptmp
            long value

        i = 0
        while PyDict_Next(other, &i, &pkey, &pval):
            ptmp = PyDict_GetItem(self, <object>pkey)
            if ptmp == NULL:
                value = 0
            else:
                value = PyInt_AsLong(<object>ptmp)
            value -= PyInt_AsLong(<object>pval)
            PyDict_SetItem(self, <object>pkey, PyInt_FromLong(value))

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

    cpdef str serialize(self):
        # Causes segmentation fault in unrelated mass calculation code?
        # cdef str result = _prepare_glycan_composition_string(self)
        cdef str result = _reformat(self)
        reduced = self._reducing_end
        if reduced is not None:
            result = "%s$%s" % (result, formula(reduced.total_composition()))
        return result


cdef str _reformat(dict self):
    cdef:
        list tokens, strings
        Py_ssize_t i, n
        PyObject* k
        PyObject* v
        object key
        int count
        tuple token

    tokens = []
    strings = []
    i = 0
    while PyDict_Next(self, &i, &k, &v):
        count = PyInt_AsLong(<object>v)
        if count != 0:
            key = <object>k
            PyList_Append(tokens, (key.mass(), str(key), (<object>v)))
    PyList_Sort(tokens)
    i = 0
    n = PyList_GET_SIZE(tokens)
    for i in range(n):
        token = <tuple>PyList_GET_ITEM(tokens, i)
        name = <object>PyTuple_GET_ITEM(token, 1)
        cnt = <object>PyTuple_GET_ITEM(token, 2)
        part = "%s:%d" % (name, cnt)
        PyList_Append(strings, part)
    return "{%s}" % '; '.join(strings)