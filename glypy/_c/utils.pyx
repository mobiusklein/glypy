cimport cython

from cpython.object cimport PyObject
from cpython.ref cimport Py_INCREF
from cpython.list cimport PyList_Append
from cpython.dict cimport PyDict_GetItem, PyDict_Size, PyDict_Next, PyDict_SetItem
from cpython.int cimport PyInt_AsLong, PyInt_FromLong


from libc.string cimport strcmp, memcpy, strlen, strncpy, strcpy
from libc.stdlib cimport malloc, free, realloc, atoi, calloc

from glypy._c.compat cimport (
    PyStr_FromStringAndSize, PyStr_AsUTF8AndSize, PyInt_FromString,
    PyStr_InternInPlace)

cdef extern from * nogil:
    int printf (const char *template, ...)
    void qsort (void *base, unsigned short n, unsigned short w, int (*cmp_func)(void*, void*))


# Parser Helpers


cpdef tuple parse_name_count(str string):
    cdef:
        size_t i
        Py_ssize_t n
        Py_ssize_t boundary
        char* cstring
        str elt
        object count

    cstring = PyStr_AsUTF8AndSize(<str>string, &n)
    elt = None
    count = 0
    for i in range(n):
        if cstring[i] == ':' and i < n:
            elt = PyStr_FromStringAndSize(cstring, i)
            count = PyInt_FromString(&cstring[i + 1], NULL, 10)
            return elt, count
    raise ValueError()


cpdef tuple get_parse_tokens(object string):
    cdef:
        size_t i
        Py_ssize_t n
        Py_ssize_t elt_start, elt_end
        list tokens
        str reduced, token
        char* cstring
    if not isinstance(string, str):
        string = str(string)
    reduced = None
    tokens = []
    cstring = PyStr_AsUTF8AndSize(<str>string, &n)
    elt_start = 1
    elt_end = -1
    if n < 2 or (n > 2 and cstring[0] != '{'):
        raise ValueError("Malformed glycan composition!")
    if n == 2:
        if cstring[0] == '{' and cstring[1] == '}':
            return [], reduced
        else:
            raise ValueError("Malformed glycan composition!")
    for i in range(1, n - 1):
        if (cstring[i] == ';' and cstring[i + 1] == ' '):
            elt_end = i
            token = PyStr_FromStringAndSize(
                &cstring[elt_start], elt_end - elt_start)
            PyList_Append(tokens, token)
            elt_start = i + 2
        elif cstring[i] == '}':
            elt_end = i
            token = PyStr_FromStringAndSize(
                &cstring[elt_start], elt_end - elt_start)
            PyList_Append(tokens, token)
            elt_start = i + 1
        elif cstring[i] == '$':
            reduced = PyStr_FromStringAndSize(&cstring[i + 1], n - (i + 1))
            break
    else:
        if cstring[i + 1] == '}':
            elt_end = i + 1
            token = PyStr_FromStringAndSize(
                &cstring[elt_start], elt_end - elt_start)
            PyList_Append(tokens, token)
    return tokens, reduced

cdef struct string_cell:
    char* string
    size_t size


# Serialization Helpers

DEF ARRAY_SIZE = 2 ** 12
cdef string_cell[ARRAY_SIZE] str_ints
cdef dict str_int_overflow_intern_cache = dict()

cdef void init_str_ints():
    cdef int i
    cdef Py_ssize_t z
    cdef str pa
    cdef char* a
    cdef char* atemp

    for i in range(ARRAY_SIZE):
        pa = str(i)
        atemp = PyStr_AsUTF8AndSize(pa, &z)
        a = <char*>malloc(sizeof(char) * z)
        strcpy(a, atemp)
        a[z] = "\0"
        str_ints[i].string = a
        str_ints[i].size = z


init_str_ints()



cdef int get_str_int(int i, string_cell* out):
    '''Get a `string_cell` for a given integer.

    For integers < 2 ** 12, go through the `str_ints` array
    fast path, otherwise use `str_int_overflow_intern_cache`.

    `str_int_overflow_intern_cache` keeps a mapping from PyInt to
    an immortalized PyStr whose internal buffer is used to serve
    the `char*` for `string_cell` instances larger than  2 ** 12 - 1.
    '''
    cdef:
        PyObject* ptemp
        object py_i
        str pa
        Py_ssize_t z
    if i < ARRAY_SIZE:
        out[0] = str_ints[i]
        return 0

    py_i = PyInt_FromLong(i)
    ptemp = PyDict_GetItem(str_int_overflow_intern_cache, py_i)
    if ptemp == NULL:
        pa = str(py_i)
        Py_INCREF(pa)
        PyDict_SetItem(str_int_overflow_intern_cache, py_i, pa)
    else:
        pa = <str>ptemp
    atemp = PyStr_AsUTF8AndSize(pa, &z)
    out.string = atemp
    out.size = z
    return 0


cdef struct monosaccharide_count_cell:
    char* string
    Py_ssize_t size
    int count
    float mass


cdef int cmp_monosaccharide_count_cell(const void* lhs, const void* rhs) nogil:
    cdef monosaccharide_count_cell* left = <monosaccharide_count_cell*>lhs
    cdef monosaccharide_count_cell* right = <monosaccharide_count_cell*>rhs

    if left.mass < right.mass:
        return -1
    if right.mass < left.mass:
        return 1
    if left.string == NULL:
        if right.string == NULL:
            return 0
        else:
            return -1
    if right.string == NULL:
        return 1
    return strcmp(left.string, right.string)


cdef dict intern_table = dict()


cpdef str _prepare_glycan_composition_string(dict composition):
    cdef size_t n = PyDict_Size(composition)
    cdef monosaccharide_count_cell* cells = <monosaccharide_count_cell*>malloc(sizeof(monosaccharide_count_cell) * n)
    cdef:
        Py_ssize_t i
        size_t j, z, k
        PyObject* pkey
        PyObject* pvalue
        PyObject* ptmp
        str name
        int value
        char* text_buffer
        string_cell int_conv

    i = 0
    j = 0
    z = 2
    while PyDict_Next(composition, &i, &pkey, &pvalue):
        value = PyInt_AsLong(<object>pvalue)
        if value == 0:
            cells[j].string = NULL
            cells[j].size = cells[j].count = 0
            cells[j].mass = -1
        else:
            name = str(<object>pkey)
            ptmp = PyDict_GetItem(intern_table, name)
            if ptmp == NULL:
                PyDict_SetItem(intern_table, name, name)
            else:
                name = <str>ptmp

            cells[j].string = PyStr_AsUTF8AndSize(<str>name, &cells[j].size)
            # The separators, of which there are two extra from the last monosaccharide
            z += 3
            z += cells[j].size
            # Number of digits
            z += (value // 10) + 1
            cells[j].count = value
            cells[j].mass = (<object>pkey).mass()
        j += 1
    qsort(<void*>cells, j, sizeof(monosaccharide_count_cell), cmp_monosaccharide_count_cell)
    if z == 2:
        free(cells)
        return "{}"
    text_buffer = <char*>malloc(sizeof(char) * z)

    k = 0
    text_buffer[k] = '{'
    k += 1
    for j in range(n):
        if cells[j].string == NULL:
            continue

        memcpy(&text_buffer[k], cells[j].string, cells[j].size)
        k += cells[j].size
        text_buffer[k] = ':'
        k += 1
        value = cells[j].count
        get_str_int(value, &int_conv)
        memcpy(&text_buffer[k], int_conv.string, int_conv.size)
        k += int_conv.size
        text_buffer[k] = ';'
        k += 1
        text_buffer[k] = ' '
        k += 1
    # We've written the "; " separator for the last monosaccharide,
    # so now we must rewind two steps and write the terminal '}' over
    # them and shrink k by 1 to reflect actual length of the string.
    if k > 1:
        text_buffer[k - 2] = '}'
        k -= 1
    else:
        text_buffer[k] = '}'
        k += 1
    free(cells)
    result = PyStr_FromStringAndSize(text_buffer, k)
    return result