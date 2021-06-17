cimport cython

from cpython.list cimport PyList_Append

from libc.string cimport strcmp, memcpy, strlen, strncpy
from libc.stdlib cimport malloc, free, realloc, atoi, calloc

from glypy._c.compat cimport (
    PyStr_FromStringAndSize, PyStr_AsUTF8AndSize)


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