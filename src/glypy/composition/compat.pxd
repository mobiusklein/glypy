cdef extern from "compat.h":
    char* PyStr_AsString(str string)
    str PyStr_FromString(char* string)
    long PyInt_AsLong(object i)
    object PyInt_FromLong(long i)
    str PyStr_Format(object format, object args)