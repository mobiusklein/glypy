cimport cython

from cpython.object cimport PyObject
from cpython.ref cimport Py_INCREF
from cpython.list cimport PyList_Append
from cpython.dict cimport PyDict_GetItem, PyDict_Size, PyDict_Next, PyDict_SetItem

from libc.string cimport strcmp, memcpy, strlen, strncpy, strcpy
from libc.stdlib cimport malloc, free, realloc, atoi, calloc

from glypy._c.compat cimport (
    PyStr_FromStringAndSize, PyStr_AsUTF8AndSize, PyInt_FromString,
    PyStr_InternInPlace, PyInt_AsLong, PyInt_FromLong)

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

DEF ARRAY_SIZE = 2 ** 8
cdef string_cell[ARRAY_SIZE] str_ints
cdef dict str_int_overflow_intern_cache = dict()

# [[[cog
# import cog
# NULL = '\0'
# for i in range(0, 2**8):
#   cog.outl(fr"""str_ints[{i}] = string_cell("{str(i)}", {len(str(i))})""")
# ]]]
str_ints[0] = string_cell("0", 1)
str_ints[1] = string_cell("1", 1)
str_ints[2] = string_cell("2", 1)
str_ints[3] = string_cell("3", 1)
str_ints[4] = string_cell("4", 1)
str_ints[5] = string_cell("5", 1)
str_ints[6] = string_cell("6", 1)
str_ints[7] = string_cell("7", 1)
str_ints[8] = string_cell("8", 1)
str_ints[9] = string_cell("9", 1)
str_ints[10] = string_cell("10", 2)
str_ints[11] = string_cell("11", 2)
str_ints[12] = string_cell("12", 2)
str_ints[13] = string_cell("13", 2)
str_ints[14] = string_cell("14", 2)
str_ints[15] = string_cell("15", 2)
str_ints[16] = string_cell("16", 2)
str_ints[17] = string_cell("17", 2)
str_ints[18] = string_cell("18", 2)
str_ints[19] = string_cell("19", 2)
str_ints[20] = string_cell("20", 2)
str_ints[21] = string_cell("21", 2)
str_ints[22] = string_cell("22", 2)
str_ints[23] = string_cell("23", 2)
str_ints[24] = string_cell("24", 2)
str_ints[25] = string_cell("25", 2)
str_ints[26] = string_cell("26", 2)
str_ints[27] = string_cell("27", 2)
str_ints[28] = string_cell("28", 2)
str_ints[29] = string_cell("29", 2)
str_ints[30] = string_cell("30", 2)
str_ints[31] = string_cell("31", 2)
str_ints[32] = string_cell("32", 2)
str_ints[33] = string_cell("33", 2)
str_ints[34] = string_cell("34", 2)
str_ints[35] = string_cell("35", 2)
str_ints[36] = string_cell("36", 2)
str_ints[37] = string_cell("37", 2)
str_ints[38] = string_cell("38", 2)
str_ints[39] = string_cell("39", 2)
str_ints[40] = string_cell("40", 2)
str_ints[41] = string_cell("41", 2)
str_ints[42] = string_cell("42", 2)
str_ints[43] = string_cell("43", 2)
str_ints[44] = string_cell("44", 2)
str_ints[45] = string_cell("45", 2)
str_ints[46] = string_cell("46", 2)
str_ints[47] = string_cell("47", 2)
str_ints[48] = string_cell("48", 2)
str_ints[49] = string_cell("49", 2)
str_ints[50] = string_cell("50", 2)
str_ints[51] = string_cell("51", 2)
str_ints[52] = string_cell("52", 2)
str_ints[53] = string_cell("53", 2)
str_ints[54] = string_cell("54", 2)
str_ints[55] = string_cell("55", 2)
str_ints[56] = string_cell("56", 2)
str_ints[57] = string_cell("57", 2)
str_ints[58] = string_cell("58", 2)
str_ints[59] = string_cell("59", 2)
str_ints[60] = string_cell("60", 2)
str_ints[61] = string_cell("61", 2)
str_ints[62] = string_cell("62", 2)
str_ints[63] = string_cell("63", 2)
str_ints[64] = string_cell("64", 2)
str_ints[65] = string_cell("65", 2)
str_ints[66] = string_cell("66", 2)
str_ints[67] = string_cell("67", 2)
str_ints[68] = string_cell("68", 2)
str_ints[69] = string_cell("69", 2)
str_ints[70] = string_cell("70", 2)
str_ints[71] = string_cell("71", 2)
str_ints[72] = string_cell("72", 2)
str_ints[73] = string_cell("73", 2)
str_ints[74] = string_cell("74", 2)
str_ints[75] = string_cell("75", 2)
str_ints[76] = string_cell("76", 2)
str_ints[77] = string_cell("77", 2)
str_ints[78] = string_cell("78", 2)
str_ints[79] = string_cell("79", 2)
str_ints[80] = string_cell("80", 2)
str_ints[81] = string_cell("81", 2)
str_ints[82] = string_cell("82", 2)
str_ints[83] = string_cell("83", 2)
str_ints[84] = string_cell("84", 2)
str_ints[85] = string_cell("85", 2)
str_ints[86] = string_cell("86", 2)
str_ints[87] = string_cell("87", 2)
str_ints[88] = string_cell("88", 2)
str_ints[89] = string_cell("89", 2)
str_ints[90] = string_cell("90", 2)
str_ints[91] = string_cell("91", 2)
str_ints[92] = string_cell("92", 2)
str_ints[93] = string_cell("93", 2)
str_ints[94] = string_cell("94", 2)
str_ints[95] = string_cell("95", 2)
str_ints[96] = string_cell("96", 2)
str_ints[97] = string_cell("97", 2)
str_ints[98] = string_cell("98", 2)
str_ints[99] = string_cell("99", 2)
str_ints[100] = string_cell("100", 3)
str_ints[101] = string_cell("101", 3)
str_ints[102] = string_cell("102", 3)
str_ints[103] = string_cell("103", 3)
str_ints[104] = string_cell("104", 3)
str_ints[105] = string_cell("105", 3)
str_ints[106] = string_cell("106", 3)
str_ints[107] = string_cell("107", 3)
str_ints[108] = string_cell("108", 3)
str_ints[109] = string_cell("109", 3)
str_ints[110] = string_cell("110", 3)
str_ints[111] = string_cell("111", 3)
str_ints[112] = string_cell("112", 3)
str_ints[113] = string_cell("113", 3)
str_ints[114] = string_cell("114", 3)
str_ints[115] = string_cell("115", 3)
str_ints[116] = string_cell("116", 3)
str_ints[117] = string_cell("117", 3)
str_ints[118] = string_cell("118", 3)
str_ints[119] = string_cell("119", 3)
str_ints[120] = string_cell("120", 3)
str_ints[121] = string_cell("121", 3)
str_ints[122] = string_cell("122", 3)
str_ints[123] = string_cell("123", 3)
str_ints[124] = string_cell("124", 3)
str_ints[125] = string_cell("125", 3)
str_ints[126] = string_cell("126", 3)
str_ints[127] = string_cell("127", 3)
str_ints[128] = string_cell("128", 3)
str_ints[129] = string_cell("129", 3)
str_ints[130] = string_cell("130", 3)
str_ints[131] = string_cell("131", 3)
str_ints[132] = string_cell("132", 3)
str_ints[133] = string_cell("133", 3)
str_ints[134] = string_cell("134", 3)
str_ints[135] = string_cell("135", 3)
str_ints[136] = string_cell("136", 3)
str_ints[137] = string_cell("137", 3)
str_ints[138] = string_cell("138", 3)
str_ints[139] = string_cell("139", 3)
str_ints[140] = string_cell("140", 3)
str_ints[141] = string_cell("141", 3)
str_ints[142] = string_cell("142", 3)
str_ints[143] = string_cell("143", 3)
str_ints[144] = string_cell("144", 3)
str_ints[145] = string_cell("145", 3)
str_ints[146] = string_cell("146", 3)
str_ints[147] = string_cell("147", 3)
str_ints[148] = string_cell("148", 3)
str_ints[149] = string_cell("149", 3)
str_ints[150] = string_cell("150", 3)
str_ints[151] = string_cell("151", 3)
str_ints[152] = string_cell("152", 3)
str_ints[153] = string_cell("153", 3)
str_ints[154] = string_cell("154", 3)
str_ints[155] = string_cell("155", 3)
str_ints[156] = string_cell("156", 3)
str_ints[157] = string_cell("157", 3)
str_ints[158] = string_cell("158", 3)
str_ints[159] = string_cell("159", 3)
str_ints[160] = string_cell("160", 3)
str_ints[161] = string_cell("161", 3)
str_ints[162] = string_cell("162", 3)
str_ints[163] = string_cell("163", 3)
str_ints[164] = string_cell("164", 3)
str_ints[165] = string_cell("165", 3)
str_ints[166] = string_cell("166", 3)
str_ints[167] = string_cell("167", 3)
str_ints[168] = string_cell("168", 3)
str_ints[169] = string_cell("169", 3)
str_ints[170] = string_cell("170", 3)
str_ints[171] = string_cell("171", 3)
str_ints[172] = string_cell("172", 3)
str_ints[173] = string_cell("173", 3)
str_ints[174] = string_cell("174", 3)
str_ints[175] = string_cell("175", 3)
str_ints[176] = string_cell("176", 3)
str_ints[177] = string_cell("177", 3)
str_ints[178] = string_cell("178", 3)
str_ints[179] = string_cell("179", 3)
str_ints[180] = string_cell("180", 3)
str_ints[181] = string_cell("181", 3)
str_ints[182] = string_cell("182", 3)
str_ints[183] = string_cell("183", 3)
str_ints[184] = string_cell("184", 3)
str_ints[185] = string_cell("185", 3)
str_ints[186] = string_cell("186", 3)
str_ints[187] = string_cell("187", 3)
str_ints[188] = string_cell("188", 3)
str_ints[189] = string_cell("189", 3)
str_ints[190] = string_cell("190", 3)
str_ints[191] = string_cell("191", 3)
str_ints[192] = string_cell("192", 3)
str_ints[193] = string_cell("193", 3)
str_ints[194] = string_cell("194", 3)
str_ints[195] = string_cell("195", 3)
str_ints[196] = string_cell("196", 3)
str_ints[197] = string_cell("197", 3)
str_ints[198] = string_cell("198", 3)
str_ints[199] = string_cell("199", 3)
str_ints[200] = string_cell("200", 3)
str_ints[201] = string_cell("201", 3)
str_ints[202] = string_cell("202", 3)
str_ints[203] = string_cell("203", 3)
str_ints[204] = string_cell("204", 3)
str_ints[205] = string_cell("205", 3)
str_ints[206] = string_cell("206", 3)
str_ints[207] = string_cell("207", 3)
str_ints[208] = string_cell("208", 3)
str_ints[209] = string_cell("209", 3)
str_ints[210] = string_cell("210", 3)
str_ints[211] = string_cell("211", 3)
str_ints[212] = string_cell("212", 3)
str_ints[213] = string_cell("213", 3)
str_ints[214] = string_cell("214", 3)
str_ints[215] = string_cell("215", 3)
str_ints[216] = string_cell("216", 3)
str_ints[217] = string_cell("217", 3)
str_ints[218] = string_cell("218", 3)
str_ints[219] = string_cell("219", 3)
str_ints[220] = string_cell("220", 3)
str_ints[221] = string_cell("221", 3)
str_ints[222] = string_cell("222", 3)
str_ints[223] = string_cell("223", 3)
str_ints[224] = string_cell("224", 3)
str_ints[225] = string_cell("225", 3)
str_ints[226] = string_cell("226", 3)
str_ints[227] = string_cell("227", 3)
str_ints[228] = string_cell("228", 3)
str_ints[229] = string_cell("229", 3)
str_ints[230] = string_cell("230", 3)
str_ints[231] = string_cell("231", 3)
str_ints[232] = string_cell("232", 3)
str_ints[233] = string_cell("233", 3)
str_ints[234] = string_cell("234", 3)
str_ints[235] = string_cell("235", 3)
str_ints[236] = string_cell("236", 3)
str_ints[237] = string_cell("237", 3)
str_ints[238] = string_cell("238", 3)
str_ints[239] = string_cell("239", 3)
str_ints[240] = string_cell("240", 3)
str_ints[241] = string_cell("241", 3)
str_ints[242] = string_cell("242", 3)
str_ints[243] = string_cell("243", 3)
str_ints[244] = string_cell("244", 3)
str_ints[245] = string_cell("245", 3)
str_ints[246] = string_cell("246", 3)
str_ints[247] = string_cell("247", 3)
str_ints[248] = string_cell("248", 3)
str_ints[249] = string_cell("249", 3)
str_ints[250] = string_cell("250", 3)
str_ints[251] = string_cell("251", 3)
str_ints[252] = string_cell("252", 3)
str_ints[253] = string_cell("253", 3)
str_ints[254] = string_cell("254", 3)
str_ints[255] = string_cell("255", 3)
# [[[end]]]

# cdef void init_str_ints():
#     cdef int i
#     cdef Py_ssize_t z
#     cdef str pa
#     cdef char* a
#     cdef char* atemp
#     print("Starting string number table build")
#     for i in range(ARRAY_SIZE):
#         pa = str(i)
#         atemp = PyStr_AsUTF8AndSize(pa, &z)
#         a = <char*>malloc(sizeof(char) * z)
#         strcpy(a, atemp)
#         a[z] = "\0"
#         str_ints[i].string = a
#         str_ints[i].size = z
#     print("Finished building string number table")

# init_str_ints()



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

    if i < ARRAY_SIZE and i >= 0:
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


cdef int cmp_monosaccharide_count_cell(const void* lhs, const void* rhs) noexcept nogil:
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
                Py_INCREF(name)
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