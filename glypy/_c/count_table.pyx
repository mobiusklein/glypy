import sys
cimport cython
from cpython cimport PyErr_SetString, PyObject_Hash
from cpython.mem cimport PyMem_Malloc, PyMem_Free, PyMem_Realloc
from cpython.ref cimport PyObject, Py_INCREF, Py_DECREF, Py_XDECREF, Py_XINCREF
from cpython.int cimport PyInt_AsLong, PyInt_FromLong
from cpython.dict cimport PyDict_Next, PyDict_SetItem

try:
    from collections import Mapping, MutableMapping
except ImportError:
    from collections.abc import Mapping, MutableMapping


cdef int initialize_count_table_bin(count_table_bin* bin, size_t size):
    bin.cells = <count_table_bin_cell*>PyMem_Malloc(sizeof(count_table_bin_cell) * size)
    if bin.cells == NULL:
        return 1
    for i in range(size):
        bin.cells[i].key = NULL
    bin.used = 0
    bin.size = size


cdef void free_count_table_bin(count_table_bin* bin):
    for i in range(bin.used):
        if bin.cells[i].key != NULL:
            Py_XDECREF(bin.cells[i].key)
            bin.cells[i].key = NULL
    PyMem_Free(bin.cells)


cdef int count_table_bin_append(count_table_bin* bin, PyObject* key, long value):
    if bin.used == bin.size - 1:
        bin.cells = <count_table_bin_cell*>PyMem_Realloc(bin.cells, sizeof(count_table_bin_cell) * bin.size * 2)
        if bin.cells == NULL:
            return 1
        bin.size *= 2
    Py_XINCREF(key)
    bin.cells[bin.used].key = key
    bin.cells[bin.used].value = value
    bin.used += 1
    return 0


cdef int count_table_bin_find(count_table_bin* bin, PyObject* query, Py_ssize_t* cell_index):
    cdef:
        object query_obj
    query_obj = <object>query
    for i in range(bin.used):
        if bin.cells[i].key == NULL:
            continue
        Py_XINCREF(bin.cells[i].key)
        if (bin.cells[i].key == (<PyObject*>query_obj)) or ((<object>bin.cells[i].key) == (query_obj)):
            Py_XDECREF(bin.cells[i].key)
            cell_index[0] = i
            return 0
        Py_XDECREF(bin.cells[i].key)
    cell_index[0] = -1
    return 0


cdef count_table* make_count_table(size_t table_size, size_t bin_size):
    cdef:
        size_t i, j
        count_table* table
    table = <count_table*>PyMem_Malloc(sizeof(count_table))
    if table == NULL:
        PyErr_SetString(MemoryError, "Could not allocate memory for count_table")
        return NULL
    table.bins = <count_table_bin*>PyMem_Malloc(sizeof(count_table_bin) * table_size)
    if table.bins == NULL:
        PyErr_SetString(MemoryError, "Could not allocate memory for count_table")
        return NULL
    for i in range(table_size):
        initialize_count_table_bin(&(table.bins[i]), bin_size)
    table.size = table_size
    return table


cdef void free_count_table(count_table* table):
    for i in range(table.size):
        free_count_table_bin(&table.bins[i])
    PyMem_Free(table.bins)
    PyMem_Free(table)


@cython.cdivision(True)
cdef int count_table_find_bin(count_table* table, PyObject* query, Py_ssize_t* bin_index) except 1:
    try:
        bin_index[0] = PyObject_Hash(<object>query) % table.size
    except TypeError:
        PyErr_SetString(TypeError, "%r is not hashable" % (<object>query, ))
        return 1
    return 0


cdef int count_table_put(count_table* table, PyObject* key, long value) except 1:
    cdef:
        Py_ssize_t bin_index, cell_index
        int status
    if key == NULL:
        return 1
    Py_XINCREF(key)
    status = count_table_find_bin(table, key, &bin_index)
    if status != 0:
        return 1
    status = count_table_bin_find(&table.bins[bin_index], key, &cell_index)
    if cell_index == -1:
        status = count_table_bin_append(&table.bins[bin_index], key, value)
        Py_XDECREF(key)
        if status != 0:
            return 1
    else:
        table.bins[bin_index].cells[cell_index].value = value
        Py_XDECREF(key)
    return 0


cdef int count_table_del(count_table* table, PyObject* key, long* value) except 1:
    cdef:
        Py_ssize_t bin_index, cell_index
        int status
    Py_XINCREF(key)
    value[0] = 0
    status = count_table_find_bin(table, key, &bin_index)
    if status != 0:
        return 1
    status = count_table_bin_find(&table.bins[bin_index], key, &cell_index)
    if cell_index == -1:
        value[0] = 0
        Py_XDECREF(key)
        return 0
    else:
        value[0] = table.bins[bin_index].cells[cell_index].value
        Py_XDECREF(table.bins[bin_index].cells[cell_index].key)
        table.bins[bin_index].cells[cell_index].key = NULL
        table.bins[bin_index].cells[cell_index].value = 0
        Py_XDECREF(key)
        return 0


cdef int count_table_get(count_table* table, PyObject* key, long* value) except 1:
    cdef:
        Py_ssize_t bin_index, cell_index
        int status
    value[0] = 0
    Py_XINCREF(key)
    status = count_table_find_bin(table, key, &bin_index)
    if status != 0:
        return 1
    status = count_table_bin_find(&table.bins[bin_index], key, &cell_index)
    Py_XDECREF(key)
    if cell_index == -1:
        value[0] = 0
        return 0
    else:
        value[0] = table.bins[bin_index].cells[cell_index].value
        return 0


cdef int count_table_increment(count_table* table, PyObject* key, long value):
    cdef:
        Py_ssize_t bin_index, cell_index
        int status
    Py_XINCREF(key)
    status = count_table_find_bin(table, key, &bin_index)
    status = count_table_bin_find(&table.bins[bin_index], key, &cell_index)
    if cell_index == -1:
        status = count_table_bin_append(&table.bins[bin_index], key, value)
        Py_XDECREF(key)
        if status != 0:
            return 1
    else:
        table.bins[bin_index].cells[cell_index].value += value
        Py_XDECREF(key)
    return 0


cdef int count_table_decrement(count_table* table, PyObject* key, long value):
    cdef:
        Py_ssize_t bin_index, cell_index
        int status
    Py_XINCREF(key)
    status = count_table_find_bin(table, key, &bin_index)
    status = count_table_bin_find(&table.bins[bin_index], key, &cell_index)
    if cell_index == -1:
        status = count_table_bin_append(&table.bins[bin_index], key, -value)
        Py_XDECREF(key)
        if status != 0:
            return 1
    else:
        table.bins[bin_index].cells[cell_index].value -= value
        Py_XDECREF(key)
    return 0


cdef Py_ssize_t count_table_count(count_table* table):
    cdef:
        Py_ssize_t count
        size_t i, j
    count = 0
    for i in range(table.size):
        for j in range(table.bins[i].used):
            count += table.bins[i].cells[j].key != NULL
    return count


cdef list count_table_keys(count_table* table):
    cdef:
        list keys
        size_t i, j
    keys = []
    count = 0
    for i in range(table.size):
        for j in range(table.bins[i].used):
            if table.bins[i].cells[j].key != NULL:
                keys.append(<object>table.bins[i].cells[j].key)
    return keys


cdef list count_table_values(count_table* table):
    cdef:
        list values
        size_t i, j
    values = []
    count = 0
    for i in range(table.size):
        for j in range(table.bins[i].used):
            if table.bins[i].cells[j].key != NULL:
                values.append(PyInt_FromLong(table.bins[i].cells[j].value))
    return values


cdef list count_table_items(count_table* table):
    cdef:
        list items
        size_t i, j
    items = []
    for i in range(table.size):
        for j in range(table.bins[i].used):
            if table.bins[i].cells[j].key != NULL:
                items.append(
                    (<object>table.bins[i].cells[j].key,
                     PyInt_FromLong(table.bins[i].cells[j].value)))
    return items


cdef void count_table_add(count_table* table_a, count_table* table_b):
    cdef:
        size_t i, j
        long value
    for i in range(table_b.size):
        for j in range(table_b.bins[i].used):
            if table_b.bins[i].cells[j].key != NULL:
                count_table_get(table_b, table_b.bins[i].cells[j].key, &value)
                count_table_increment(table_a, table_b.bins[i].cells[j].key, value)


cdef void count_table_subtract(count_table* table_a, count_table* table_b):
    cdef:
        size_t i, j
        long value, temp
    for i in range(table_b.size):
        for j in range(table_b.bins[i].used):
            if table_b.bins[i].cells[j].key != NULL:
                count_table_get(table_b, table_b.bins[i].cells[j].key, &value)
                count_table_decrement(table_a, table_b.bins[i].cells[j].key, value)


cdef void count_table_scale(count_table* table, long value):
    cdef:
        size_t i, j
    for i in range(table.size):
        for j in range(table.bins[i].used):
            if table.bins[i].cells[j].key != NULL:
                table.bins[i].cells[j].value *= value


cdef void count_table_update(count_table* table_a, count_table* table_b):
    cdef:
        size_t i, j
        long value
    for i in range(table_b.size):
        for j in range(table_b.bins[i].used):
            if table_b.bins[i].cells[j].key != NULL:
                count_table_get(table_b, table_b.bins[i].cells[j].key, &value)
                count_table_put(table_a, table_b.bins[i].cells[j].key, value)


cdef void count_table_clear(count_table* table):
    cdef:
        size_t i, j
    for i in range(table.size):
        for j in range(table.bins[i].used):
            Py_XDECREF(table.bins[i].cells[j].key)
            table.bins[i].cells[j].key = NULL
            table.bins[i].cells[j].value = 0


cdef count_table* count_table_copy(count_table* table_a):
    cdef:
        count_table* dup
        size_t i, j
    dup = make_count_table(table_a.size, 2)
    for i in range(table_a.size):
        for j in range(table_a.bins[i].used):
            if table_a.bins[i].cells[j].key != NULL:
                count_table_bin_append(
                    &dup.bins[i], table_a.bins[i].cells[j].key,
                    table_a.bins[i].cells[j].value)
    return dup


cdef dict count_table_to_dict(count_table* table_a):
    cdef:
        dict dup
        size_t i, j
    dup = dict()
    for i in range(table_a.size):
        for j in range(table_a.bins[i].used):
            if table_a.bins[i].cells[j].key != NULL:
                PyDict_SetItem(
                    dup, <object>table_a.bins[i].cells[j].key,
                    PyInt_FromLong(table_a.bins[i].cells[j].value))
    return dup


cdef bint count_table_equals(count_table* table_a, count_table* table_b):
    cdef:
        size_t i, j
        long value
    for i in range(table_a.size):
        for j in range(table_a.bins[i].used):
            if table_a.bins[i].cells[j].key != NULL:
                count_table_get(table_b, table_a.bins[i].cells[j].key, &value)
                if table_a.bins[i].cells[j].value != value:
                    return False
    for i in range(table_b.size):
        for j in range(table_b.bins[i].used):
            if table_b.bins[i].cells[j].key != NULL:
                count_table_get(table_a, table_b.bins[i].cells[j].key, &value)
                if table_b.bins[i].cells[j].value != value:
                    return False
    return True


cdef class CountTableIterator(object):

    def __init__(self, table):
        self.table = (<CountTable>table).table
        self.bin_index = 0
        self.cell_index = 0
        self.next_key = NULL
        self.next_value = 0
        self.advance()

    @staticmethod
    cdef CountTableIterator _create(CountTable table):
        cdef CountTableIterator it = CountTableIterator.__new__(CountTableIterator)
        it.table = table.table
        it.bin_index = 0
        it.cell_index = 0
        it.next_key = NULL
        it.next_value = 0
        it.advance()
        return it

    cdef bint has_more(self):
        return self.next_key != NULL

    cdef int _locate_next_value(self, PyObject** key, long* value):
        if not self.bin_index < self.table.size:
            key[0] = NULL
            return 1
        # start from the previous iteration's bin
        for i in range(self.cell_index, self.table.bins[self.bin_index].used):
            self.cell_index = i
            if self.table.bins[self.bin_index].cells[i].key != NULL:
                key[0] = self.table.bins[self.bin_index].cells[i].key
                value[0] = self.table.bins[self.bin_index].cells[i].value
                self.cell_index += 1
                return 0
        # move onto the next bins
        for i in range(self.bin_index + 1, self.table.size):
            self.bin_index = i
            for j in range(0, self.table.bins[i].used):
                self.cell_index = j
                if self.table.bins[i].cells[j].key != NULL:
                    key[0] = self.table.bins[self.bin_index].cells[j].key
                    value[0] = self.table.bins[self.bin_index].cells[j].value
                    self.cell_index += 1
                    return 0
        self.bin_index = self.table.size
        key[0] = NULL
        return 1

    cdef int get_next_value(self, PyObject** key, long* value):
        cdef:
            PyObject* next_key
            long next_value
            int status
        if self.next_key == NULL:
            key[0] = NULL
            value[0] = 0
            return 1
        else:
            key[0] = self.next_key
            value[0] = self.next_value
            status = self._locate_next_value(&next_key, &next_value)
            self.next_key = next_key
            self.next_value = next_value
        return 0

    cdef void advance(self):
        cdef:
            PyObject* next_key
            long next_value
            int status
        status = self._locate_next_value(&next_key, &next_value)
        self.next_key = next_key
        self.next_value = next_value

    def __iter__(self):
        return self

    def __next__(self):
        cdef:
            PyObject* key
            PyObject* next_key
            long value, next_value
            int status
        if self.next_key == NULL:
            raise StopIteration()
        else:
            key = self.next_key
            value = self.next_value
            status = self._locate_next_value(&next_key, &next_value)
            self.next_key = next_key
            self.next_value = next_value
        return (<object>key, value)

    def _test(self):
        cdef:
            PyObject* key
            long count
            int pos, j
            list result
        pos = 0
        j = 0
        result = []
        while self.has_more():
            pos = self.get_next_value(&key, &count)
            if pos != 0:
                break
            result.append((<object>key, count))
            j += 1
        return result



cdef class CountTable(object):

    @staticmethod
    cdef CountTable _create():
        cdef CountTable inst
        inst = CountTable.__new__(CountTable)
        inst._initialize_table()
        return inst

    @staticmethod
    cdef CountTable _create_from(CountTable table):
        cdef CountTable inst
        inst = CountTable.__new__(CountTable)
        inst.table = count_table_copy(table.table)
        return inst        

    def __init__(self, obj=None, **kwargs):
        self._initialize_table()
        if obj is not None:
            self.update(obj)
        if kwargs:
            self.update(kwargs)

    cdef int _initialize_table(self) except 1:
        self.table = make_count_table(6, 2)
        if self.table == NULL:
            return 1
        return 0

    cpdef update(self, obj):
        if isinstance(obj, CountTable):
            self._update_from_count_table(<CountTable>obj)
        if isinstance(obj, dict):
            self._update_from_dict(<dict>obj)
        else:
            for k, v in obj.items():
                self[k] = v

    cdef void _update_from_dict(self, dict d):
        cdef:
            PyObject *key
            PyObject *value
            Py_ssize_t pos
        pos = 0
        while PyDict_Next(d, &pos, &key, &value):
            self.setitem(<object>key, PyInt_AsLong(int(<object>value)))

    cdef void _update_from_count_table(self, CountTable other):
        count_table_update(self.table, other.table)

    cpdef CountTable copy(self):
        cdef CountTable inst = CountTable._create_from(self)
        return inst

    def __reduce__(self):
        return self.__class__, (self._to_dict(), )

    def __dealloc__(self):
        free_count_table(self.table)

    def __getitem__(self, object key):
        cdef long value = self.getitem(key)
        return PyInt_FromLong(value)

    def __setitem__(self, object key, object value):
        self.setitem(key, PyInt_AsLong(int(value)))

    def __delitem__(self, object key):
        self.pop(key)

    def __len__(self):
        return count_table_count(self.table)

    def __contains__(self, key):
        cdef long value = self.getitem(key)
        return value != 0

    def __iter__(self):
        return iter(self.keys())

    def __add__(self, other):
        cdef:
            CountTable dup
        if isinstance(other, CountTable):
            if isinstance(self, CountTable):
                dup = (<CountTable>self).copy()
                dup._add_from(<CountTable>other)
                return dup
            elif isinstance(self, dict):
                dup = (<CountTable>other).copy()
                dup._add_from_dict(<dict>self)
                return dup
            else:
                self = CountTable(self)
                dup = (<CountTable>other).copy()
                dup._add_from(self)
                return dup
        else:
            if isinstance(other, dict):
                dup = (<CountTable>self).copy()
                dup._add_from_dict(<dict>other)
                return dup
            else:
                dup = (<CountTable>self).copy()
                dup._add_from(CountTable(other))
                return dup

    def __iadd__(self, other):
        if isinstance(other, CountTable):
            self._add_from(<CountTable>other)
        elif isinstance(other, dict):
            self._add_from_dict(<dict>other)
        else:
            self._add_from(CountTable(other))
        return self

    def __sub__(self, other):
        cdef:
            CountTable dup
        if isinstance(other, CountTable):
            if isinstance(self, CountTable):
                dup = (<CountTable>self).copy()
                dup._subtract_from(<CountTable>other)
                return dup
            elif isinstance(self, dict):
                dup = (<CountTable>other).copy()
                dup._subtract_from_dict(<dict>self)
                return dup
            else:
                self = CountTable(self)
                dup = (<CountTable>other).copy()
                dup._subtract_from(self)
                return dup
        else:
            if isinstance(other, dict):
                dup = (<CountTable>self).copy()
                dup._subtract_from_dict(<dict>other)
                return dup
            else:
                dup = (<CountTable>self).copy()
                dup._subtract_from(CountTable(other))
                return dup

    def __isub__(self, other):
        if isinstance(other, CountTable):
            self._subtract_from(<CountTable>other)
        elif isinstance(other, dict):
            self._subtract_from_dict(<dict>other)
        else:
            self._subtract_from(CountTable(other))
        return self

    def __mul__(self, other):
        cdef:
            CountTable dup
        if isinstance(other, CountTable):
            dup = (<CountTable>other).copy()
            dup._scale_by(self)
            return dup
        elif isinstance(self, CountTable):
            dup = (<CountTable>self).copy()
            dup._scale_by(other)
            return dup
        else:
            return NotImplemented

    def __imul__(self, other):
        self._scale_by(other)
        return self

    def __neg__(self):
        return self._scale_by(-1)

    cdef void _add_from(self, CountTable other):
        count_table_add(self.table, other.table)

    cdef int _add_from_dict(self, dict other) except 1:
        cdef:
            PyObject *key
            PyObject *value
            Py_ssize_t pos
            long v
        pos = 0
        try:
            while PyDict_Next(other, &pos, &key, &value):
                v = PyInt_AsLong(int(<object>value))
                self.setitem(
                    <object>key,
                    v + self.getitem(<object>key))
        except TypeError:
            PyErr_SetString(TypeError, "%r must be a number" % (<object>value, ))
            return 1
        return 0

    cdef void _subtract_from(self, CountTable other):
        count_table_subtract(self.table, other.table)

    cdef int _subtract_from_dict(self, dict other) except 1:
        cdef:
            PyObject *key
            PyObject *value
            Py_ssize_t pos
            long v
        pos = 0
        try:
            while PyDict_Next(other, &pos, &key, &value):
                v = PyInt_AsLong(int(<object>value))
                self.setitem(
                    <object>key,
                    v - self.getitem(<object>key))
        except TypeError:
            PyErr_SetString(TypeError, "%r must be a number" % (<object>value, ))
            return 1
        return 0

    cdef void _scale_by(self, long value):
        count_table_scale(self.table, value)

    cpdef dict _to_dict(self):
        return count_table_to_dict(self.table)

    def __repr__(self):
        return "{self.__class__.__name__}({content})".format(
            self=self, content=self._to_dict())

    cdef bint equal_to(self, CountTable other):
        return count_table_equals(self.table, other.table)

    def __eq__(self, other):
        if isinstance(self, CountTable):
            if isinstance(other, CountTable):
                return count_table_equals(self.table, (<CountTable>other).table)
            else:
                return self._to_dict() == other
        else:
            return self == other._to_dict()

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        cdef:
            Py_hash_t hash_value
            size_t i, j

        hash_value = 0
        for i in range(self.table.size):
            for j in range(self.table.bins[i].used):
                if self.table.bins[i].cells[j].key != NULL:
                    hash_value ^= hash(<object>self.table.bins[i].cells[j].key)
                    hash_value ^= self.table.bins[i].cells[j].value
        hash_value ^= (hash_value >> 11) ^ (hash_value >> 25)
        hash_value = hash_value * 69069 + 907133923
        return hash_value

    cpdef list keys(self):
        return count_table_keys(self.table)

    cpdef list values(self):
        return count_table_values(self.table)

    cpdef list items(self):
        return count_table_items(self.table)

    cpdef clear(self):
        count_table_clear(self.table)

    cpdef setdefault(self, key, value):
        if self.getitem(key) == 0:
            self.setitem(key, value)

    cdef long getitem(self, object key) except *:
        cdef long value
        cdef int status 
        cdef PyObject* pkey = <PyObject*>key
        Py_INCREF(key)
        status = count_table_get(self.table, pkey, &value)
        Py_DECREF(key)
        if status != 0:
            return 0
        return value        

    cdef int setitem(self, object key, long value) except 1:
        cdef PyObject* pkey = <PyObject*>key
        Py_INCREF(key)
        cdef int status = count_table_put(self.table, pkey, value)
        Py_DECREF(key)
        if status != 0:
            return 1
        return 0

    cdef long delitem(self, object key) except *:
        cdef PyObject* pkey = <PyObject*>key
        cdef long value
        Py_INCREF(key)
        cdef int status = count_table_del(self.table, pkey, &value)
        Py_DECREF(key)
        if status != 0:
            return 0
        return value

    cdef void increment(self, object key, long value):
        cdef PyObject* pkey = <PyObject*>key
        Py_INCREF(key)
        cdef int status = count_table_increment(self.table, pkey, value)
        Py_DECREF(key)

    cdef void decrement(self, object key, long value):
        cdef PyObject* pkey = <PyObject*>key
        Py_INCREF(key)
        cdef int status = count_table_decrement(self.table, pkey, value)
        Py_DECREF(key)

    cpdef long pop(self, object key, object default=None) except *:
        cdef long value = self.delitem(key)
        if value == 0:
            return default
        return PyInt_FromLong(value)

    cpdef get(self, object key, object default=None):
        cdef long value = self.getitem(key)
        if value == 0:
            return default
        return PyInt_FromLong(value)


Mapping.register(CountTable)
MutableMapping.register(CountTable)


def main():
    cdef long val
    cdef str string
    cdef object pint
    cdef PyObject* key
    string = "spam"
    pint = 12
    key = <PyObject*>string
    cdef count_table* table = make_count_table(6, 2)
    count_table_put(table, key, 42)
    count_table_get(table, key, &val)
    assert val == 42
    key = <PyObject*>pint
    count_table_put(table, key, 67)
    count_table_del(table, key, &val)
    assert val == 67
    val = 1
    count_table_get(table, key, &val)
    free_count_table(table)
    pytable = CountTable()
    Py_INCREF(<object>key)
    pytable[<object>key] = 2
    pytable[int(12)] = 11
    pytable[255.2] = 99
    assert len(pytable) == 2
    del pytable
    pytable = CountTable._create()
    print pytable
