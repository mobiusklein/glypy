cimport cython


@cython.freelist(100000)
cdef class MoleculeBase(object):
    node_type = object()

    def __cinit__(self, *args, **kwargs):
        self._order = 0

    cpdef order(self, deep=False):
        '''
        Return the "graph theory" order of this molecule

        Returns
        -------
        int
        '''
        return self._order

    cpdef has_undefined_linkages(self):
        return True

    def copy(self, *args, **kwargs):
        return self.clone(*args, **kwargs)


cdef class SaccharideBase(MoleculeBase):
    node_type = object()


cdef class SaccharideCollection(SaccharideBase):
    cpdef _derivatized(self, substituent, id_base):
        pass

    cpdef _strip_derivatization(self):
        pass


cdef class SubstituentBase(MoleculeBase):
    node_type = object()


cdef class ModificationBase(MoleculeBase):
    node_type = object()