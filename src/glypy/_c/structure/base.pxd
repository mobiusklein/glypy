

cdef class MoleculeBase(object):

    cdef:
        public int _order

    cpdef order(self, deep=*)
    cpdef has_undefined_linkages(self)


cdef class SaccharideBase(MoleculeBase):
    pass


cdef class SaccharideCollection(SaccharideBase):
    cpdef _derivatized(self, substituent, id_base)

    cpdef _strip_derivatization(self)


cdef class SubstituentBase(MoleculeBase):
    pass


cdef class ModificationBase(MoleculeBase):
    pass
