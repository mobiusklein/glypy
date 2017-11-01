
cdef class CComposition(dict):
    cdef object _mass
    cdef tuple _mass_args
    cpdef CComposition clone(self)
    cpdef double calc_mass(self, int average=?, charge=?, dict mass_data=?) except -1
    cpdef _from_formula(self, str formula, dict mass_data)
    cpdef _from_dict(self, comp)
    cdef long getitem(self, str elem)
    cdef void setitem(self, str elem, long val)

    @staticmethod
    cdef CComposition _create(CComposition inst)


cdef: 
    str _atom
    str _formula
    str _isotope_string

    object isotope_pattern
    object formula_pattern

    cdef str _parse_isotope_string(str label, int* isotope_num)
    cdef str _make_isotope_string(str element_name, int isotope_num)


cpdef double calculate_mass(CComposition composition=?, str formula=?, int average=?, charge=?, mass_data=?) except -1
