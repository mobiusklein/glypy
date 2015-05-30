
cdef class CComposition(dict):
    cpdef CComposition clone(self)
    cpdef double calc_mass(self, int average=?, charge=?, dict mass_data=?) except -1
    cpdef _from_formula(self, str formula, dict mass_data)
    cpdef _from_dict(self, comp)


cdef: 
    dict std_mol_comp
    str _atom
    str _formula
    str _isotope_string

    object isotope_pattern
    object formula_pattern

    cdef inline tuple _parse_isotope_string(str label)
    cdef inline str _make_isotope_string(str element_name, int isotope_num)


cpdef inline double calculate_mass(CComposition composition=?, str formula=?, int average=?, charge=?, mass_data=?) except -1
