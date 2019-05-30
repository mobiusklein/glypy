from glypy._c.structure.base cimport SaccharideBase
from glypy.composition.ccomposition cimport CComposition
from glypy.utils.cenum cimport EnumValue


cdef class Monosaccharide(SaccharideBase):
    cdef:
        public long id
        public EnumValue _anomer
        public EnumValue _superclass
        public int ring_start
        public int ring_end
        public object links
        public object substituent_links
        public object modifications
        public CComposition composition
        public object _reducing_end
        public object _degree
        public bint _checked_for_reduction


