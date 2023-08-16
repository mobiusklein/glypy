# cython: boundscheck=False

# Credit to Pyteomics - http://pythonhosted.org/pyteomics - for majority of design
import re
from .mass_dict import nist_mass
from .base import ChemicalCompositionError, composition_factory

cimport cython

from cpython.ref cimport PyObject
from cpython.dict cimport (PyDict_GetItem, PyDict_SetItem, PyDict_Next,
                           PyDict_Keys, PyDict_Update, PyDict_DelItem, PyDict_Size)
from cpython.int cimport PyInt_AsLong, PyInt_Check, PyInt_FromLong

from glypy.composition.compat cimport PyStr_Format

from cpython.float cimport PyFloat_AsDouble
from cpython.tuple cimport PyTuple_GetItem
from cpython.list cimport PyList_GET_ITEM, PyList_GET_SIZE

# Forward Declaration
cdef:
    str _atom = r'([A-Z][a-z+]*)(?:\[(\d+)\])?([+-]?\d+)?'
    str _formula = r'^({})*$'.format(_atom)
    str _isotope_string = r'^([A-Z][a-z+]*)(?:\[(\d+)\])?$'

    object isotope_pattern = re.compile(_isotope_string)
    object formula_pattern = re.compile(_formula)


@cython.boundscheck(False)
cdef str _parse_isotope_string(str label, int* isotope_num):
    cdef:
        # int isotope_num = 0
        int i = 0
        int in_bracket = False
        # str element_name
        str current
        list name_parts = []
        list num_parts = []
        #Isotope result
    for i in range(len(label)):
        current = label[i]
        if in_bracket:
            if current == "]":
                break
            num_parts.append(current)
        elif current == "[":
            in_bracket = True
        else:
            name_parts.append(current)
    element_name = (''.join(name_parts))
    if len(num_parts) > 0:
        isotope_num[0] = (int(''.join(num_parts)))
    else:
        isotope_num[0] = 0
    return element_name

def parse_isotope_string(str label):
    cdef:
        int iso
        str elem
    elem = _parse_isotope_string(label, &iso)
    return (elem, iso)

cdef str _make_isotope_string(str element_name, int isotope_num):
    """Form a string label for an isotope."""
    cdef:
        tuple parts
    if isotope_num == 0:
        return element_name
    else:
        parts = (element_name, isotope_num)
        return <str>PyStr_Format('%s[%d]', parts)


cpdef CComposition composition_sum(list compositions):
    cdef:
        size_t i, n
        CComposition accumulator, current
        str elem
        long cnt
        PyObject *pkey
        PyObject *pvalue
        Py_ssize_t ppos = 0

    n = PyList_GET_SIZE(compositions)
    if n == 1:
        accumulator = <CComposition>PyList_GET_ITEM(compositions, 0)
        return accumulator.copy()
    elif n == 0:
        accumulator = CComposition()
        return accumulator
    else:
        accumulator = (<CComposition>PyList_GET_ITEM(compositions, 0)).copy()
        for i in range(1, n):
            current = <CComposition>PyList_GET_ITEM(compositions, i)
            ppos = 0

            while(PyDict_Next(current, &ppos, &pkey, &pvalue)):
                elem = <str>pkey
                cnt = accumulator.getitem(elem)
                accumulator.setitem(elem, cnt + PyInt_AsLong(<object>pvalue))
        return accumulator


@cython.c_api_binop_methods(True)
cdef class CComposition(dict):
    """A Composition object stores a chemical composition of a
    substance. Basically it is a dict object, in which keys are the names
    of chemical elements and values contain integer numbers of
    corresponding atoms in a substance.

    The main improvement over dict is that Composition objects allow
    addition and subtraction.

    If ``formula`` is not specified, the constructor will look at the first
    positional argument and try to build the object from it. Without
    positional arguments, a Composition will be constructed directly from
    keyword arguments.

    Parameters
    ----------
    formula : str, optional
        A string with a chemical formula. All elements must be present in
        `mass_data`.
    mass_data : dict, optional
        A dict with the masses of chemical elements (the default
        value is :py:data:`nist_mass`). It is used for formulae parsing only.
    """

    sum = staticmethod(composition_sum)

    def __str__(self):   # pragma: no cover
        return 'Composition({})'.format(dict.__repr__(self))

    def __repr__(self):  # pragma: no cover
        return str(self)

    cdef CComposition add_from(self, other):
        cdef:
            str elem
            long cnt
            PyObject *pkey
            PyObject *pvalue
            Py_ssize_t ppos = 0

        while(PyDict_Next(other, &ppos, &pkey, &pvalue)):
            elem = <str>pkey
            cnt = self.getitem(elem)
            self.setitem(elem, cnt + PyInt_AsLong(<object>pvalue))

        self._mass_args = None
        return self

    def __iadd__(CComposition self, other):
        self.add_from(other)
        return self


    def __add__(self, other):
        cdef:
            CComposition result
        if not isinstance(self, CComposition):
            other, self = self, other
        result = CComposition._create(self)
        result.add_from(other)
        return result

    cdef CComposition subtract_from(self, other):
        cdef:
            str elem
            long cnt
            PyObject *pkey
            PyObject *pvalue
            Py_ssize_t ppos = 0

        while(PyDict_Next(other, &ppos, &pkey, &pvalue)):
            elem = <str>pkey
            cnt = self.getitem(elem)
            self.setitem(elem, cnt - PyInt_AsLong(<object>pvalue))

        self._mass_args = None
        return self

    def __isub__(self, other):
        return self.subtract_from(other)

    def __sub__(self, other):
        cdef:
            CComposition result
        if not isinstance(self, CComposition):
            self = CComposition(self)
        result = CComposition._create(self)
        result.subtract_from(other)
        return result

    def __reduce__(self):
        return composition_factory, (list(self),), self.__getstate__()

    def __getstate__(self):
        return dict(self)

    def __setstate__(self, d):
        self._from_dict(d)
        self._mass = None
        self._mass_args = None

    def __mul__(self, other):
        cdef:
            CComposition prod
            int rep, v
            str k
            PyObject *pkey
            PyObject *pvalue
            Py_ssize_t ppos = 0

        if isinstance(other, CComposition):
            self, other = other, self

        # prod = self.clone()
        # prod = CComposition()
        prod = CComposition._create(None)

        if not isinstance(other, int):
            raise ChemicalCompositionError(
                'Cannot multiply Composition by non-integer', other)

        rep = PyInt_AsLong(other)
        while(PyDict_Next(self, &ppos, &pkey, &pvalue)):
            k = <str>pkey
            v = PyInt_AsLong(<object>pvalue)
            prod.setitem(k, v * rep)

        return prod

    def __imul__(self, other):
        cdef:
            str k
            int v, rep
        if not isinstance(other, int):
            raise ChemicalCompositionError(
                'Cannot multiply Composition by non-integer',
                other)
        rep = other
        for k, v in list(self.items()):
            self.setitem(k, v * rep)
        return self

    def __eq__(self, other):
        cdef:
            Py_ssize_t self_size, other_size
            PyObject *pkey
            PyObject *pvalue
            PyObject *pinterm
            long self_value, other_value
        if not isinstance(other, dict):
            return False
        self_size = PyDict_Size(self)
        other_size = PyDict_Size(other)

        # if one dict is larger, it may contain more keys which may be zero,
        # which will be ignored.
        if other_size > self_size:
            other_size = 0
            while(PyDict_Next(other, &other_size, &pkey, &pvalue)):
                other_value = PyInt_AsLong(<object>pvalue)
                pinterm = PyDict_GetItem(self, <object>pkey)
                if pinterm == NULL:
                    self_value = 0
                else:
                    self_value = PyInt_AsLong(<object>pinterm)
                if self_value != other_value:
                    return False
        else:
            self_size = 0
            while(PyDict_Next(self, &self_size, &pkey, &pvalue)):
                self_value = PyInt_AsLong(<object>pvalue)
                pinterm = PyDict_GetItem(other, <object>pkey)
                if pinterm == NULL:
                    other_value = 0
                else:
                    other_value = PyInt_AsLong(<object>pinterm)
                if self_value != other_value:
                    return False
        return True

    def __ne__(self, other):
        cdef:
            Py_ssize_t self_size, other_size
            PyObject *pkey
            PyObject *pvalue
            PyObject *pinterm
            long self_value, other_value
        if not isinstance(other, dict):
            return True
        self_size = PyDict_Size(self)
        other_size = PyDict_Size(other)
        if other_size > self_size:
            other_size = 0
            while(PyDict_Next(other, &other_size, &pkey, &pvalue)):
                other_value = PyInt_AsLong(<object>pvalue)
                pinterm = PyDict_GetItem(self, <object>pkey)
                if pinterm == NULL:
                    self_value = 0
                else:
                    self_value = PyInt_AsLong(<object>pinterm)
                if self_value != other_value:
                    return True
        else:
            self_size = 0
            while(PyDict_Next(self, &self_size, &pkey, &pvalue)):
                self_value = PyInt_AsLong(<object>pvalue)
                pinterm = PyDict_GetItem(other, <object>pkey)
                if pinterm == NULL:
                    other_value = 0
                else:
                    other_value = PyInt_AsLong(<object>pinterm)
                if self_value != other_value:
                    return True
        return False

    def __neg__(self):
        return self * -1

    @staticmethod
    cdef CComposition _create(CComposition inst):
        cdef:
            CComposition result

        result = CComposition.__new__(CComposition)
        result._mass = None
        result._mass_args = None
        if inst is None:
            return result
        PyDict_Update(result, inst)
        return result


    # Override the default behavior, if a key is not present
    # do not initialize it to 0.
    def __missing__(self, str key):
        return 0

    def __setitem__(self, str key, object value):
        cdef long int_value = PyInt_AsLong(round(value))
        if int_value:  # Will not occur on 0 as 0 is falsey AND an integer
            self.setitem(key, int_value)
        elif key in self:
            del self[key]
        self._mass_args = None

    def copy(self):
        cdef CComposition dup
        dup = CComposition._create(None)
        PyDict_Update(dup, self)
        return dup

    cdef long getitem(self, str elem):
        cdef:
            PyObject* resobj
            long count
        resobj = PyDict_GetItem(self, elem)
        if (resobj == NULL):
            return 0
        count = PyInt_AsLong(<object>resobj)
        return count

    cdef void setitem(self, str elem, long val):
        PyDict_SetItem(self, elem, val)
        self._mass_args = None

    cpdef CComposition clone(self):
        return CComposition._create(self)

    def update(self, *args, **kwargs):
        dict.update(self, *args, **kwargs)
        self._mass_args = None

    @cython.boundscheck(False)
    cpdef _from_formula(self, str formula, dict mass_data):
        cdef:
            str elem, isotope, number
        if '(' in formula:
            self._from_formula_parens(formula, mass_data)
        elif not formula_pattern.match(formula):
            raise ChemicalCompositionError('Invalid formula: ' + formula)
        else:
            for elem, isotope, number in re.findall(_atom, formula):
                if not elem in mass_data:
                    raise ChemicalCompositionError('Unknown chemical element: ' + elem)
                self[_make_isotope_string(elem, int(isotope) if isotope else 0)
                        ] += int(number) if number else 1

    @cython.boundscheck(True)
    def _from_formula_parens(self, formula, mass_data):
        # Parsing a formula backwards.
        cdef:
            Py_ssize_t prev_chem_symbol_start, i
            int seek_mode, group_coef
            str parse_stack
            list resolve_stack


        prev_chem_symbol_start = len(formula)
        i = len(formula) - 1

        seek_mode = 0
        parse_stack = ""
        resolve_stack = []
        group_coef = 1

        while i >= 0:
            if seek_mode < 1:
                if (formula[i] == ")"):
                    seek_mode += 1
                    if i + 1 == prev_chem_symbol_start:
                        group_coef = 1
                    elif formula[i + 1].isdigit():
                        group_coef = int(formula[i + 1:prev_chem_symbol_start])
                    i -= 1
                    continue
                # Read backwards until a non-number character is met.
                if (formula[i].isdigit() or formula[i] == '-'):
                    i -= 1
                    continue

                else:
                    # If the number of atoms is omitted then it is 1.
                    if i + 1 == prev_chem_symbol_start:
                        num_atoms = 1
                    else:
                        try:
                            num_atoms = int(formula[i + 1:prev_chem_symbol_start])
                        except ValueError:
                            raise ChemicalCompositionError(
                                'Badly-formed number of atoms: %s' % formula)

                    # Read isotope number if specified, else it is undefined (=0).
                    if formula[i] == ']':
                        brace_pos = formula.rfind('[', 0, i)
                        if brace_pos == -1:
                            raise ChemicalCompositionError(
                                'Badly-formed isotope number: %s' % formula)
                        try:
                            isotope_num = int(formula[brace_pos + 1:i])
                        except ValueError:
                            raise ChemicalCompositionError(
                                'Badly-formed isotope number: %s' % formula)
                        i = brace_pos - 1
                    else:
                        isotope_num = 0

                    # Match the element name to the mass_data.
                    element_found = False
                    # Sort the keys from longest to shortest to workaround
                    # the overlapping keys issue
                    for element_name in sorted(mass_data, key=len, reverse=True):
                        if formula.endswith(element_name, 0, i + 1):
                            isotope_string = _make_isotope_string(
                                element_name, isotope_num)
                            self[isotope_string] += num_atoms
                            i -= len(element_name)
                            prev_chem_symbol_start = i + 1
                            element_found = True
                            break

                    if not element_found:
                        raise ChemicalCompositionError(
                            'Unknown chemical element in the formula: %s' % formula)
            else:
                ch = formula[i]
                parse_stack += ch
                i -= 1
                if(ch == "("):
                    seek_mode -= 1
                    if seek_mode == 0:

                        resolve_stack.append(Composition(
                                             # Omit the last character, then reverse the parse
                                             # stack string.
                                             formula=parse_stack[:-1][::-1],
                                             mass_data=mass_data)
                                             * group_coef)
                        prev_chem_symbol_start = i + 1
                        seek_mode = False
                        parse_stack = ""
                elif(formula[i] == ")"):
                    seek_mode += 1
                else:
                    # continue to accumulate tokens
                    pass

        # Unspool the resolve stack, adding together the chunks
        # at this level. __add__ operates immutably, so must manually
        # loop through each chunk.
        for chunk in resolve_stack:
            for elem, cnt in chunk.items():
                self[elem] += cnt

    cpdef _from_dict(self, comp):
        '''
        Directly overwrite this object's keys with the values in
        `comp` without checking their type.
        '''
        PyDict_Update(self, comp)

    cpdef double calc_mass(self, int average=False, charge=None, dict mass_data=nist_mass) except -1:
        cdef object mdid
        mdid = id(mass_data)
        if self._mass_args is not None and average is self._mass_args[0]\
                and charge == self._mass_args[1] and mdid == self._mass_args[2]:
            return self._mass
        else:
            self._mass_args = (average, charge, mdid)
            self._mass = calculate_mass(composition=self, average=average, charge=charge, mass_data=mass_data)
            return self._mass

    property mass:
        def __get__(self):
            return self.calc_mass()

    def __init__(self, args=None, formula=None, dict mass_data=None, **kwargs):
        #dict.__init__(self)
        if mass_data is None:
            mass_data = nist_mass

        if formula is not None:
            self._from_formula(formula, mass_data)

        elif args:
            if isinstance(args, dict):
                self._from_dict(args)
            elif isinstance(args, str):
                try:
                    self._from_formula(args, mass_data)
                except ChemicalCompositionError:
                    raise ChemicalCompositionError(
                        'Could not create a Composition object from '
                        'string: "{}": not a '
                        'formula'.format(args))
        else:
            self._from_dict(kwargs)
        self._mass = None
        self._mass_args = None

Composition = CComposition


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef double calculate_mass(CComposition composition=None, str formula=None, int average=False, charge=None, mass_data=None) except -1:
    """Calculates the monoisotopic mass of a chemical formula or CComposition object.

    Parameters
    ----------
    composition : CComposition
        A Composition object with the elemental composition of a substance. Exclusive with `formula`
    formula: str
        A string describing a chemical composition. Exclusive with `composition`
    average : bool, optional
        If :py:const:`True` then the average mass is calculated. Note that mass
        is not averaged for elements with specified isotopes. Default is
        :py:const:`False`.
    charge : int, optional
        If not 0 then m/z is calculated: the mass is increased
        by the corresponding number of proton masses and divided
        by z.
    mass_data : dict, optional
        A dict with the masses of the chemical elements (the default
        value is :py:data:`nist_mass`).

    Returns
    -------
        mass : float
    """
    cdef:
        int old_charge, isotope_num, isotope, quantity
        double mass, isotope_mass, isotope_frequency
        long _charge
        str isotope_string, element_name
        dict mass_provider
        list key_list
        PyObject* interm
        Py_ssize_t iter_pos = 0

    if mass_data is None:
        mass_provider = nist_mass
    else:
        mass_provider = mass_data

    if composition is None:
        if formula is not None:
            composition = CComposition(formula)
        else:
            raise ChemicalCompositionError("Must provide a composition or formula argument")
    else:
        if formula is not None:
            raise ChemicalCompositionError("Must provide a composition or formula argument, but not both")

    # Get charge.
    if charge is None:
        charge = composition.getitem('H+')
    else:
        if charge != 0 and composition.getitem('H+') != 0:
            raise ChemicalCompositionError("Charge is specified both by the number of protons and parameters")
    _charge = PyInt_AsLong(charge)
    old_charge = composition.getitem('H+')
    composition.setitem('H+', charge)

    # Calculate mass.
    mass = 0.0
    key_list = PyDict_Keys(composition)
    for iter_pos in range(len(key_list)):
        isotope_string = <str>PyList_GET_ITEM(key_list, iter_pos)
        # element_name, isotope_num = _parse_isotope_string(isotope_string)
        element_name = _parse_isotope_string(isotope_string, &isotope_num)

        # Calculate average mass if required and the isotope number is
        # not specified.
        if (not isotope_num) and average:
            for isotope in mass_provider[element_name]:
                if isotope != 0:
                    quantity = <int>composition.getitem(element_name)
                    isotope_mass = <double>mass_provider[element_name][isotope][0]
                    isotope_frequency = <double>mass_provider[element_name][isotope][1]

                    mass += quantity * isotope_mass * isotope_frequency
        else:
            interim = PyDict_GetItem(mass_provider, element_name)
            interim = PyDict_GetItem(<dict>interim, isotope_num)
            isotope_mass = PyFloat_AsDouble(<object>PyTuple_GetItem(<tuple>interim, 0))

            mass += (composition.getitem(isotope_string) * isotope_mass)

    # Calculate m/z if required.
    if _charge != 0:
        mass /= abs(_charge)

    if old_charge != 0:
        composition.setitem('H+', old_charge)
    else:
        PyDict_DelItem(composition, "H+")
    return mass
