#cython: boundscheck=False, debug=True, profile=True

# Credit to Pyteomics - http://pythonhosted.org/pyteomics - for majority of design
import re
from .mass_dict import nist_mass
from .base import ChemicalCompositionError, composition_factory

cimport cython
from cpython.ref cimport PyObject
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem, PyDict_Next, PyDict_Keys, PyDict_Update
from cpython.int cimport PyInt_AsLong, PyInt_Check
from cpython.string cimport PyString_Format
from cpython.float cimport PyFloat_AsDouble
from cpython.tuple cimport PyTuple_GetItem


# Forward Declaration
cdef: 
    dict std_mol_comp = {}
    str _atom = r'([A-Z][a-z+]*)(?:\[(\d+)\])?([+-]?\d+)?'
    str _formula = r'^({})*$'.format(_atom)
    str _isotope_string = r'^([A-Z][a-z+]*)(?:\[(\d+)\])?$'

    object isotope_pattern = re.compile(_isotope_string)
    object formula_pattern = re.compile(_formula)


@cython.boundscheck(False)
cdef inline str _parse_isotope_string(str label, int* isotope_num):
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


cdef inline str _make_isotope_string(str element_name, int isotope_num):
    """Form a string label for an isotope."""
    cdef:
        tuple parts
    if isotope_num == 0:
        return element_name
    else:
        parts = (element_name, isotope_num)
        return <str>PyString_Format('%s[%d]', parts)


cdef class CComposition(dict):
    '''Represent arbitrary elemental compositions'''
    def __str__(self):   # pragma: no cover
        return 'Composition({})'.format(dict.__repr__(self))

    def __repr__(self):  # pragma: no cover
        return str(self)

    def __iadd__(CComposition self, other):
        cdef:
            str elem
            long cnt
        for elem, cnt in other.items():
            self.setitem(elem, self.getitem(elem) + cnt)
        return self

    def __add__(self, other):
        cdef:
            str elem
            long cnt
            CComposition result
        if not isinstance(self, CComposition):
            other, self = self, other
        result = CComposition(self)
        for elem, cnt in other.items():
            result.setitem(elem, result.getitem(elem) + cnt)
        return result


    def __isub__(self, other):
        cdef:
            str elem
            long cnt
        for elem, cnt in other.items():
            self.setitem(elem, self.getitem(elem) - cnt)
        return self

    def __sub__(self, other):
        cdef:
            str elem
            long cnt
            CComposition result
        if not isinstance(self, CComposition):
            self = CComposition(self)
        result = CComposition(self)
        for elem, cnt in other.items():
            result.setitem(elem, result.getitem(elem) - cnt)
        return result

    def __reduce__(self):
        return composition_factory, (list(self),), self.__getstate__()

    def __getstate__(self):
        return dict(self)

    def __setstate__(self, d):
        self._from_dict(d)


    def __mul__(self, other):
        cdef:
            CComposition prod = CComposition()
            int rep, v
            str k

        if isinstance(other, CComposition):
            self, other = other, self
        
        if not isinstance(other, int):
            raise ChemicalCompositionError(
                'Cannot multiply Composition by non-integer',
                other)
        rep = other
        for k, v in self.items():
            prod.setitem(k, v * rep)
        return prod


    def __richcmp__(self, other, int code):
        if code == 2:
            if not isinstance(other, dict):
                return False
            self_items = set([i for i in self.items() if i[1]])
            other_items = set([i for i in other.items() if i[1]])
            return self_items == other_items
        else:
            return NotImplemented

    def __neg__(self):
        return self * -1

    # Override the default behavior, if a key is not present
    # do not initialize it to 0.
    def __missing__(self, str key):
        return 0

    def __setitem__(self, str key, int value):
        if value:  # Will not occur on 0 as 0 is falsey AND an integer
            super(Composition, self).__setitem__(key, value)
        elif key in self:
            del self[key]

    def copy(self):
        return CComposition(self)

    cdef inline long getitem(self, str elem):
        cdef:
            PyObject* resobj
            long count
        resobj = PyDict_GetItem(self, elem)
        if (resobj == NULL):
            return 0
        count = PyInt_AsLong(<object>resobj)
        return count

    cdef inline void setitem(self, str elem, long val):
        PyDict_SetItem(self, elem, val)

    cpdef CComposition clone(self):
        return CComposition(self)

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
        prev_chem_symbol_start = len(formula)
        i = len(formula) - 1

        seek_mode = 0
        parse_stack = ""
        resolve_stack = []
        group_coef = None

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
        # self.update(comp)
        PyDict_Update(self, comp)


    cpdef double calc_mass(self, int average=False, charge=None, dict mass_data=nist_mass) except -1:
        return calculate_mass(self, average=average, charge=charge, mass_data=mass_data)

    property mass:
        def __get__(self):
            return self.calc_mass()

    def __init__(self, *args, **kwargs):
        """
        A Composition object stores a chemical composition of a
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
        mol_comp : dict, optional
            A dict with the elemental composition of the standard molecules (the
            default value is std_mol_comp).
        mass_data : dict, optional
            A dict with the masses of chemical elements (the default
            value is :py:data:`nist_mass`). It is used for formulae parsing only.
        """
        dict.__init__(self)
        cdef:
            dict mol_comp, mass_data
            str kwa
            set kw_sources
        mol_comp = kwargs.get('mol_comp', std_mol_comp)
        mass_data = kwargs.get('mass_data', nist_mass)

        kw_sources = set(
            ('formula',))
        kw_given = kw_sources.intersection(kwargs)
        if len(kw_given) > 1:
            raise ChemicalCompositionError('Only one of {} can be specified!\n\
                Given: {}'.format(', '.join(kw_sources),
                                  ', '.join(kw_given)))

        elif kw_given:
            kwa = kw_given.pop()
            if kwa == "formula":
                self._from_formula(kwargs[kwa], mass_data)
        # can't build from kwargs
        elif args:
            if isinstance(args[0], dict):
                self._from_dict(args[0])
            elif isinstance(args[0], str):
                try:
                    self._from_formula(args[0], mass_data)
                except ChemicalCompositionError:
                    raise ChemicalCompositionError(
                        'Could not create a Composition object from '
                        'string: "{}": not a valid sequence or '
                        'formula'.format(args[0]))
        else:
            self._from_dict(kwargs)


Composition = CComposition

std_mol_comp.update({
    # Amino Acids
    'A':   Composition({'H': 5, 'C': 3, 'O': 1, 'N': 1}),
    'C':   Composition({'H': 5, 'C': 3, 'S': 1, 'O': 1, 'N': 1}),
    'D':   Composition({'H': 5, 'C': 4, 'O': 3, 'N': 1}),
    'E':   Composition({'H': 7, 'C': 5, 'O': 3, 'N': 1}),
    'F':   Composition({'H': 9, 'C': 9, 'O': 1, 'N': 1}),
    'G':   Composition({'H': 3, 'C': 2, 'O': 1, 'N': 1}),
    'H':   Composition({'H': 7, 'C': 6, 'N': 3, 'O': 1}),
    'I':   Composition({'H': 11, 'C': 6, 'O': 1, 'N': 1}),
    'K':   Composition({'H': 12, 'C': 6, 'N': 2, 'O': 1}),
    'L':   Composition({'H': 11, 'C': 6, 'O': 1, 'N': 1}),
    'M':   Composition({'H': 9, 'C': 5, 'S': 1, 'O': 1, 'N': 1}),
    'N':   Composition({'H': 6, 'C': 4, 'O': 2, 'N': 2}),
    'P':   Composition({'H': 7, 'C': 5, 'O': 1, 'N': 1}),
    'Q':   Composition({'H': 8, 'C': 5, 'O': 2, 'N': 2}),
    'R':   Composition({'H': 12, 'C': 6, 'N': 4, 'O': 1}),
    'S':   Composition({'H': 5, 'C': 3, 'O': 2, 'N': 1}),
    'T':   Composition({'H': 7, 'C': 4, 'O': 2, 'N': 1}),
    'V':   Composition({'H': 9, 'C': 5, 'O': 1, 'N': 1}),
    'W':   Composition({'C': 11, 'H': 10, 'N': 2, 'O': 1}),
    'Y':   Composition({'H': 9, 'C': 9, 'O': 2, 'N': 1}),

    # Protein Sequence Terminals (modX format)
    'H-':  Composition({'H': 1}),
    '-OH': Composition({'O': 1, 'H': 1}),

    # Glycans
    'Hex':    Composition({'H': 12, 'C': 6, 'O': 6}),
    'Pen':    Composition({'H': 10, 'C': 5, 'O': 5}),
    'HexNAc': Composition({'H': 13, 'C': 8, 'O': 5, 'N': 1}),
    'NeuAc':  Composition({'H': 17, 'C': 11, 'O': 8, 'N': 1}),
    'NeuGc':  Composition({'H': 17, 'C': 11, 'O': 9, 'N': 1}),
})


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef inline double calculate_mass(CComposition composition=None, str formula=None, int average=False, charge=None, mass_data=None) except -1:
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
        PyObject* interm

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
    for isotope_string in composition:
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
        mass /= _charge

    composition.setitem('H+', old_charge)
    return mass
