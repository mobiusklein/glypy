# Credit to Pyteomics - http://pythonhosted.org/pyteomics - for majority of design
import re
import math
from collections import defaultdict
from .mass_dict import nist_mass
from .base import ChemicalCompositionError, composition_factory

import os

use_cython = bool(int(os.environ.get("USE_CYTHON_COMPOSITION", 1)))

try:
    from .ccomposition import CComposition, calculate_mass as ccalculate_mass
except ImportError:  # pragma: no cover
    use_cython = False

# Forward Declaration

_isotope_string = r'^([A-Z][a-z+]*)(?:\[(\d+)\])?$'
_atom = r'([A-Z][a-z+]*)(?:\[(\d+)\])?([+-]?\d+)?'
_formula = r'^({})*$'.format(_atom)
formula_pattern = re.compile(_formula)


def _make_isotope_string(element_name, isotope_num):
    """Form a string label for an isotope."""
    if isotope_num == 0:
        return element_name
    else:
        return '%s[%d]' % (element_name, isotope_num)


def _parse_isotope_string(label):
    """Parse an string with an isotope label and return the element name and
    the isotope number.

    >>> _parse_isotope_string('C')
    ('C', 0)
    >>> _parse_isotope_string('C[12]')
    ('C', 12)
    """
    element_name, num = re.match(_isotope_string, label).groups()
    isotope_num = int(num) if num else 0
    return element_name, isotope_num


class PComposition(defaultdict):
    '''Represent arbitrary elemental compositions'''
    def __str__(self):   # pragma: no cover
        return 'Composition({})'.format(dict.__repr__(self))

    def __repr__(self):  # pragma: no cover
        return str(self)

    def __iadd__(self, other):
        for elem, cnt in (other.items()):
            self[elem] += cnt
        self._mass_args = None
        return self

    def __add__(self, other):
        result = self.clone()
        for elem, cnt in other.items():
            result[elem] += cnt
        return result

    def __radd__(self, other):
        return self + other

    def __isub__(self, other):
        for elem, cnt in other.items():
            self[elem] -= cnt
        self._mass_args = None
        return self

    def __sub__(self, other):
        result = self.clone()
        for elem, cnt in other.items():
            result[elem] -= cnt
        return result

    def __rsub__(self, other):
        return (self - other) * (-1)

    def __mul__(self, other):
        if not isinstance(other, int):
            raise ChemicalCompositionError(
                'Cannot multiply Composition by non-integer',
                other)

        prod = {}
        for k, v in self.items():
            prod[k] = v * other

        return Composition(prod)

    def __rmul__(self, other):
        return self * other

    def __eq__(self, other):
        if not isinstance(other, dict):
            return False
        self_items = set([i for i in self.items() if i[1]])
        other_items = set([i for i in other.items() if i[1]])
        return self_items == other_items

    def __neg__(self):
        return -1 * self

    def __reduce__(self):
        return composition_factory, (list(self),), self.__getstate__()

    def __getstate__(self):
        return dict(self)

    def __setstate__(self, d):
        self._from_dict(d)

    # Override the default behavior, if a key is not present
    # do not initialize it to 0.
    def __missing__(self, key):
        return 0

    def __setitem__(self, key, value):
        if isinstance(value, float):
            value = int(round(value))
        elif not isinstance(value, (int)):
            raise ChemicalCompositionError(
                'Only integers allowed as values in Composition, got {}.'.format(type(value).__name__))
        if value:  # Will not occur on 0 as 0 is falsey AND an integer
            super(PComposition, self).__setitem__(key, value)
        elif key in self:
            del self[key]
        self._mass_args = None

    def clone(self):
        return PComposition(self)

    copy = clone

    def update(self, *args, **kwargs):
        dict.update(self, *args, **kwargs)
        self._mass_args = None

    def _from_formula(self, formula, mass_data):
        if '(' in formula:
            return self._from_formula_parens(formula, mass_data)
        if not formula_pattern.match(formula):
            raise ChemicalCompositionError('Invalid formula: ' + formula)
        for elem, isotope, number in re.findall(_atom, formula):
            if elem not in mass_data:
                raise ChemicalCompositionError('Unknown chemical element: ' + elem)
            self[_make_isotope_string(elem, int(isotope) if isotope else 0)
                 ] += int(number) if number else 1

    def _from_formula_parens(self, formula, mass_data):  # pragma: no cover
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

    def _from_dict(self, comp):
        # for isotope_string, num_atoms in comp.items():
        #     element_name, isotope_num = _parse_isotope_string(
        #         isotope_string)

        #     # Remove explicitly undefined isotopes (e.g. X[0]).
        #     self[_make_isotope_string(element_name, isotope_num)] = (
        #         num_atoms)
        self.update(comp)

    def calc_mass(self, average=False, charge=None, mass_data=None):
        mdid = id(mass_data)
        if self._mass_args is not None and average is self._mass_args[0]\
                and charge == self._mass_args[1] and mdid == self._mass_args[2]:
            return self._mass
        else:
            self._mass_args = (average, charge, mdid)
            self._mass = pcalculate_mass(composition=self, average=average, charge=charge, mass_data=mass_data)
            return self._mass

    @property
    def mass(self):
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
        mass_data : dict, optional
            A dict with the masses of chemical elements (the default
            value is :py:data:`nist_mass`). It is used for formulae parsing only.
        """
        defaultdict.__init__(self, int)
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
            getattr(self, '_from_' + kwa)(kwargs[kwa], mass_data)
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

        self._mass = None
        self._mass_args = None


def pcalculate_mass(composition=None, formula=None, average=False, charge=None, mass_data=None):
    """Calculates the monoisotopic mass of a chemical formula or
    Composition object.

    Parameters
    ----------
    composition : Composition, optional
        A Composition object with the elemental composition of a substance.
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

    mass_data = mass_data or nist_mass
    # Make a deep copy of `composition` keyword argument.
    composition = composition if composition else PComposition(formula)

    # Get charge.
    old_charge = composition['H+']
    if charge is not None:
        if old_charge != 0:
            raise ChemicalCompositionError(
                'Charge is specified both by the number of protons and '
                '`charge` in kwargs')
        composition['H+'] = charge

    # Calculate mass.
    mass = 0.0
    for isotope_string in composition:
        element_name, isotope_num = _parse_isotope_string(isotope_string)
        # Calculate average mass if required and the isotope number is
        # not specified.
        if (not isotope_num) and average:
            for isotope in mass_data[element_name]:
                if isotope != 0:
                    mass += (composition[element_name]
                             * mass_data[element_name][isotope][0]
                             * mass_data[element_name][isotope][1])
        else:
            mass += (composition[isotope_string]
                     * mass_data[element_name][isotope_num][0])

    # Calculate m/z if required.
    if charge:
        mass /= abs(charge)
    if old_charge != 0:
        composition["H+"] = old_charge
    elif charge:
        del composition["H+"]
    return mass


# Checks to see if the Cython versions are available
if use_cython:  # pragma: no cover
    Composition = CComposition
    calculate_mass = ccalculate_mass
else:  # pragma: no cover
    Composition = PComposition
    calculate_mass = pcalculate_mass


def most_probable_isotopic_composition(*args, **kwargs):
    """Calculate the most probable isotopic composition of a
    chemical formula or |Composition| object.

    For each element, only two most abundant isotopes are considered.

    Parameters
    ----------
    formula : str, optional
        A string with a chemical formula.
    composition : |Composition|, optional
        A :py:class:`Composition` object with the elemental composition of a
        substance.
    elements_with_isotopes : list of str
        A list of elements to be considered in isotopic distribution
        (by default, every element has a isotopic distribution).
    mass_data : dict, optional
        A dict with the masses of chemical elements (.

    Returns
    -------
    out: tuple (|Composition|, float)
        A tuple with the most probable isotopic composition and its
        relative abundance.
    """

    composition = (dict(kwargs['composition']) if 'composition' in kwargs
                   else Composition(*args, **kwargs))

    # Removing isotopes from the composition.
    for isotope_string in composition:
        element_name, isotope_num = _parse_isotope_string(isotope_string)
        if isotope_num:
            composition[element_name] += composition.pop(isotope_string)

    mass_data = kwargs.get('mass_data', nist_mass)
    elements_with_isotopes = kwargs.get('elements_with_isotopes')
    isotopic_composition = Composition()

    for element_name in composition:
        if element_name == 'H+':
            continue
        if (not elements_with_isotopes or (element_name in elements_with_isotopes)):
            # Take the two most abundant isotopes.
            first_iso, second_iso = sorted(
                [(i[0], i[1][1]) for i in mass_data[element_name].items() if i[0]],
                key=lambda x: -x[1])[:2]

            # Write the number of isotopes of the most abundant type.
            first_iso_str = _make_isotope_string(element_name, first_iso[0])
            isotopic_composition[first_iso_str] = round(int(math.ceil(composition[element_name])) * first_iso[1])

            # Write the number of the second isotopes.
            second_iso_str = _make_isotope_string(element_name, second_iso[0])
            isotopic_composition[second_iso_str] = round(
                (composition[element_name] - isotopic_composition[first_iso_str]))
        else:
            isotopic_composition[element_name] = composition[element_name]

    return (isotopic_composition,
            isotopic_composition_abundance(
                composition=isotopic_composition,
                mass_data=mass_data))


def isotopic_composition_abundance(*args, **kwargs):
    """Calculate the relative abundance of a given isotopic composition
    of a molecule.

    Parameters
    ----------
    formula : str, optional
        A string with a chemical formula.
    composition : Composition, optional
        A Composition object with the isotopic composition of a substance.
    mass_data : dict, optional
        A dict with the masses of chemical elements (the default
        value is :py:data:`nist_mass`).

    Returns
    -------
    relative_abundance : float
        The relative abundance of a given isotopic composition.
    """

    composition = (Composition(kwargs['composition'])
                   if 'composition' in kwargs
                   else Composition(*args, **kwargs))

    isotopic_composition = defaultdict(dict)

    # Check if there are default and non-default isotopes of the same
    # element and rearrange the elements.
    for element in composition:
        if element == 'H+':
            continue
        element_name, isotope_num = _parse_isotope_string(element)

        # If there is already an entry for this element and either it
        # contains a default isotope or newly added isotope is default
        # then raise an exception.
        if ((element_name in isotopic_composition) and (isotope_num == 0 or 0 in isotopic_composition[element_name])):
            raise ChemicalCompositionError(
                'Please specify the isotopic states of all atoms of '
                '%s or do not specify them at all.' % element_name)
        else:
            isotopic_composition[element_name][isotope_num] = (
                composition[element])

    # Calculate relative abundance.
    mass_data = kwargs.get('mass_data', nist_mass)
    num1, num2, denom = 1., 1., 1.
    for element_name, isotope_dict in isotopic_composition.items():
        num1 *= math.factorial(sum(isotope_dict.values()))
        for isotope_num, isotope_content in isotope_dict.items():
            denom *= math.factorial(isotope_content)
            if isotope_num:
                num2 *= (mass_data[element_name][isotope_num][1] ** isotope_content)

    return num2 * (num1 / denom)


def formula(composition):
    keys = sorted(composition)
    return ''.join("%s%d" % (k, composition[k]) for k in keys)
