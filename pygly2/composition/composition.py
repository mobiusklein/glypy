# Credit to Pyteomics - http://pythonhosted.org/pyteomics - for majority of design

from collections import defaultdict
from .mass_dict import nist_mass


class ChemicalCompositionError(Exception):
    pass


# Forward Declaration
std_mol_comp = {}


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
    if label.endswith(']'):
        isotope_num = int(label[label.find('[') + 1:-1])
        element_name = label[:label.find('[')]
    else:
        isotope_num = 0
        element_name = label
    return (element_name, isotope_num)


class Composition(defaultdict):
    '''Represent arbitrary elemental compositions'''
    def __str__(self):   # pragma: no cover
        return 'Composition({})'.format(dict.__repr__(self))

    def __repr__(self):  # pragma: no cover
        return str(self)

    def __add__(self, other):
        result = self.copy()
        for elem, cnt in other.items():
            result[elem] += cnt
        return result

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        result = self.copy()
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

    # Override the default behavior, if a key is not present
    # do not initialize it to 0.
    def __missing__(self, key):
        return 0

    def __setitem__(self, key, value):
        if isinstance(value, float):
            value = int(round(value))
        elif not isinstance(value, int):
            raise ChemicalCompositionError(
                'Only integers allowed as values in \
                Composition, got {}.'.format(type(value).__name__))
        if value:  # Will not occur on 0 as 0 is falsey AND an integer
            super(Composition, self).__setitem__(key, value)
        elif key in self:
            del self[key]

    def copy(self):
        return Composition(self)

    # def _from_parsed_sequence(self, parsed_sequence, **kwargs):
    #     raise NotImplementedError()

    # def _from_split_sequence(self, split_sequence, **kwargs):
    #     raise NotImplementedError()

    # def _from_sequence(self, sequence, **kwargs):
    #     raise NotImplementedError()

    def _from_formula(self, formula, mass_data):
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
                        bra_pos = formula.rfind('[', 0, i)
                        if bra_pos == -1:
                            raise ChemicalCompositionError(
                                'Badly-formed isotope number: %s' % formula)
                        try:
                            isotope_num = int(formula[bra_pos + 1:i])
                        except ValueError:
                            raise ChemicalCompositionError(
                                'Badly-formed isotope number: %s' % formula)
                        i = bra_pos - 1
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

    def calc_mass(self, *args, **kwargs):
        kwargs["composition"] = self
        return calculate_mass(*args, **kwargs)

    @property
    def mass(self):
        return calculate_mass(composition=self)

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
        defaultdict.__init__(self, int)

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
            getattr(self, '_from_' + kwa)(kwargs[kwa],
                                          mass_data if kwa == 'formula' else mol_comp)
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

std_ion_comp = {
    'M':        Composition(formula=''),
    'a':        Composition(formula='H-2O-1' + 'C-1O-1'),
    'a-H2O':    Composition(formula='H-2O-1' + 'C-1O-1' + 'H-2O-1'),
    'a-NH3':    Composition(formula='H-2O-1' + 'C-1O-1' + 'N-1H-3'),
    'b':        Composition(formula='H-2O-1'),
    'b-H2O':    Composition(formula='H-2O-1' + 'H-2O-1'),
    'b-NH3':    Composition(formula='H-2O-1' + 'N-1H-3'),
    'c':        Composition(formula='H-2O-1' + 'NH3'),
    'c-H2O':    Composition(formula='H-2O-1' + 'NH3' + 'H-2O-1'),
    'c-NH3':    Composition(formula='H-2O-1'),
    'x':        Composition(formula='H-2O-1' + 'CO2'),
    'x-H2O':    Composition(formula='H-2O-1' + 'CO2' + 'H-2O-1'),
    'x-NH3':    Composition(formula='H-2O-1' + 'CO2' + 'N-1H-3'),
    'y':        Composition(formula=''),
    'y-H2O':    Composition(formula='H-2O-1'),
    'y-NH3':    Composition(formula='N-1H-3'),
    'z':        Composition(formula='H-2O-1' + 'ON-1H-1'),
    'z-H2O':    Composition(formula='H-2O-1' + 'ON-1H-1' + 'H-2O-1'),
    'z-NH3':    Composition(formula='H-2O-1' + 'ON-1H-1' + 'N-1H-3'),
}


def calculate_mass(*args, **kwargs):
    """Calculates the monoisotopic mass of a polypeptide defined by a
    sequence string, parsed sequence, chemical formula or
    Composition object.

    Parameters
    ----------
    formula : str, optional
        A string with a chemical formula.
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
    mol_comp : dict, optional
        A dict with the elemental composition of the commonly used molecules (the
        default value is std_mol_comp).
    mass_data : dict, optional
        A dict with the masses of the chemical elements (the default
        value is :py:data:`nist_mass`).
    ion_comp : dict, optional
        A dict with the relative elemental compositions of peptide ion
        fragments (default is :py:data:`std_ion_comp`).

    Returns
    -------
        mass : float
    """

    mass_data = kwargs.get('mass_data', nist_mass) or nist_mass
    ion_comp = kwargs.get('ion_comp', std_ion_comp)
    # Make a deep copy of `composition` keyword argument.
    composition = (Composition(kwargs['composition'])
                   if 'composition' in kwargs
                   else Composition(*args, **kwargs))

    if 'ion_type' in kwargs:
        composition += ion_comp[kwargs['ion_type']]

    # Get charge.
    charge = composition['H+']
    if 'charge' in kwargs:
        if charge:
            raise ChemicalCompositionError(
                'Charge is specified both by the number of protons and '
                '`charge` in kwargs')
        charge = kwargs['charge']
        composition['H+'] = charge

    # Calculate mass.
    mass = 0.0
    average = kwargs.get('average', False)
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
        mass /= charge
    return mass
