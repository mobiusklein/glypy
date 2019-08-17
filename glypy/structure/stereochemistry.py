'''
Defined by `Extended Stereocode <http://monosaccharidedb.org/notation.action?topic=stereocode#ext_stereocode>`_

Symbol  Description
h   "head or tail group", CH2OH group at a terminal position
d   DEOXY core modification at non-terminal position
m   DEOXY core modification at terminal position ("methyl" group)
a   ACID core modification
o   aldehyde group
k   KETO core modification at non-terminal position
e   EN + deoxy core modifications
n   EN core modification without DEOXY core modification
E   EN core modification with unknown deoxygenation status
y   YN core modification at non-terminal position
s   SP2 core modification
t   SP core modification (always at terminal position)
1   "L-Configuration" carbon atom
2   "D-Configuration" carbon atom
x   unknown configuration (D or L) carbon atom
'''

import warnings
from six import string_types as basestring
from six.moves import range

from glypy.utils.enum import EnumValue

from .constants import Stereocoding, Modification, Anomer, Configuration, Stem, UnknownPosition, NoPosition


def sdecode(code):
    return [Stereocoding[c] for c in code]


def invert(code):
    out = []
    L = 1
    D = 2
    for c in code:
        if c == L:
            out.append(D)
        elif c == D:
            out.append(L)
        else:
            out.append(c)
    return sdecode(out)


reference_stereomap = {
    ('alpha', 'l', (('man',), 'hex')): '122110',
    ('alpha', 'd', (('man',), 'hex')): '211220',
    ('beta', 'd', (('man',), 'hex')): '111220',
    ('beta', 'l', (('man',), 'hex')): '222110',
    ('uncyclized', 'd', (('gal',), 'hex')): '021120',
    ('alpha', 'l', (('glc',), 'hex')): '112110',
    ('beta', 'l', (('ara',), 'pen')): '22110',
    ('beta', 'l', (('alt',), 'hex')): '221110',
    ('alpha', 'l', (('gro', 'man'), 'hep')): '2112210',
    ('beta', 'd', (('gro', 'gal'), 'non')): '010211220',
    ('alpha', 'd', (('gal',), 'hex')): '221120',
    ('beta', 'd', (('gul',), 'hex')): '122120',
    ('alpha', 'd', (('ido',), 'hex')): '212120',
    ('alpha', 'd', (('rib',), 'hex')): '220220',
    ('beta', 'd', (('all',), 'hex')): '122220',
    ('alpha', 'd', (('all',), 'hex')): '222220',
    ('beta', 'd', (('xyl',), 'hex')): '121020',
    ('beta', 'd', (('gal',), 'hex')): '121120',
    ('alpha', 'd', (('gro', 'gal'), 'non')): '020211220',
    ('alpha', 'd', (('gro', 'alt'), 'hep')): '2122220',
    ('beta', 'd', (('thr',), 'hex')): '100120',
    ('alpha', 'l', (('ara',), 'pen')): '12110',
    ('alpha', 'l', (('lyx',), 'hex')): '102210',
    ('alpha', 'd', (('ara',), 'hex')): '021220',
    ('beta', 'l', (('glc',), 'hex')): '212110',
    ('beta', 'l', (('ara',), 'hex')): '202110',
    ('uncyclized', 'd', (('ery',), 'pen')): '00220',
    ('beta', 'd', (('rib',), 'hex')): '102220',
    ('uncyclized', 'd', (('glc',), 'hex')): '021220',
    ('alpha', 'l', (('tal',), 'hex')): '122210',
    ('alpha', 'd', (('glc',), 'hex')): '221220',
    ('alpha', 'd', (('gro', 'man'), 'hep')): '2112220',
    ('alpha', 'd', (('rib',), 'pen')): '22220',
    ('beta', 'd', (('lyx',), 'hex')): '101120',
    ('alpha', 'd', (('man',), 'oct')): '02011220',
    ('beta', 'l', (('gul',), 'hex')): '211210',
    ('beta', 'l', (('rib',), 'hex')): '201110',
    ('alpha', 'd', (('xyl',), 'pen')): '22120',
    ('alpha', 'l', (('gal',), 'hex')): '112210',
    ('alpha', 'l', (('ido',), 'hex')): '121210',
    ('beta', 'l', (('ido',), 'hex')): '221210',
    ('beta', 'l', (('gal',), 'hex')): '212210',
    ('uncyclized', 'l', (('ara',), 'pen')): '02110',
    ('alpha', 'l', (('rib',), 'hex')): '101110',
    ('uncyclized', 'd', (('thr',), 'pen')): '00120',
    ('beta', 'l', (('thr',), 'hex')): '221000',
    ('beta', 'd', (('xyl',), 'pen')): '12120',
    ('alpha', 'l', (('gul',), 'hex')): '111210',
    ('beta', 'd', (('glc',), 'hex')): '121220',
    ('beta', 'l', (('tal',), 'hex')): '222210',
    ('alpha', 'l', (('ara',), 'hex')): '102110',
    ('beta', 'd', (('ara',), 'hex')): '011220',
    ('beta', 'd', (('man',), 'oct')): '01011220',
    ('alpha', 'd', (('man',), 'oct')): '02011220',
}

reference_stereomap = {k: sdecode(v) for k, v in reference_stereomap.items()}


def get_stereocode_key(monosaccharide, anomer=None, configuration=None):
    if anomer is None:
        anomer = monosaccharide.anomer
    if configuration is None:
        configuration = monosaccharide.configuration[0]
    key = (anomer, configuration,
           (monosaccharide.stem, monosaccharide.superclass))
    return key


def _update_stereocode_basic(code, monosaccharide):
    for i in range(len(code)):
        modifications = monosaccharide.modifications[i + 1]
        for mod in modifications:
            if mod.name == Modification.d:
                code[i] = 0
            elif mod.name == Modification.a.name:
                code[i] = 0
    return code


def _update_stereocode_extended(code, monosaccharide):
    null_position = (UnknownPosition, NoPosition)
    if monosaccharide.ring_start not in null_position and monosaccharide.ring_end not in null_position:
        ring_range = range(monosaccharide.ring_start, monosaccharide.ring_end + 1)
    else:
        ring_range = []
    ring_dimensions = range(1, monosaccharide.superclass.value + 1)
    termini = set(ring_dimensions) - set(ring_range)
    terminal = max(termini) - 1
    code[terminal] = Stereocoding.h
    for i in range(len(code)):
        modifications = monosaccharide.modifications[i + 1]
        for mod in modifications:
            if mod.name == Modification.d.name:
                if (i + 1) in ring_range:
                    code[i] = Stereocoding.d
                else:
                    code[i] = Stereocoding.m
            elif mod.name == Modification.a:
                code[i] = Stereocoding.a
    return code


class Stereocode(object):
    def __init__(self, code):
        if isinstance(code, basestring):
            code = sdecode(code)
        self.code = list(code)

    def __getitem__(self, i):
        return self.code[i]

    def __setitem__(self, i, v):
        if not isinstance(v, EnumValue):
            self.code[i] = Stereocoding[str(v)]
        else:
            self.code[i] = Stereocoding[v.name]

    def __iter__(self):
        return iter(self.code)

    def __len__(self):
        return len(self.code)

    def basic_encoding(self):
        encoding = {
            Stereocoding.x: 'x',
            Stereocoding.h: '0',
            Stereocoding.a: '0',
            Stereocoding.d: '0'
        }

        return ''.join(encoding.get(c, c.value) for c in self.code)

    def extended_encoding(self):
        return ''.join(c.value if c.value is not None else 'x' for c in self.code)

    def __str__(self):
        return self.extended_encoding()

    def __repr__(self):
        return "Stereocode(%r)" % (self.extended_encoding(),)


def _resolve_unknown_anomer(monosaccharide):
    configuration = monosaccharide.configuration[0]
    has_known_configuration = (configuration != Configuration.x)
    if not has_known_configuration:
        configuration = Configuration.d
    a_key = get_stereocode_key(monosaccharide, Anomer.a, configuration)
    b_key = get_stereocode_key(monosaccharide, Anomer.b, configuration)
    a_code = reference_stereomap[a_key]
    b_code = reference_stereomap[b_key]
    fill_value = Stereocoding[None]
    if monosaccharide.anomer == Anomer.uncyclized:
        fill_value = Stereocoding.h
    substituted_code = [ai if ai == bi else fill_value for ai, bi in zip(a_code, b_code)]
    if not has_known_configuration:
        translate_unknown_configuration = {
            Stereocoding.L: Stereocoding.LD,
            Stereocoding.D: Stereocoding.DL,

        }
        substituted_code = [
            translate_unknown_configuration.get(i, i) for i in substituted_code
        ]
    return substituted_code


def stereocode(monosaccharide):
    key = get_stereocode_key(monosaccharide)
    try:
        if monosaccharide.stem[0] == Stem.x:
            base_code = Stereocode(Stereocoding[None] for i in range(monosaccharide.superclass.value))
        else:
            try:
                base_code = Stereocode(reference_stereomap[key])
            except KeyError:
                base_code = _resolve_unknown_anomer(monosaccharide)
    except KeyError:
        warnings.warn("Could not locate a stereocode template for %r based upon key %s" % (
            monosaccharide, key), stacklevel=3)
        base_code = Stereocode(Stereocoding[None] for i in range(monosaccharide.superclass.value))
    _update_stereocode_extended(base_code, monosaccharide)
    return base_code
