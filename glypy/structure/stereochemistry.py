import warnings
from six import string_types as basestring

from glypy.utils.enum import EnumValue

from .constants import Stereocoding, Modification


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
}

reference_stereomap = {k: sdecode(v) for k, v in reference_stereomap.items()}


def get_stereocode_key(monosaccharide):
    key = (monosaccharide.anomer, monosaccharide.configuration[0],
           (monosaccharide.stem, monosaccharide.superclass))

    return key


def _update_stereocode_basic(code, monosaccharide):
    for i in range(len(code)):
        site = monosaccharide[i + 1]
        modifications = site['modifications']
        for mod in modifications:
            print mod, i
            if mod.name == Modification.d:
                code[i] = 0
            elif mod.name == Modification.a.name:
                code[i] = 0
    return code


def _update_stereocode_extended(code, monosaccharide):
    ring_range = range(monosaccharide.ring_start, monosaccharide.ring_end + 1)
    ring_dimensions = range(1, monosaccharide.superclass.value + 1)
    termini = set(ring_dimensions) - set(ring_range)
    terminal = max(termini) - 1
    code[terminal] = Stereocoding.h
    for i in range(len(code)):
        site = monosaccharide[i + 1]
        modifications = site['modifications']
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


def stereocode(monosaccharide):
    key = get_stereocode_key(monosaccharide)
    try:
        base_code = Stereocode(reference_stereomap[key])
    except KeyError:
        warnings.warn("Could not locate a stereocode template for %r based upon key %s" % (
            monosaccharide, key))
        print "Could not locate a stereocode template for %r based upon key %s" % (
            monosaccharide, key)
        base_code = Stereocode(Stereocoding[None] for i in range(monosaccharide.superclass.value))
    _update_stereocode_extended(base_code, monosaccharide)
    return base_code
