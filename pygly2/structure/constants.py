from ..utils.enum import EnumMeta


class SuperClass(object):
    __metaclass__ = EnumMeta
    tri = 3
    tet = 4
    pen = 5
    hex = 6
    hep = 7
    oct = 8
    non = 9
    missing = None
    x = missing


class Stem(object):
    __metaclass__ = EnumMeta
    gro = 1
    ery = 2
    rib = 3
    ara = 4
    all = 5
    alt = 6
    glc = 7
    man = 8
    tre = 9
    xyl = 10
    lyx = 11
    gul = 12
    ido = 13
    gal = 14
    tal = 15
    thr = 16
    missing = None
    x = missing


class Configuration(object):
    '''Optical Stereomers'''
    __metaclass__ = EnumMeta
    d = 1
    l = 2
    missing = None
    x = missing


class Modification(object):
    __metaclass__ = EnumMeta
    d = 1
    keto = 2
    en = 3
    a = 4
    aldi = 5
    sp2 = 6
    sp = 7
    geminal = 8
    _reserve = 9
    _cleave = 10
    _reduce = 11


ReducingEnd = Modification._reduce


class Anomer(object):
    __metaclass__ = EnumMeta
    alpha = 1
    beta = 2
    uncyclized = 3
    missing = None
    x = missing
