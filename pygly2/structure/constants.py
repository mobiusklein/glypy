from ..utils.enum import Enum


class SuperClass(Enum):
    '''Corresponds to the number of carbon atoms in the carbohydrate backbone
    of the |Monosaccharide|

    Is an |Enum|
    '''
    tri = 3
    tet = 4
    pen = 5
    hex = 6
    hep = 7
    oct = 8
    non = 9
    x = None


class Stem(Enum):
    '''Corresponds to the bond formation pattern between the carbon atoms in the
    carbohydrate backbone of the |Monosaccharide|

    Is an |Enum|
    '''
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
    x = None


class Configuration(Enum):
    '''
    Corresponds to the optical stereomeric state of the |Monosaccharide|

    Is an |Enum|
    '''
    d = 1
    l = 2
    x = None


class Modification(Enum):
    '''
    Corresponds to discrete composition shifts of the |Monosaccharide| which
    are simple enough to not constitute a distinct object to represent like |Substituent|.

    Is an |Enum|
    '''
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

#: An alias for :attr:`Modification._reduce`
ReducingEnd = Modification._reduce


class Anomer(Enum):
    '''
    Corresponds to the type of linkage found at the anomeric carbon of this |Monosaccharide|

    Is an |Enum|
    '''
    alpha = 1
    beta = 2
    uncyclized = 3
    x = None
