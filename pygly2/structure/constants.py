from ..utils.enum import Enum


class SuperClass(Enum):
    '''Corresponds to the number of carbon atoms in the carbohydrate backbone
    of the |Monosaccharide|

    Is an |Enum|
    '''
    sug = 2
    tri = 3
    tet = 4
    pen = 5
    hex = 6
    hep = 7
    oct = 8
    non = 9
    dec = 10
    s11 = 11
    s12 = 12
    s13 = 13
    s14 = 14
    s15 = 15
    s16 = 16
    s17 = 17
    s18 = 18
    s19 = 19
    s20 = 20
    x = None

SuperClass.sug.add_name("Sugose")
SuperClass.tri.add_name("Triose")
SuperClass.tet.add_name("Tetrose")
SuperClass.pen.add_name("Pentose")
SuperClass.hex.add_name("Hexose")
SuperClass.hep.add_name("Heparin")
SuperClass.oct.add_name("Octose")
SuperClass.non.add_name("Nonose")
SuperClass.dec.add_name("Decose")


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

#: An alias for :attr:`Modification.aldi`
Modification.aldi.add_name("Alditol")
Modification.a.add_name("Acidic")
Modification.d.add_name("Deoxygenated")
Modification.keto.add_name("Ketone")
Modification.en.add_name("DoubleBond")
Modification.geminal.add_name("Geminal")
Modification.sp.add_name("SP")
Modification.sp2.add_name("SP2")
ReducingEnd = Modification.aldi


class Anomer(Enum):
    '''
    Corresponds to the type of linkage found at the anomeric carbon of this |Monosaccharide|

    Is an |Enum|
    '''
    alpha = 1
    beta = 2
    uncyclized = 3
    x = None


class RingType(Enum):
    '''
    Corresponds to the type of ring structure of this |Monosaccharide|. Pyranose rings are
    five-member rings including one Oxygen and Furanose rings are four member rings including
    one Oxygen.

    Is an |Enum|
    '''
    pyranose = 6
    furanose = 5
    open = 0
    x = None
