from ..utils.enum import Enum


class SuperClass(Enum):
    '''Corresponds to the number of carbon atoms in the carbohydrate backbone
    of the |Monosaccharide|

    Is an |Enum|
    '''
    #: Disaccharide - 2 monosaccharides
    sug = 2
    #: Triose - 3 monosaccharides
    tri = 3
    #: Tetrose - 4 monosaccharides
    tet = 4
    #: Pentose - 5 monosaccharides
    pen = 5
    #: Hexose - 6 monosaccharides
    hex = 6
    #: Heptose - 7 monosaccharides
    hep = 7
    #: Octose - 8 monosaccharides
    oct = 8
    #: Nonose - 9 monosaccharides
    non = 9
    #: Decose - 10 monosaccharides
    dec = 10
    #: 11 monosaccharides
    s11 = 11
    #: 12 monosaccharides
    s12 = 12
    #: 13 monosaccharides
    s13 = 13
    #: 14 monosaccharides
    s14 = 14
    #: 15 monosaccharides
    s15 = 15
    #: 16 monosaccharides
    s16 = 16
    #: 17 monosaccharides
    s17 = 17
    #: 18 monosaccharides
    s18 = 18
    #: 19 monosaccharides
    s19 = 19
    #: 20 monosaccharides
    s20 = 20
    #: Unknown
    x = None

#: Sugose alias of `sug`
SuperClass.sug.add_name("Sugose")
#: Triose alias of `tri`
SuperClass.tri.add_name("Triose")
#: Tetrose alias of `tet`
SuperClass.tet.add_name("Tetrose")
#: Pentose alias of `pen`
SuperClass.pen.add_name("Pentose")
#: Hexose alias of `hex`
SuperClass.hex.add_name("Hexose")
#: Heparin alias of `hep`
SuperClass.hep.add_name("Heparin")
#: Octose alias of `oct`
SuperClass.oct.add_name("Octose")
#: Nonose alias of `non`
SuperClass.non.add_name("Nonose")
#: Decose alias of `dec`
SuperClass.dec.add_name("Decose")


class Stem(Enum):
    '''Corresponds to the bond formation pattern between the carbon atoms in the
    carbohydrate backbone of the |Monosaccharide|

    Is an |Enum|
    '''
    #: Glyceraldehyde
    gro = 1
    #: Erythrose
    ery = 2
    #: Ribose
    rib = 3
    #: Arabinose
    ara = 4
    #: Allose
    all = 5
    #: Altrose
    alt = 6
    #: Glucose
    glc = 7
    #: Mannose
    man = 8
    #: Threose
    tre = 9
    #: Xylose
    xyl = 10
    #: Lyxose
    lyx = 11
    #: Gulose
    gul = 12
    #: Idose
    ido = 13
    #: Galactose
    gal = 14
    #: Talose
    tal = 15
    #: Unknown
    x = None

    thr = 16


class Configuration(Enum):
    '''
    Corresponds to the optical stereomeric state of the |Monosaccharide|

    Is an |Enum|
    '''
    #: D Configuration
    d = 1
    #: L Configuration
    l = 2
    x = None


class Modification(Enum):
    '''
    Corresponds to discrete composition shifts of the |Monosaccharide| which
    are simple enough to not constitute a distinct object to represent like |Substituent|.

    Is an |Enum|
    '''
    #: Alditol
    d = 1
    #: Acidic
    keto = 2
    #: Deoxygenated
    en = 3
    #: Ketone
    a = 4
    #: DoubleBond
    aldi = 5
    #: Geminal
    sp2 = 6
    #: SP
    sp = 7
    #: SP2
    geminal = 8

#: alias of `aldi`
Modification.aldi.add_name("Alditol")
#: alias of `a`
Modification.a.add_name("Acidic")
#: alias of `d`
Modification.d.add_name("Deoxygenated")
#: alias of `keto`
Modification.keto.add_name("Ketone")
#: alias of `en`
Modification.en.add_name("DoubleBond")
#: alias of `geminal`
Modification.geminal.add_name("Geminal")
#: alias of `sp`
Modification.sp.add_name("SP")
#: alias of `sp2`
Modification.sp2.add_name("SP2")


class Anomer(Enum):
    '''
    Corresponds to the type of linkage found at the anomeric carbon of this |Monosaccharide|

    Is an |Enum|
    '''
    #: Alpha linkage
    alpha = 1
    #: Beta linkage
    beta = 2
    #: Uncyclized open chain
    uncyclized = 3
    #: Unknown
    x = None


class RingType(Enum):
    '''
    Corresponds to the type of ring structure of this |Monosaccharide|. Pyranose rings are
    five-member rings including one Oxygen and Furanose rings are four member rings including
    one Oxygen.

    Is an |Enum|
    '''
    #: Six member ring
    pyranose = 6
    #: Five member ring
    furanose = 5
    #: Open chain
    open = 0
    #: Unknown
    x = None
