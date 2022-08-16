from glypy.utils.enum import Enum


class SuperClass(Enum):
    '''Corresponds to the number of carbon atoms in the carbohydrate backbone
    of the |Monosaccharide|

    Is an |Enum|
    '''
    #: 2 carbons
    sug = 2
    #: Triose - 3 carbons
    tri = 3
    #: Tetrose - 4 carbons
    tet = 4
    #: Pentose - 5 carbons
    pen = 5
    #: Hexose - 6 carbons
    hex = 6
    #: Heptose - 7 carbons
    hep = 7
    #: Octose - 8 carbons
    oct = 8
    #: Nonose - 9 carbons
    non = 9
    #: Decose - 10 carbons
    dec = 10
    #: 11 carbons
    s11 = 11
    #: 12 carbons
    s12 = 12
    #: 13 carbons
    s13 = 13
    #: 14 carbons
    s14 = 14
    #: 15 carbons
    s15 = 15
    #: 16 carbons
    s16 = 16
    #: 17 carbons
    s17 = 17
    #: 18 carbons
    s18 = 18
    #: 19 carbons
    s19 = 19
    #: 20 carbons
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
#: Heptose alias of `hep`
SuperClass.hep.add_name("Heptose")
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


Stem.gro.add_name("Glyceraldehyde")
Stem.ery.add_name("Erythrose")
Stem.rib.add_name("Ribose")
Stem.ara.add_name("Arabinose")
Stem.all.add_name("Allose")
Stem.alt.add_name("Altose")
Stem.glc.add_name("Glucose")
Stem.man.add_name("Mannose")
Stem.tre.add_name("Threose")
Stem.xyl.add_name("Xylose")
Stem.lyx.add_name("Lyxose")
Stem.gul.add_name("Gulose")
Stem.ido.add_name("Idose")
Stem.gal.add_name("Galactose")
Stem.tal.add_name("Talose")
Stem.x.add_name("Unknown")


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


Configuration.d.add_name("Dextro")
Configuration.l.add_name("Levo")
Configuration.x.add_name("Unknown")


class Modification(Enum):
    '''
    Corresponds to discrete composition shifts of the |Monosaccharide| which
    are simple enough to not constitute a distinct object to represent like |Substituent|.

    Is an |Enum|
    '''
    #: Deoxygenated
    d = 1
    #: Ketone
    keto = 2
    #: DoubleBond
    en = 3
    #: Acidic
    a = 4
    #: Alditol
    aldi = 5
    #: SP2
    sp2 = 6
    #: SP
    sp = 7
    #: Geminal
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


Anomer.beta.add_name("b")
Anomer.alpha.add_name("a")
Anomer.uncyclized.add_name("o")
Anomer.uncyclized.add_name("open-chain")


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


class Stereocoding(Enum):
    x = None
    h = 'h'
    L = '1'
    D = '2'
    LD = '3'
    DL = '4'
    d = 'd'
    m = 'm'
    a = 'a'
    o = 'o'
    k = 'k'
    e = 'e'
    n = 'n'
    E = 'E'
    y = 'y'
    s = 's'
    t = 't'


Stereocoding.h.add_name('0')


UnknownPosition = -1
NoPosition = None


class LinkageType(Enum):
    backbone_oxygen = 0
    backbone_hydrogen = 1
    other = 2
    unknown = None


LinkageType.unknown.add_name('x')
