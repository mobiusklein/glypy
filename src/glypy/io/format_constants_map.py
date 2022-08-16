from glypy.composition import Composition
from glypy.structure import constants

anomer_map = {
    'a': constants.Anomer.alpha,
    'b': constants.Anomer.beta,
    'o': constants.Anomer.uncyclized,
    'x': constants.Anomer.x
}


SuperClass = constants.SuperClass
superclass_map = {
    "sug".upper(): SuperClass["sug"],
    "tri".upper(): SuperClass["tri"],
    "tet".upper(): SuperClass["tet"],
    "pen".upper(): SuperClass["pen"],
    "hex".upper(): SuperClass["hex"],
    "hep".upper(): SuperClass["hep"],
    "oct".upper(): SuperClass["oct"],
    "non".upper(): SuperClass["non"],
    "dec".upper(): SuperClass["dec"],
    "s11".upper(): SuperClass["s11"],
    "s12".upper(): SuperClass["s12"],
    "s13".upper(): SuperClass["s13"],
    "s14".upper(): SuperClass["s14"],
    "s15".upper(): SuperClass["s15"],
    "s16".upper(): SuperClass["s16"],
    "s17".upper(): SuperClass["s17"],
    "s18".upper(): SuperClass["s18"],
    "s19".upper(): SuperClass["s19"],
    "s20".upper(): SuperClass["s20"],
    "x".upper(): SuperClass["x"],
}


link_replacement_composition_map = {
    "o": Composition(H=1),  # The source Carbon's Hydroxyl group loses its Hydrogen
    "d": Composition("OH"),  # The Hydroxyl group of the target Carbon is lost
    "h": Composition(H=1),  # The Hydrogen of the target Carbon is lost
    "n": Composition(H=1),  # Non-sugar unit recieves bond
    "x": Composition({})  # Unknown
}


linkage_type_map = {
    'o': constants.LinkageType.backbone_oxygen,
    'd': constants.LinkageType.backbone_oxygen,
    'h': constants.LinkageType.backbone_hydrogen,
    'n': constants.LinkageType.other,
    'x': constants.LinkageType.x
}


Modification = constants.Modification
modification_map = {
    'd': Modification.Deoxygenated,
    'a': Modification.Acidic,
    'aldi': Modification.Alditol,
    'keto': Modification.Ketone,
    'en': Modification.DoubleBond,
    'sp': Modification.SP,
    'sp2': Modification.SP2,
    'geminal': Modification.Geminal
}
