from ..composition import Composition
from ..structure import constants

anomer_map = {
    'a': constants.Anomer.alpha,
    'b': constants.Anomer.beta,
    'o': constants.Anomer.uncyclized,
    'x': constants.Anomer.missing
}

superclass_map = {k.upper(): v for k, v in constants.SuperClass}


link_replacement_composition_map = {
    "o": Composition(H=1),  # The source Carbon's Hydroxyl group loses its Hydrogen
    "d": Composition(O=1, H=1),  # The Hydroxyl group of the target Carbon is lost
    "h": Composition(H=1),  # The Hydrogen of the target Carbon is lost
    "n": Composition(H=1),  # Non-sugar unit recieves bond
    "x": None  # Unknown
}