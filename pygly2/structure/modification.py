'''

'''

from .base import ModificationBase
from ..composition.structure_composition import modification_compositions
from ..io.format_constants_map import modification_map

glycoct_map = {v: k for k, v in modification_map.items()}


def resolve_composition_rule(name):
    try:
        return modification_compositions[name]
    except KeyError:
        return modification_compositions[glycoct_map[name]]


class Modification(ModificationBase):
    def __init__(self, name, position, composition=None):
        composition = composition or resolve_composition_rule(name)(position)
        self.name = name
        self.position = position
        self.composition = composition

    def __repr__(self):  # pragma: no cover
        return "<Modification {name}@{position}>".format(
            name=self.name, position=self.position)

    def __eq__(self, other):
        return (self.name == other) or\
            (self.name == other.name and self.composition == other.composition)

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.name)

    def to_glycoct(self):
        return glycoct_map.get(self.name, self.name)
