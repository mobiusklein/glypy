'''

'''

from .base import ModificationBase
from glypy.composition.structure_composition import modification_compositions
from glypy.io.format_constants_map import modification_map
from glypy.utils import uid

glycoct_map = {v: k for k, v in modification_map.items()}


def resolve_composition_rule(name):
    try:
        return modification_compositions[name]
    except KeyError:
        return modification_compositions[glycoct_map[name]]


class Modification(ModificationBase):
    __slots__ = ('name', 'position', 'composition', 'id')

    def __init__(self, name, position, composition=None, id=None):
        composition = composition or resolve_composition_rule(name)(position)
        if id is None:
            id = uid()
        self.name = name
        self.position = position
        self.composition = composition
        self.id = id

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

    def clone(self, prop_id=True):
        return self.__class__(self.name, self.position, self.composition, id=self.id if prop_id else None)

    def __reduce__(self):
        return self.__class__, (self.name, self.position, self.composition, self.id)

    def site_count(self):
        return 1


class MultiSiteModification(Modification):
    __slots__ = ('positions',)

    def __init__(self, name, positions, composition=None, id=None):
        super(MultiSiteModification, self).__init__(
            name, positions[0], composition=composition, id=id)
        self.positions = tuple(positions)

    def clone(self, prop_id=True):
        return self.__class__(self.name, self.positions, self.composition, id=self.id if prop_id else None)

    def __reduce__(self):
        return self.__class__, (self.name, self.positions, self.composition, self.id)

    def site_count(self):
        return len(self.positions)

    def __repr__(self):  # pragma: no cover
        return "<MultiSiteModification {name}@{positions}>".format(
            name=self.name, positions=self.positions)
