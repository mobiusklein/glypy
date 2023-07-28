import re
from collections import namedtuple
from typing import List, Tuple, Union, DefaultDict

from six import string_types as basestring

from glypy.composition import Composition
from glypy.structure import Monosaccharide, Substituent, Link, substituent as _substituent
from glypy.io.tree_builder_utils import try_int

from .carbon_descriptors import CarbonDescriptors
from .substituent_conversion import alin_to_substituent, substituent_to_alin


_HYDROXYL = Composition("OH")

_NodeTypeSpec = namedtuple("NodeTypeSpec", [
    "carbon_descriptor",
    "substituents"
])


def _mapspec_position_parse(token: str):
    if token == '?':
        return -1
    return int(token.split("n")[0])



def translate_mapspec_position(map_spec: str, i: int):
    return [_mapspec_position_parse(tok) for tok in map_spec[:i].split('-')]


class NodeTypeSpec(_NodeTypeSpec):
    @classmethod
    def parse(cls, text, version_number):
        if version_number < 2.0:
            raise TypeError("Cannot parse type version earlier than 2.0")
        parts = text.split("_")
        skeleton_anomer = parts[0]
        substituents = []
        if len(parts) > 1:
            ring_specification_or_substituents = parts[1]
            if re.match(r"[0-9\-1?]+-[0-9\-1?]+", ring_specification_or_substituents):
                ring_specification = ring_specification_or_substituents
            else:
                ring_specification = '?-?'
                substituents.append(ring_specification_or_substituents)
        else:
            ring_specification = '?-?'
        substituents.extend(parts[2:])
        try:
            skeleton, anomer = skeleton_anomer.split("-")
        except Exception:
            skeleton = skeleton_anomer
            anomer = (-1, None)
        ring_start, ring_end = map(try_int, ring_specification.split("-"))
        anomeric_position, anomer_type = anomer
        return cls(CarbonDescriptors(skeleton, anomer_type, anomeric_position, ring_start, ring_end),
                   cls.translate_substituents(substituents))

    @classmethod
    def from_monosaccharide(cls, monosaccharide: Monosaccharide):
        descriptors = CarbonDescriptors.from_monosaccharide(monosaccharide)
        substituents = []
        subsituent_to_positions = DefaultDict(list)
        substituent_by_id = {}
        for position, link in monosaccharide.substituent_links.items():
            dest: Substituent = link.to(monosaccharide)
            subsituent_to_positions[dest.id].append((position, link.parent_loss == _HYDROXYL))
            substituent_by_id[dest.id] = dest
        for k, link_specs in subsituent_to_positions.items():
            dest = substituent_by_id[k]
            if len(link_specs) == 1:
                position, deoxy = link_specs[0]
            else:
                position = [spec[0] for spec in link_specs]
                # Doesn't account for the case where multiple attachment points are deoxy
                deoxy = any(spec[1] for spec in link_specs)
            substituents.append((position, dest.name, deoxy))
        substituents = sorted(
            substituents, key=lambda x: x[0] if isinstance(x[0], int) else x[0][0], reverse=True)
        return cls(descriptors, substituents)

    @classmethod
    def translate_substituents(self, substituents: List[str]) -> List[Tuple[Union[int, List[int]],
                                                                            str, bool]]:
        result = []
        for map_spec in substituents:
            start_i = map_spec.index("*")
            positions = translate_mapspec_position(map_spec, start_i)
            map_code = map_spec[start_i + 1:]
            # Doesn't account for the case where multiple attachment points are deoxy
            record, deoxy = alin_to_substituent(map_code)
            result.append((positions[0] if len(positions) == 1 else positions, record.name, deoxy))
        return result

    def to_monosaccharide(self):
        base: Monosaccharide = self.carbon_descriptor.to_base_type()
        child_loss = Composition("H")
        substituents_ordered = sorted(self.substituents, key=lambda x: x[0] if isinstance(x[0], int) else x[0][0])
        for position, subst, deoxy in substituents_ordered:
            if deoxy:
                parent_loss = Composition("OH")
            else:
                parent_loss = Composition("H")
            if isinstance(position, list):
                subst = Substituent(subst)
                for pi in position:
                    Link(
                        base, subst, parent_position=pi, child_position=1, parent_loss=parent_loss,
                        child_loss=child_loss
                    )
            else:
                base.add_substituent(subst, position, parent_loss=parent_loss, child_loss=child_loss,
                                     max_occupancy=100)
        return base

    def to_res(self):
        mods = []
        for position, substituent, deoxy in self.substituents:
            alin = substituent_to_alin(substituent, deoxy)
            if isinstance(position, list):
                position = [str(p) if p != -1 else '?' for p in position]
                position = '-'.join(position)
            else:
                if position == -1:
                    position = "?"
            mods.append("%s*%s" % (position, alin))
        encoded = '_'.join([self.carbon_descriptor.to_backbone_code()] + mods)
        return encoded

    def __str__(self):
        return self.to_res()

    def __repr__(self):
        return super(NodeTypeSpec, self).__repr__()

    def __eq__(self, other):
        if other is None:
            return False
        if isinstance(other, basestring):
            return str(self) == other
        if self.carbon_descriptor != other.carbon_descriptor:
            return False
        elif self.substituents != other.substituents:
            return False
        return True

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.carbon_descriptor)
