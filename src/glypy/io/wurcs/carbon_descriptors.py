# Much code here is derived from https://github.com/glycoinfo/glycocttowurcs
# though the representation of glycans in Eurocarbdb.MolecularFramework
# may not map 1:1.
try:
    from collections.abc import Sequence
except ImportError:
    from collections import Sequence

import warnings

from six import string_types as basestring

from .basetype_conversion import (
    descriptors_to_base_type)

from glypy.structure.monosaccharide import Monosaccharide, ReducedEnd
from glypy.structure.constants import SuperClass, Anomer, Modification, Stem, Configuration, UnknownPosition
from glypy import OrderedMultiMap


anomer_map = {
    Anomer.beta: 'b',
    Anomer.alpha: 'a',
    Anomer.uncyclized: 'o',
    Anomer.x: 'x'
}


class CarbonDescriptors(Sequence):
    def __init__(self, descriptors, anomer, anomeric_position, ring_start, ring_end):
        self.descriptors = tuple(descriptors)
        self.anomer = Anomer[anomer]
        self.anomeric_position = self._translate_position(anomeric_position)
        self.ring_start = ring_start if ring_start is not None else UnknownPosition
        self.ring_end = ring_end if ring_end is not None else UnknownPosition

    def _translate_position(self, position):
        if position == '?':
            position = -1
        elif position == -1:
            position = '?'
        else:
            position = int(position)
        return position

    def __len__(self):
        return len(self.descriptors)

    def __eq__(self, other):
        if other is None:
            return False
        if isinstance(other, basestring):
            return str(self) == other
        if self.descriptors != other.descriptors:
            return False
        elif self.anomer != other.anomer:
            return False
        elif self.anomeric_position != other.anomeric_position:
            return False
        elif self.ring_start != other.ring_start:
            return False
        elif self.ring_end != other.ring_end:
            return False
        return True

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.descriptors)

    def __getitem__(self, i):
        return self.descriptors[i]

    def __iter__(self):
        return iter(self.descriptors)

    def to_d_stereoform(self, code):
        out = []
        is_l_stereoform = code[-1] == '3'
        for site in code:
            if is_l_stereoform:
                if site == '3':
                    out.append('4')
                elif site == '4':
                    out.append('3')
                else:
                    out.append(site)
            else:
                out.append(site)
        return out

    def to_base_type(self):
        '''Convert the :class:`CarbonDescriptors` into a
        :class:`~.Monosaccharide`, not including substituents.

        Returns
        -------
        :class:`~.Monosaccharide`
        '''
        superclass = SuperClass[len(self)]
        carbon_coding = list(map(str, self))
        modifications = OrderedMultiMap()
        is_reduced = False
        # translate stereocode into generic carbon code
        for i, site in enumerate(carbon_coding):
            if site == '1':
                carbon_coding[i] = '3'
            elif site == '2':
                carbon_coding[i] = '4'
        start = 1
        stems = []
        configurations = []
        anomer = self.anomer
        ring_start = self.ring_start
        ring_end = self.ring_end
        if carbon_coding[0] == carbon_coding[-1] == 'h':
            anomer = Anomer.uncyclized
            ring_start = 0
            ring_end = 0
            is_reduced = True
        # if the stereosites are all defined
        if 'x' not in carbon_coding:
            # incrementally walk along the carbon sequence
            while start < superclass.value:
                # consider ring stereoforms of up to four carbons ahead, preferring longer
                # stereosequences,
                for i in range(4, 0, -1):
                    # extract the raw stereosequence
                    raw_chunk = carbon_coding[start:start + i]
                    # convert the stereosequence to D configuration and
                    # convert to a string for hash lookup
                    chunk = ''.join(self.to_d_stereoform(raw_chunk))
                    try:
                        # if the look up is successful
                        stem_name = descriptors_to_base_type[chunk]
                        # save the mapped stem name
                        stems.append(stem_name)
                        # infer the chirality of the ring from the last
                        # stereosite
                        conf = Configuration.x
                        if raw_chunk[-1] == '3':
                            conf = Configuration.l
                        elif raw_chunk[-1] == '4':
                            conf = Configuration.d
                        configurations.append(conf)
                        # start the lookup process again from the next starting
                        # location
                        start += len(raw_chunk)
                        break
                    except KeyError:
                        continue
                else:
                    # if no stereosequence could be detected, if the start position
                    # is a stereosite, then we may have a grolene trilose component
                    if chunk in ('3', '4'):
                        stems.append(descriptors_to_base_type['x'])
                        # infer the chirality of the ring from the last
                        # stereosite
                        conf = Configuration.x
                        if raw_chunk[-1] == '3':
                            conf = Configuration.l
                        elif raw_chunk[-1] == '4':
                            conf = Configuration.d
                        configurations.append(conf)
                    start += 1
        else:
            # This cannot handle unspecified nonulonic acids and other modified but unspecified
            # monosaccharides with multiple chiral centers well.
            stems.append(None)
            if carbon_coding[0] in ('u', 'h'):
                configurations.append(None)
            else:
                warnings.warn("Cannot infer chirality from %r" % (str(self),))
                configurations.append(None)
            # Guess if the monosaccharide is large enough to have a second chiral center, because
            # no other rule seems obvious. This could produce incorrect monosaccharide compositions?
            if len(carbon_coding) > 6:
                stems.append(None)
                configurations.append(None)

        anomeric_position = None
        double_bonds = []
        for i, site in enumerate(self):
            if site in ('a', 'u', 'U'):
                anomeric_position = i + 1
                if anomeric_position == 2:
                    modifications[anomeric_position] = Modification.keto
            if site in ('E', 'F'):
                double_bonds.append(i + 1)
            if site in ('d', 'm'):
                modifications[i + 1] = Modification.Deoxygenated
            if site == 'A':
                modifications[i + 1] = Modification.Acidic
        for site in double_bonds[::2]:
            modifications[i] = Modification.en
        stems = [Stem[x] for x in stems[::-1]]
        configurations = configurations[::-1]

        base = Monosaccharide(
            anomer,
            configurations,
            stems,
            superclass,
            ring_start,
            ring_end,
            modifications, reduced=ReducedEnd() if is_reduced else None)
        return base

    @classmethod
    def from_monosaccharide(cls, monosaccharide):
        '''Create a :class:`CarbonDescriptors` from a given
        :class:`~.Monosaccharide`.

        Parameters
        ----------
        monosaccharide: :class:`~.Monosaccharide`
            The monosaccharide to describe

        Returns
        -------
        :class:`CarbonDescriptors`
        '''
        code = ['x'] * monosaccharide.superclass.value
        stereocode = monosaccharide.stereocode
        code = [str(x.value) if x.value is not None else 'x' for x in stereocode]
        code[0] = 'u'
        code[-1] = 'h'
        if monosaccharide.anomer == 'uncyclized':
            code[0] = 'h'
            code[-1] = 'h'
        anomer = monosaccharide.anomer
        anomeric_position = monosaccharide.ring_start
        anomeric_sites = []
        is_aldose = True
        # encode the modifications onto the carbon descriptor code
        for position, modification in monosaccharide.modifications.items():
            is_terminal = (position == 1 or position == monosaccharide.superclass.value)
            if modification == Modification.Acidic:
                if not is_terminal:
                    raise ValueError("Cannot add a carboxylic acid group to a non-terminal carbon")
                if position == 1:
                    is_aldose = False
                code[position - 1] = 'A'
            elif modification == Modification.Deoxygenated:
                if position == 1:
                    is_aldose = False
                if is_terminal:
                    code[position - 1] = 'm'
                else:
                    code[position - 1] = 'd'
            elif modification == Modification.Ketone:
                is_aldose = False
                # code[position] = 'o'
                anomeric_sites.append(position)
            elif modification == Modification.en:
                code[position - 1] = 'E'
            elif modification == Modification.Alditol:
                is_aldose = False
                if position != 1:
                    raise ValueError("\"aldi\" must occur on the first carbon")
        if is_aldose:
            anomeric_sites.append(1)
        anomeric_position = anomeric_sites[0]
        # if the anomeric position is fully defined and the monosaccharide is cyclic
        if monosaccharide.ring_start not in (UnknownPosition, 0):
            code[anomeric_position - 1] = 'a'
        # if the anomeric position is partially undefined, the carbon code is 'u'
        elif monosaccharide.ring_start == UnknownPosition:
            code[anomeric_position - 1] = 'u'
        if monosaccharide.ring_start == UnknownPosition:
            anomeric_position = "?"
        return cls(code, anomer, anomeric_position, monosaccharide.ring_start, monosaccharide.ring_end)

    def to_backbone_code(self):
        '''Convert :class:`CarbonDescriptors` into a string representation
        matching the ``<BackboneCode>`` pattern from WURCS2.0

        Returns
        -------
        :class:`str`
        '''
        parts = []
        # carbon descriptors
        parts.append(''.join([i for i in self]))
        # if the anomer is completely undefined, do not include it
        if not (self.anomeric_position == -1 and self.anomer == Anomer.x):
            parts.append("-%s%s" % (self._translate_position(self.anomeric_position),
                                    anomer_map[self.anomer]))
        # if the ring is neither undefined nor open, include it
        if (self.ring_start != -1 and self.ring_end != -1) and (self.ring_start != 0 and self.ring_end != 0):
            parts.append("_%s-%s" % tuple(map(self._translate_position, (self.ring_start, self.ring_end))))
        return ''.join(parts)

    def __str__(self):
        return self.to_backbone_code()

    def __repr__(self):
        descriptors = ''.join(map(str, self))
        template = ("{self.__class__.__name__}({descriptors!r}, {self.anomer.name!r}, "
                    "{self.anomeric_position}, {self.ring_start}, {self.ring_end})")
        return template.format(self=self, descriptors=descriptors)
