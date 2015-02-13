from .base import SubstituentBase
from .structure_composition import substituent_compositions
from .link import Link

from ..composition import Composition, calculate_mass
from ..utils import enum
from ..utils.multimap import OrderedMultiMap


class SubstituentEnum(object):
    __metaclass__ = enum.EnumMeta
    acetyl = 6
    amino = 3
    anhydro = 8
    bromo = 19
    chloro = 22
    diphosphoethanolamine = 14
    ethanolamine = 30
    ethyl = 23
    fluoro = 18
    formyl = 16
    glycolyl = 20
    hydroxymethyl = 34
    iodo = 29
    lactone = 36
    methyl = 10
    n_acetyl = 1
    n_amidino = 28
    n_formyl = 25
    n_glycolyl = 5
    n_methyl = 11
    n_succinate = 12
    n_sulfate = 7
    ndimethyl = 21
    phosphate = 9
    phosphocholine = 39
    phosphoethanolamine = 4
    pyrophosphate = 15
    pyruvate = 13
    r_carboxyethyl = 37
    r_lactate = 32
    r_pyruvate = 27
    s_carboxyethyl = 38
    s_lactate = 31
    s_pyruvate = 26
    succinate = 24
    sulfate = 2
    thio = 17
    triphosphate = 35
    x_lactate = 33


class Substituent(SubstituentBase):
    '''
    Represents a non-saccharide molecule commonly found bound to saccharide units.
    '''

    def __init__(self, name, links=None, composition=None):
        if links is None:
            links = OrderedMultiMap()
        self.name = name
        self.links = links
        if composition is None:
            composition = substituent_compositions[self.name]
        self.composition = composition

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = value.replace("-", "_")

    def to_glycoct(self):
        return "s:{0}".format(self.name.replace("_", "-"))

    def __repr__(self):  # pragma: no cover
        return "<Substituent {name}>".format(name=self._name)

    def __eq__(self, other):
        return (other is not None) and (self.name == other.name) and (self.composition == other.composition)

    def __ne__(self, other):
        return not self == other

    def is_occupied(self, position):
        if position > 2 or position < 1:
            raise IndexError("Position out of range")
        return len(self.links[position])

    def add_substituent(self, substitent, position=-1, max_occupancy=1,
                        child_position=-1, parent_loss=None, child_loss=None):
        if self.is_occupied(position) > max_occupancy:
            raise IndexError("Site is already occupied")
        if parent_loss is None:
            parent_loss = Composition(H=1)
        if child_loss is None:
            child_loss = Composition(H=1)
        Link(parent=self, child=substitent,
             parent_position=position, child_position=child_position,
             parent_loss=parent_loss, child_loss=child_loss)
        return self

    def drop_substiuent(self, position, substituent=None, refund=True):
        link_obj = None
        for substituent_link in self.links[position]:
            if substituent_link.child == substituent or substituent is None:
                link_obj = substituent_link
                break
        if link_obj is None:
            raise IndexError("No matching substituent found at {position}".format(position=position))

        link_obj.break_link(refund=refund)
        return self

    def mass(self, average=False, charge=0, mass_data=None):
        '''
        Calculates the total mass of `self` and all nodes returned by :meth:`children`.

        Parameters
        ----------
        average: bool, optional, defaults to False
            Whether or not to use the average isotopic composition when calculating masses.
            When ``average == False``, masses are calculated using monoisotopic mass.
        charge: int, optional, defaults to 0
            If charge is non-zero, m/z is calculated, where m is the theoretical mass, and z is `charge`
        mass_data: dict, optional, defaults to `None`
            If mass_data is None, standard NIST mass and isotopic abundance data are used. Otherwise the
            contents of mass_data are assumed to contain elemental mass and isotopic abundance information.

        Returns
        -------
        :class:`float`

        See also
        --------
        :func:`pygly2.composition.composition.calculate_mass`
        '''
        mass = calculate_mass(self.composition, average=average, charge=charge, mass_data=mass_data)
        for link_pos, child in self.children():
            mass += child.mass(average=average, charge=charge, mass_data=mass_data)
        return mass

    def clone(self):
        substituent = Substituent(self.name)
        for pos, link in self.links.items():
            if link.is_child(self):
                continue
            sub = link.to(self)
            dup = sub.clone()
            Link(substituent, dup, link.parent_position, link.child_position,
                 link.parent_loss, link.child_loss)
        return substituent

    def total_composition(self):
        '''
        Computes the sum of the composition of `self` and each of its linked
        :class:`~.substituent.Substituent`s

        Returns
        -------
        :class:`~pygly2.composition.Composition`
        '''
        comp = self.composition
        for pos, sub in self.children():
            comp = comp + sub.total_composition()
        return comp

    def children(self):
        '''
        Returns an iterator over the :class:`Monosaccharide`s which are considered
        the descendants of ``self``.
        '''
        for pos, link in self.links.items():
            if link.is_child(self):
                continue
            yield (pos, link.child)

    def parents(self):
        '''
        Returns an iterator over the :class:`Monosaccharide`s which are considered
        the ancestors of ``self``.
        '''
        for pos, link in self.links.items():
            if link.is_parent(self):
                continue
            yield (pos, link.parent)
