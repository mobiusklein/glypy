from uuid import uuid4

from .base import SubstituentBase
from ..composition.structure_composition import substituent_compositions
from .link import Link

from ..composition import Composition, calculate_mass
from ..utils.multimap import OrderedMultiMap


class DerivatizePathway(object):

    def __init__(self, can_nh_derivatize=False, is_nh_derivatizable=False):
        self.can_nh_derivatize = can_nh_derivatize
        self.is_nh_derivatizable = is_nh_derivatizable

    def __repr__(self):  # pragma: no cover
        return "<DerivatizePathway {}>".format(self.__dict__)

    @classmethod
    def register(cls, name, can_nh_derivatize, is_nh_derivatizable):
        derivatize_info[name.replace("-", "_")] = DerivatizePathway(can_nh_derivatize, is_nh_derivatizable)


derivatize_info = {
    "acetyl": DerivatizePathway(True, False),
    "amino": DerivatizePathway(True, False),
    "anhydro": DerivatizePathway(True, False),
    "bromo": DerivatizePathway(True, False),
    "chloro": DerivatizePathway(True, False),
    "diphosphoethanolamine": DerivatizePathway(True, False),
    "ethanolamine": DerivatizePathway(True, False),
    "ethyl": DerivatizePathway(True, False),
    "fluoro": DerivatizePathway(True, False),
    "formyl": DerivatizePathway(True, False),
    "glycolyl": DerivatizePathway(True, False),
    "hydroxymethyl": DerivatizePathway(True, False),
    "iodo": DerivatizePathway(True, False),
    "lactone": DerivatizePathway(True, False),
    "methyl": DerivatizePathway(True, False),
    "phosphate": DerivatizePathway(True, False),
    "phosphocholine": DerivatizePathway(True, False),
    "phosphoethanolamine": DerivatizePathway(True, False),
    "pyrophosphate": DerivatizePathway(True, False),
    "pyruvate": DerivatizePathway(True, False),
    "succinate": DerivatizePathway(True, False),
    "sulfate": DerivatizePathway(True, False),
    "thio": DerivatizePathway(True, False),
    "triphosphate": DerivatizePathway(True, False),

    "x_lactate": DerivatizePathway(True, False),
    "r_carboxyethyl": DerivatizePathway(True, False),
    "r_lactate": DerivatizePathway(True, False),
    "r_pyruvate": DerivatizePathway(True, False),
    "s_carboxyethyl": DerivatizePathway(True, False),
    "s_lactate": DerivatizePathway(True, False),
    "s_pyruvate": DerivatizePathway(True, False),

    "n_acetyl": DerivatizePathway(True, True),
    "n_amidino": DerivatizePathway(True, True),
    "n_formyl": DerivatizePathway(True, True),
    "n_glycolyl": DerivatizePathway(True, True),
    "n_methyl": DerivatizePathway(True, True),
    "n_succinate": DerivatizePathway(True, True),
    "n_sulfate": DerivatizePathway(True, True),
    "n_dimethyl": DerivatizePathway(True, True),
}


class Substituent(SubstituentBase):

    '''
    Represents a non-saccharide molecule commonly found bound to saccharide units.
    '''

    def __init__(self, name, links=None, composition=None, id=None, can_nh_derivatize=None, is_nh_derivatizable=None):
        if links is None:
            links = OrderedMultiMap()
        self.name = name
        self.links = links
        if composition is None:
            composition = substituent_compositions[self.name]
        self.composition = composition
        self.id = id or uuid4().int
        try:
            if can_nh_derivatize is is_nh_derivatizable is None:
                self.can_nh_derivatize = derivatize_info[self.name].can_nh_derivatize
                self.is_nh_derivatizable = derivatize_info[self.name].is_nh_derivatizable
            else:
                self.can_nh_derivatize = can_nh_derivatize or False
                self.is_nh_derivatizable = is_nh_derivatizable or False
        except KeyError:
            raise KeyError("{} does not have defined derivatization rules. Please specify them.".format(self.name))

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        '''
        Translate the name of the substituent from the common dash-separated notation
        to a valid identifier, replacing - with _.

        Parameters
        ----------
        value: str
            The name being set

        '''
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
        '''
        Check to see if `position` is occupied. Unlike |Monosaccharide|, |Substituent| objects
        can only have two attachment sites at this time.

        Parameters
        ----------
        position: int

        Returns
        -------
        int:
            Number of links at `position`

        Raises
        ------
        IndexError:
            If `position` > 2 or < 1
        '''
        if position > 2 or position < 1:
            raise IndexError("Position out of range")
        return len(self.links[position])

    def add_substituent(self, substitent, position=2, max_occupancy=1,
                        child_position=1, parent_loss=None, child_loss=None):
        '''
        Adds a :class:`~pygly2.structure.substituent.Substituent` and associated :class:`~pygly2.structure.link.Link`
        to :attr:`links` at the site given by ``position``. This new substituent is included when calculating mass
        with substituents included

        Parameters
        ----------
        substituent: str or Substituent
            The substituent to add. If passed a |str|, it will be
            translated into an instance of |Substituent|
        position: int or 'x'
            The location to add the |Substituent| link to :attr:`links`. Defaults to 2
        child_position: int
            The location to add the link to in `substituent`'s :attr:`links`. Defaults to 1. Substituent
            indices are currently not checked.
        max_occupancy: int, optional
            The maximum number of items acceptable at ``position``. Defaults to :const:`1`
        parent_loss: Composition or str
            The elemental composition removed from ``self``
        child_loss: Composition or str
            The elemental composition removed from ``substituent``

        Raises
        ------
        ValueError:
            ``position`` is occupied by more than ``max_occupancy`` elements
        '''

        if self.is_occupied(position) > max_occupancy:
            raise ValueError("Site is already occupied")
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
            raise IndexError(
                "No matching substituent found at {position}".format(position=position))

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
        mass = calculate_mass(
            self.composition, average=average, charge=charge, mass_data=mass_data)
        for link_pos, child in self.children():
            mass += child.mass(average=average,
                               charge=charge, mass_data=mass_data)
        return mass

    def clone(self, prop_id=True):
        '''
        Duplicates this |Substituent| object, recursively copying all children
        as well.

        Parameters
        ----------
        prop_id: bool
            Whether or not to propagate :attr:`id` to the clone.


        Returns
        -------
        Substituent

        See Also
        --------
        :meth:`.structure.Monosaccharide.clone`
        '''
        substituent = Substituent(self.name)
        for pos, link in self.links.items():
            if link.is_child(self):
                continue
            sub = link.to(self)
            dup = sub.clone()
            Link(substituent, dup, link.parent_position, link.child_position,
                 link.parent_loss, link.child_loss)
            substituent.id = self.id if prop_id else uuid4().int
        return substituent

    def order(self):
        return len(self.links)

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
        Returns an iterator over the objects which are considered
        the ancestors of ``self``.
        '''
        for pos, link in self.links.items():
            if link.is_parent(self):
                continue
            yield (pos, link.parent)
