from glypy.composition.structure_composition import substituent_compositions
from .base import SubstituentBase
from .link import Link

from glypy.utils import uid
from glypy.composition import Composition, calculate_mass
from glypy.utils.multimap import OrderedMultiMap


class DerivatizePathway(object):

    def __init__(self, can_nh_derivatize=False, is_nh_derivatizable=False, name=None):
        self.can_nh_derivatize = can_nh_derivatize
        self.is_nh_derivatizable = is_nh_derivatizable
        if name is not None:
            derivatize_info[name] = self

    def __repr__(self):  # pragma: no cover
        return "<DerivatizePathway {}>".format(self.__dict__)

    @classmethod
    def register(cls, name, can_nh_derivatize, is_nh_derivatizable):
        derivatize_info[name.replace("-", "_")] = DerivatizePathway(can_nh_derivatize, is_nh_derivatizable)


attachment_composition_info = {
    "sulfate": Composition("H"),
    "methyl": Composition("H"),
    "n_acetyl": Composition("OH"),
    "n_glycolyl": Composition("OH"),
    "amino": Composition("OH")
}
default_attachment_composition = Composition("H")

derivatize_info = {
    "acetyl": DerivatizePathway(True, False),
    "amino": DerivatizePathway(True, False),
    "anhydro": DerivatizePathway(True, False),
    "bromo": DerivatizePathway(True, False),
    "chloro": DerivatizePathway(True, False),
    "diphospho_ethanolamine": DerivatizePathway(True, False),
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
    "phospho_ethanolamine": DerivatizePathway(True, False),
    "pyrophosphate": DerivatizePathway(True, False),
    "pyruvate": DerivatizePathway(True, False),
    "succinate": DerivatizePathway(True, False),
    "sulfate": DerivatizePathway(True, False),
    "thio": DerivatizePathway(True, False),
    "triphosphate": DerivatizePathway(True, False),

    "(x)_lactate": DerivatizePathway(True, False),
    "(r)_carboxyethyl": DerivatizePathway(True, False),
    "(r)_lactate": DerivatizePathway(True, False),
    "(r)_pyruvate": DerivatizePathway(True, False),
    "(s)_carboxyethyl": DerivatizePathway(True, False),
    "(s)_lactate": DerivatizePathway(True, False),
    "(s)_pyruvate": DerivatizePathway(True, False),

    "n_acetyl": DerivatizePathway(True, True),
    "n_amidino": DerivatizePathway(True, True),
    "n_formyl": DerivatizePathway(True, True),
    "n_glycolyl": DerivatizePathway(True, True),
    "n_methyl": DerivatizePathway(True, True),
    "n_succinate": DerivatizePathway(True, True),
    "n_sulfate": DerivatizePathway(True, True),
    "n_dimethyl": DerivatizePathway(True, True),

    "phospho_choline": DerivatizePathway(True, False),
}


class Substituent(SubstituentBase):
    '''
    Represents a non-saccharide molecule commonly found bound to saccharide units.

    Attributes
    ----------
    name: |str|
        The name of the substituent, used to uniquely identify it.
    links: |OrderedMultiMap|
        All links to all molecules connected to this one.
    composition: |Composition|
        The chemical makeup of this molecule.
    attachment_composition: |Composition|
        The default cost of attaching this substituent to a |Monosaccharide|
    id: |int|
        A unique identifier number for this molecule.

    can_nh_derivatize: |bool|
        Whether this substituent will derivatize at an amine group.
    is_nh_derivatizable: |bool|
        Whether this substituent contains a derivatizable amine group.
    _derivatize: |bool|
        Whether this substituent was added by a derivatization process.

    _order: |int|
        The number of connections to this molecule. Mutated internally by |Link| objects,
        not for external use. See :meth:`order`.

    '''

    def __init__(self, name, links=None, composition=None, id=None,
                 can_nh_derivatize=None, is_nh_derivatizable=None, derivatize=False,
                 attachment_composition=None):
        if links is None:
            links = OrderedMultiMap()
        self.name = name
        self.links = links
        if composition is None:
            composition = substituent_compositions[self._name]
        elif composition is not None and self._name not in substituent_compositions:
            substituent_compositions[self._name] = composition.clone()
        self.composition = composition
        self.id = id or uid()
        self._order = self.order()
        try:
            if can_nh_derivatize is is_nh_derivatizable is None:
                derivatize_pathway = derivatize_info[self.name]
                self.can_nh_derivatize = derivatize_pathway.can_nh_derivatize
                self.is_nh_derivatizable = derivatize_pathway.is_nh_derivatizable
            else:
                self.can_nh_derivatize = can_nh_derivatize or False
                self.is_nh_derivatizable = is_nh_derivatizable or False
        except KeyError:
            self.can_nh_derivatize = can_nh_derivatize or False
            self.is_nh_derivatizable = is_nh_derivatizable or False
        self._derivatize = derivatize
        self.attachment_composition = attachment_composition if attachment_composition is not None\
            else attachment_composition_info.get(self.name, default_attachment_composition)

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

    def __repr__(self):  # pragma: no cover
        return "<Substituent {name}>".format(name=self._name)

    def __hash__(self):
        return hash((self.id, self.name))

    def __eq__(self, other):
        return (other is not None) and (self.name == other.name) and (self.composition == other.composition)

    def __ne__(self, other):
        return not self == other

    def open_attachment_sites(self):
        return [1, 2], 0

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
        Adds a :class:`~glypy.structure.substituent.Substituent` and associated :class:`~glypy.structure.link.Link`
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
            The elemental composition removed from ``self``. Defaults to ``H1``.
        child_loss: Composition or str
            The elemental composition removed from ``substituent``. Defaults to ``H1``.

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

    def drop_substituent(self, position, substituent=None, refund=True):
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
        :func:`glypy.composition.composition.calculate_mass`
        '''
        if charge == 0:
            mass = calculate_mass(
                self.composition, average=average, charge=0, mass_data=mass_data)
            for link_pos, child in self.children():
                mass += child.mass(average=average,
                                   charge=0, mass_data=mass_data)
        else:
            mass = self.total_composition().calc_mass(average=average, charge=charge, mass_data=mass_data)
        return mass

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._order = state.get("_order", len(self.links))
        self._derivatize = state.get("_derivatize", False)

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

        substituent = Substituent(self.name,
                                  can_nh_derivatize=self.can_nh_derivatize,
                                  is_nh_derivatizable=self.is_nh_derivatizable,
                                  id=self.id if prop_id else None,
                                  derivatize=self._derivatize,
                                  attachment_composition=self.attachment_composition)
        for pos, link in self.links.items():
            if link.is_child(self):
                continue
            sub = link.to(self)
            dup = sub.clone()
            link.clone(substituent, dup)
        return substituent

    def order(self):
        return len(self.links)

    def total_composition(self):
        '''
        Computes the sum of the composition of `self` and each of its linked
        :class:`~.substituent.Substituent`s

        Returns
        -------
        :class:`~glypy.composition.Composition`
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
        result = []
        for pos, link in self.links.items():
            if link.is_child(self):
                continue
            result.append((pos, link.child))
        return result

    def parents(self):
        '''
        Returns an iterator over the objects which are considered
        the ancestors of ``self``.
        '''
        result = []
        for pos, link in self.links.items():
            if link.is_parent(self):
                continue
            result.append((pos, link.parent))
        return result

    def attachment_composition_loss(self):
        return self.attachment_composition.clone()
