from uuid import uuid4

from ..composition import Composition
from .base import SaccharideBase, SubstituentBase

default_parent_loss = Composition(O=1, H=1)
default_child_loss = Composition(H=1)


class Link(object):

    '''
    Represents the linkage between two molecules, described as an edge in a graph with optional
    directedness semantics.

    Attributes
    ----------
    parent: :class:`~pygly2.structure.monosaccharide.Monosaccharide` or :class:`~.substituent.Substituent`
    child: :class:`~pygly2.structure.monosaccharide.Monosaccharide` or :class:`~.substituent.Substituent`
    parent_position: int
        Position of link on :attr:`parent`
    child_position: int
        position of link on :attr:`child`
    parent_loss: |Composition|
        Elemental composition lost from the :attr:`parent` when forming this bond
    child_loss: |Composition|
        Elemental composition lost from the :attr:`child` when forming this bond

    '''

    def __init__(self, parent, child, parent_position=-1, child_position=-1,
                 parent_loss=None, child_loss=None, id=None, attach=True):
        '''
        Defines a bond between `parent` and `child` between the molecule positions specified.
        The bond may represent a partial loss of elemental composition from the parent and/or child
        molecules.

        Instantiating the |Link| object will automatically attach it to its endpoints, mutating them
        unless `attach=False`. If not attached on instantiation, the bond can be realized by calling
        :meth:`Link.apply()` at a later time.

        Parameters
        ----------
        parent: :class:`Monosaccharide` or :class:`Substituent`
        child: :class:`Monosaccharide` or :class:`Substituent`
        parent_position: int
            The position on the parent to attach to Defaults to -1
        child_position: int
            The position on the child to attach to. Defaults to -1
        parent_loss: :class:`Composition` or str
            The elemental composition deducted from the parent when the bond is applied
        child_loss: :class:`Composition` or str
            The elemental composition deducted from the child when the bond is applied
        id: int
            A locally unique identifier within a graph. If |None|, uuid4 is used to generate one. Defaults to |None|
        attach: bool
            Whether to immediately attach the |Link| object to the `parent` and `child` molecules on instantiation
            by using :meth:`Link.apply`
        '''

        if isinstance(parent_loss, basestring):
            parent_loss = Composition(formula=parent_loss)
        if isinstance(child_loss, basestring):
            child_loss = Composition(formula=child_loss)

        self.parent = parent
        self.child = child
        self.parent_position = parent_position
        self.child_position = child_position
        self.parent_loss = parent_loss
        self.child_loss = child_loss
        self.id = id or uuid4().int
        self.label = None

        if attach:
            self.apply()

    def apply(self):
        '''
        Actually alter :attr:`parent` and :attr:`child`'s state. Deduct their respective losses from their
        :attr:`composition` and insert a reference to :obj:`self` into their :attr:`links` or :attr:`substituent_links`
        as appropriate.

        See also
        --------
        :meth:`Link.break_link`
        :meth:`Link.reconnect`
        :meth:`Link.refund`

        '''
        self.parent.composition -= (self.parent_loss or default_parent_loss)

        self.child.composition -= (self.child_loss or default_child_loss)
        if self.is_substituent_link():
            self.parent.substituent_links[self.parent_position] = self
        else:
            self.parent.links[self.parent_position] = self
        self.child.links[self.child_position] = self

    def to(self, mol):
        '''
        Traverse the link from `mol` to its adjacent node. If `mol` is not `self.parent` or
        `self.child`, instead raises an error. :meth:`__getitem__` is an alias of `to`

        Parameters
        ----------
        mol: :class:`Monosaccharide` or :class:`Substituent`

        Returns
        -------
        :class:`Monosaccharide` or :class:`Substituent`

        Raises
        ------
        KeyError:
            If `mol` is not the parent or child of this |Link|
        '''
        hmol = id(mol)
        if hmol == id(self.parent):
            return self.child
        elif hmol == id(self.child):
            return self.parent
        else:
            raise KeyError("Could not find connection for {0}".format(mol))

    def _glycoct_sigils(self):
        '''
        Helper method for determining which GlycoCT symbols and losses to present
        '''
        parent_loss_str = 'x'
        child_loss_str = 'x'
        if self.child_loss == Composition(O=1, H=1):
            child_loss_str = "d"
            parent_loss_str = "o"
        elif self.parent_loss == Composition(O=1, H=1):
            child_loss_str = 'o'
            parent_loss_str = 'd'

        if self.child_loss == Composition(
                H=1) and isinstance(self.child, SubstituentBase):
            child_loss_str = "n"
            if self.parent_loss == Composition(O=1, H=1):
                parent_loss_str = "d"
            else:
                parent_loss_str = "o"

        if self.child_loss is None:
            child_loss_str = 'x'
        if self.parent_loss is None:
            parent_loss_str = 'x'

        return parent_loss_str, child_loss_str

    def to_glycoct(self, ix, parent_ix, child_ix):
        '''
        Serializes `self` as a Condensed GlycoCT LIN entry. Depends upon
        the textual indices of its parent and child lines.

        See also
        --------
        :meth:`pygly2.structure.monosaccharide.Monosaccharide.to_glycoct`
        :meth:`pygly2.structure.glycan.Glycan.to_glycoct`
        '''

        parent_loss_str, child_loss_str = self._glycoct_sigils()
        rep = "{ix}:{parent_ix}{parent_loss}({parent_position}+{child_position}){child_ix}{child_loss}"
        return rep.format(
            ix=ix,
            parent_ix=parent_ix,
            parent_loss=parent_loss_str,
            parent_position=self.parent_position,
            child_ix=child_ix,
            child_loss=child_loss_str,
            child_position=self.child_position)

    #: Alias for :meth:`to`
    __getitem__ = to

    def __iter__(self):
        '''
        Yields the `parent` node, then the `child` node
        '''
        yield self.parent
        yield self.child

    def is_substituent_link(self):
        '''
        If :attr:`child` is a |Substituent| and :attr:`parent` is a |Monosaccharide|, then `self`
        is a *substituent_link*

        Returns
        -------
        bool
        '''
        return isinstance(self.child, SubstituentBase) and isinstance(
            self.parent, SaccharideBase)

    def is_parent(self, mol):
        '''
        Test `mol` for :func:`id` equality with :attr:`parent`

        Returns
        -------
        bool
        '''
        return id(mol) == id(self.parent)

    def is_child(self, mol):
        '''
        Test `mol` for :func:`id` equality with :attr:`child`

        Returns
        -------
        bool
        '''
        return id(mol) == id(self.child)

    def clone(self, parent, child, prop_id=True, attach=True):
        '''
        Given two new objects `parent` and `child`, duplicate all other information in `self`,
        creating a new `Link` object between `parent` and `child` with the same properties as
        `self`.

        Returns
        -------
        Link
        '''
        return Link(parent, child,
                    parent_position=self.parent_position,
                    child_position=self.child_position,
                    parent_loss=self.parent_loss,
                    child_loss=self.child_loss,
                    id=self.id if prop_id else uuid4().int,
                    attach=attach)

    def __eq__(self, other):
        '''
        Performs deep equality testing between `self` and `other`
        '''
        if (other is None or not isinstance(other, Link)):
            return False
        res = (self.parent == other.parent)
        if res:
            res = (self.child == other.child)
        if res:
            res = (self.parent_position == other.parent_position)
        if res:
            res = (self.child_position == other.child_position)
        if res:
            res = (self.parent_loss == other.parent_loss)
        if res:
            res = (self.child_loss == other.child_loss)
        return res

    def __ne__(self, other):
        return not (self == other)

    def break_link(self, refund=False):
        '''
        This function reverses :meth:`Link.apply`, removing the reference to :obj:`self` from
        both :attr:`parent` and :attr:`child`'s store of links. If ``refund`` is :const:`True`,
        :meth:`Link.refund` will be called, restoring the lost elemental composition from this
        bond.

        Returns
        -------
        :class:`~.monosaccharide.Monosaccharide` or :class:`~.substituent.Substituent` parent
        :class:`~.monosaccharide.Monosaccharide` or :class:`~.substituent.Substituent` child
        '''
        if self.is_substituent_link():
            self.parent.substituent_links.pop(self.parent_position, self)
        else:
            self.parent.links.pop(self.parent_position, self)
        self.child.links.pop(self.child_position, self)
        if refund:
            self.refund()
        return (self.parent, self.child)

    def reconnect(self, refund=False):
        '''
        The opposite of :meth:`Link.break_link`, add `self` to the appropriate places on
        :attr:`parent` and :attr:`child`

        Parameters
        ----------
        refund: bool
            Should :meth:`Link.refund` be called? Defaults to |False|

        '''
        if self.is_substituent_link():
            self.parent.substituent_links[self.parent_position] = self
        else:
            self.parent.links[self.parent_position] = self
            sorted(
                self.parent.links[self.parent_position], key=lambda x: x.child.order())

        self.child.links[self.child_position] = self

        if refund:
            self.refund()

    def refund(self):
        '''
        Returns the lost elemental composition caused by :meth:`apply`. Adds back :attr:`parent_loss`
        and :attr:`child_loss` to the :attr:`composition` of :attr:`parent` and :attr:`child` respectively
        '''
        self.parent.composition += (self.parent_loss or default_parent_loss)
        self.child.composition += (self.child_loss or default_child_loss)

    def is_attached(self):
        '''
        Test to see if `self` is present in :attr:`parent` and :attr:`child` link structures.
        Presences indicates the link is attached.

        Returns
        -------
        bool
        '''
        if self.is_substituent_link():
            parent = self in self.parent.substituent_links[self.parent_position]
        else:
            parent = self in self.parent.links[self.parent_position]
        child = self in self.child.links[self.child_position]
        return parent and child

    def __repr__(self):  # pragma: no cover
        parent_loss_str, child_loss_str = self._glycoct_sigils()
        rep = "({parent.id}){parent_loss}({parent_position}+{child_position}){child_loss}({child.id})"
        return rep.format(
            parent=self.parent,
            parent_loss=parent_loss_str,
            parent_position=self.parent_position,
            child=self.child,
            child_loss=child_loss_str,
            child_position=self.child_position)


def glycocidic_bond(parent, child, parent_position, child_position):
    '''A convenient shortcut for constructing glycans'''
    link = Link(parent, child, parent_position=parent_position,
                child_position=child_position,
                parent_loss=Composition(formula="H"),
                child_loss=Composition(formula="OH"))
    return link


class LinkMaskContext(object):  # pragma: no cover
    '''
    A context manager for masking and unmasking |Link| objects on a residue
    '''
    def __init__(self, residue, attach=False):
        self.residue = residue
        self.links = [link for link in residue.links.values()]
        self.attach = attach

    def __enter__(self):
        self.mask() if not self.attach else self.unmask()

    def __exit__(self, exc_type, exc_value, traceback):
        self.mask() if self.attach else self.unmask()

    def mask(self):
        [link.break_link(refund=True) for link in self.links]

    def unmask(self):
        [link.apply() for link in self.links]
