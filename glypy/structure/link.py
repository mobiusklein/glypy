import itertools
from uuid import uuid4
from collections import Iterable

from glypy.composition import Composition
from .base import SaccharideBase, SubstituentBase

default_parent_loss = Composition(O=1, H=1)
default_child_loss = Composition(H=1)


class Link(object):

    '''
    Represents the linkage between two molecules, described as an edge in a graph with optional
    directedness semantics.

    Attributes
    ----------
    parent: :class:`~glypy.structure.monosaccharide.Monosaccharide` or :class:`~.substituent.Substituent`
    child: :class:`~glypy.structure.monosaccharide.Monosaccharide` or :class:`~.substituent.Substituent`
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
        self._attached = False
        if attach:
            self.apply()

    def apply(self):
        '''
        Actually alter :attr:`parent` and :attr:`child`'s state. Deduct their respective losses from their
        :attr:`composition` and insert a reference to :obj:`self` into their :attr:`links` or :attr:`substituent_links`
        as appropriate.

        Sets :attr:`_attached` to |True|

        See also
        --------
        :meth:`Link.break_link`
        :meth:`Link.reconnect`
        :meth:`Link.refund`

        '''
        assert not self.is_attached(),("Cannot apply an already attached link")
        self.parent.composition -= (self.parent_loss or default_parent_loss)

        self.child.composition -= (self.child_loss or default_child_loss)
        if self.is_substituent_link():
            self.parent.substituent_links[self.parent_position] = self
        else:
            self.parent.links[self.parent_position] = self
        self.child.links[self.child_position] = self
        self.parent._order += 1
        self.child._order += 1
        self._attached = True

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
        if mol is (self.parent):
            return self.child
        elif mol is (self.child):
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
                H=1) and (self.child.node_type is SubstituentBase.node_type):
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
        :meth:`glypy.structure.monosaccharide.Monosaccharide.to_glycoct`
        :meth:`glypy.structure.glycan.Glycan.to_glycoct`
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
        # return isinstance(self.child, SubstituentBase) and isinstance(
        #     self.parent, SaccharideBase)
        return (self.child.node_type is SubstituentBase.node_type) and\
            (self.parent.node_type is SaccharideBase.node_type)

    def is_parent(self, mol):
        '''
        Test `mol` for :func:`id` equality with :attr:`parent`

        Returns
        -------
        bool
        '''
        return (mol) is (self.parent)

    def is_child(self, mol):
        '''
        Test `mol` for :func:`id` equality with :attr:`child`

        Returns
        -------
        bool
        '''
        return (mol) is (self.child)

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
        res = self._flat_equality(other)
        if res:
            res = (self.parent.id == other.parent.id) and\
                  (self.parent == other.parent)
        return res

    def _flat_equality(self, other):
        res = True
        if res:
            res = (self.parent_position == other.parent_position)
        if res:
            res = (self.child_position == other.child_position)
        if res:
            res = (self.parent_loss == other.parent_loss)
        if res:
            res = (self.child_loss == other.child_loss)
        if res:
            res = (self.child.id == other.child.id) and\
                  (self.child == other.child)
        return res

    def __ne__(self, other):
        return not (self == other)

    def break_link(self, refund=False):
        '''
        This function reverses :meth:`Link.apply`, removing the reference to :obj:`self` from
        both :attr:`parent` and :attr:`child`'s store of links. If ``refund`` is :const:`True`,
        :meth:`Link.refund` will be called, restoring the lost elemental composition from this
        bond.

        Sets :attr:`_attached` to |False|

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
        self.parent._order -= 1
        self.child._order -= 1
        self._attached = False
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
            # sorted(
            #     self.parent.links[self.parent_position], key=lambda x: x.child.order())

        self.child.links[self.child_position] = self

        if refund:
            self.refund()
        self.parent._order += 1
        self.child._order += 1
        self._attached = True

    def refund(self):
        '''
        Returns the lost elemental composition caused by :meth:`apply`. Adds back :attr:`parent_loss`
        and :attr:`child_loss` to the :attr:`composition` of :attr:`parent` and :attr:`child` respectively
        '''
        self.parent.composition += (self.parent_loss or default_parent_loss)
        self.child.composition += (self.child_loss or default_child_loss)

    def is_attached(self, deep=False):
        '''
        Test to see if `self` is present in :attr:`parent` and :attr:`child` link structures.
        Presences indicates the link is attached.

        Returns
        -------
        |bool|
        '''
        if not deep:
            return self._attached
        parent = False
        child = False
        if self.is_substituent_link():
            parent_list = self.parent.substituent_links[self.parent_position]
        else:
            parent_list = self.parent.links[self.parent_position]
        for lin in parent_list:
            if lin.id == self.id:
                if lin == self:
                    parent = True
                    break
        for lin in self.child.links[self.child_position]:
            if lin.id == self.id:
                if lin == self:
                    child = True
                    break
        assert self._attached == (parent and child)
        res = self._attached = (parent and child)
        return res

    def try_break_link(self, refund=False):  # pragma: no cover
        '''
        Try to break the link if it is attached, otherwise return |False|
        '''
        if self.is_attached():
            return self.break_link(refund=refund)
        return False

    def try_apply(self):  # pragma: no cover
        '''
        Try to apply the link if it is not attached, otherwise return |False|
        '''
        if not self.is_attached():
            return self.apply()
        return False

    def __repr__(self):  # pragma: no cover
        parent_loss_str, child_loss_str = self._glycoct_sigils()
        rep = "({parent.id}){parent_loss}({parent_position}+{child_position}){child_loss}({child.id})"
        rep = rep.format(
            parent=self.parent,
            parent_loss=parent_loss_str,
            parent_position=self.parent_position,
            child=self.child,
            child_loss=child_loss_str,
            child_position=self.child_position)
        if self.is_attached():
            rep += "[x]"
        else:
            rep += "[ ]"
        return rep


class LinkMaskContext(object):
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
        """
        For each |Link| captured, if that link is attched, call its
        :meth:`break_link` method with `refund=True`
        """
        for link in self.links:
            if link.is_attached():
                link.break_link(refund=True)

    def unmask(self):
        """
        For each |Link| captured, if that link is not attched, call its
        :meth:`apply` method
        """
        for link in self.links:
            if not link.is_attached():
                link.apply()


class AmbiguousLink(Link):
    def __init__(self, parent, child, parent_position=(-1,), child_position=(-1,),
                 parent_loss=None, child_loss=None, id=None, attach=True):
        if not isinstance(parent, (list, tuple)):
            parent = [parent]
        if not isinstance(child, (list, tuple)):
            child = [child]
        if not isinstance(parent_position, Iterable):
            parent_position = [parent_position]
        if not isinstance(child_position, Iterable):
            child_position = [child_position]

        self.parent_choices = list(parent)
        self.parent_position_choices = list(parent_position)

        self.child_choices = list(child)
        self.child_position_choices = list(child_position)

        super(AmbiguousLink, self).__init__(
            self.parent_choices[0],
            self.child_choices[0],
            self.parent_position_choices[0],
            self.child_position_choices[0],
            parent_loss, child_loss, id, attach)

    def iterconfiguration(self):  # pragma: no cover
        configurations = itertools.product(self.parent_choices, self.child_choices,
                                           self.parent_position_choices,
                                           self.child_position_choices)
        for parent_o, child_o, parent_position_o, child_position_o in configurations:
            self.break_link(refund=True)
            self.parent = parent_o
            self.child = child_o
            self.parent_position = parent_position_o
            self.child_position = child_position_o
            self.apply()
            yield parent_o, child_o, parent_position_o, child_position_o

    def clone(self, parent, child, prop_id=True, attach=True):
        if not isinstance(parent, (list, tuple)):
            parent = [parent]
        if not isinstance(child, (list, tuple)):
            child = [child]
        link = AmbiguousLink(
            parent, child,
            self.parent_position_choices,
            self.child_position_choices,
            self.parent_loss, self.child_loss, self.id if prop_id else None,
            attach=False)
        link.parent_position = self.parent_position
        link.child_position = self.child_position
        if attach:
            link.apply()

        return link

    def break_link(self, refund=True):
        res = super(AmbiguousLink, self).break_link(refund=refund)
        assert not self.is_attached(deep=True)
        return res

    def to_glycoct(self, ix, parent_ix, child_ix):  # pragma: no cover
        '''
        Serializes `self` as a Condensed GlycoCT LIN entry. Depends upon
        the textual indices of its parent and child lines.

        See also
        --------
        :meth:`glypy.structure.monosaccharide.Monosaccharide.to_glycoct`
        :meth:`glypy.structure.glycan.Glycan.to_glycoct`
        '''

        parent_loss_str, child_loss_str = self._glycoct_sigils()
        rep = "{ix}:{parent_ix}{parent_loss}({parent_position}+{child_position}){child_ix}{child_loss}"
        return rep.format(
            ix=ix,
            parent_ix=parent_ix,
            parent_loss=parent_loss_str,
            parent_position='|'.join(map(str, self.parent_position_choices)),
            child_ix=child_ix,
            child_loss=child_loss_str,
            child_position='|'.join(map(str, self.child_position_choices)))

    def __repr__(self):
        rep = super(AmbiguousLink, self).__repr__()
        return "(A)" + rep

    def find_open_position(self, attach=True):
        if attach:
            self.break_link(refund=True)
        try:
            open_parent_sites, _ = self.parent.open_attachment_sites()
            parent_site = iter(set(open_parent_sites) & set(self.parent_position_choices)).next()
            open_child_sites, _ = self.child.open_attachment_sites()
            child_site = iter(set(open_child_sites) & set(self.child_position_choices)).next()
            self.parent_position = parent_site
            self.child_position = child_site
            if attach:
                self.apply()
        except StopIteration:  # pragma: no cover
            raise ValueError("Could not find a valid configurations on current parent/child pair")
