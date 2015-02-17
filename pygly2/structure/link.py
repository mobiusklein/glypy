from uuid import uuid4

from ..composition import Composition
from .base import SaccharideBase, SubstituentBase

default_parent_loss = Composition(O=1, H=1)
default_child_loss = Composition(H=1)


class Link(object):
    '''
    Represents the linkage between two residues.

    Attributes
    ----------
    parent: :class:`~pygly2.structure.monosaccharide.Monosaccharide` or :class:`~.substituent.Substituent`
    child: :class:`~pygly2.structure.monosaccharide.Monosaccharide` or :class:`~.substituent.Substituent`
    '''
    def __init__(self, parent, child, parent_position=-1, child_position=-1,
                 parent_loss=None, child_loss=None, id=None, **kwargs):

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
        self.data = kwargs

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

        self.parent.composition = self.parent.composition - (self.parent_loss or default_parent_loss)
        self.child.composition = self.child.composition - (self.child_loss or default_child_loss)
        if self.is_substituent_link():
            self.parent.substituent_links[self.parent_position] = self
        else:
            self.parent.links[self.parent_position] = self
        self.child.links[self.child_position] = self

    def to(self, mol):
        hmol = id(mol)
        if hmol == id(self.parent):
            return self.child
        elif hmol == id(self.child):
            return self.parent
        else:
            raise KeyError("Could not find connection for {0}".format(mol))

    def _glycoct_sigils(self):
        parent_loss_str = 'x'
        child_loss_str = 'x'
        if self.child_loss == Composition(O=1, H=1):
            child_loss_str = "d"
            parent_loss_str = "o"
        elif self.parent_loss == Composition(O=1, H=1):
            child_loss_str = 'o'
            parent_loss_str = 'd'

        if self.child_loss == Composition(H=1) and isinstance(self.child, SubstituentBase):
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

    __getitem__ = to

    def __iter__(self):
        yield self.parent
        yield self.child

    def is_substituent_link(self):
        return isinstance(self.child, SubstituentBase) and isinstance(self.parent, SaccharideBase)

    def is_parent(self, mol):
        return id(mol) == id(self.parent)

    def is_child(self, mol):
        return id(mol) == id(self.child)

    def clone(self, parent, child):
        return Link(parent, child,
                    parent_position=self.parent_position,
                    child_position=self.child_position,
                    parent_loss=self.parent_loss,
                    child_loss=self.child_loss,
                    id=self.id, **self.data)

    def __eq__(self, other):
        if (other is None):
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

    def break_link(self, refund=False, reorient_fn=None):
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
        if reorient_fn is not None:
            reorient_fn(self)
        if self.is_substituent_link():
            self.parent.substituent_links.pop(self.parent_position, self)
        else:
            self.parent.links.pop(self.parent_position, self)
        self.child.links.pop(self.child_position, self)
        if refund:
            self.refund()
        return (self.parent, self.child)

    def reconnect(self, refund=False, reorient_fn=None):
        if reorient_fn is not None:
            reorient_fn(self)
        if self.is_substituent_link():
            self.parent.substituent_links[self.parent_position] = self
        else:
            self.parent.links[self.parent_position] = self
            sorted(self.parent.links[self.parent_position], key=lambda x: x.child.order())

        self.child.links[self.child_position] = self

        if refund:
            self.refund()

    def refund(self):
        '''
        Returns the lost elemental composition caused by :meth:`apply`. Adds back :attr:`parent_loss`
        and :attr:`child_loss` to the :attr:`composition` of :attr:`parent` and :attr:`child` respectively
        '''
        self.parent.composition = self.parent.composition + (self.parent_loss or default_parent_loss)
        self.child.composition = self.child.composition + (self.child_loss or default_child_loss)

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
