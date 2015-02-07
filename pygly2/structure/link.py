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
                 parent_loss=None, child_loss=None, id=None):
        self.parent = parent
        self.child = child
        self.parent_position = parent_position
        self.child_position = child_position
        self.parent_loss = parent_loss
        self.child_loss = child_loss
        self.id = id

        self.apply()

    def apply(self):
        self.parent.composition = self.parent.composition - (self.parent_loss or default_parent_loss)
        self.child.composition = self.child.composition - (self.child_loss or default_child_loss)
        if self.is_substituent_link():
            self.parent.substituent_links[self.parent_position] = self
        else:
            self.parent.links[self.parent_position] = self
        self.child.links[self.child_position] = self

    def to(self, mol):
        hmol = hash(mol)
        if hmol == hash(self.parent):
            return self.child
        elif hmol == hash(self.child):
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
        rep = "{ix}:{parent_ix}{parent_loss}({parent_position}+{child_position}){child_ix}{child_loss};"
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
        return hash(mol) == hash(self.parent)

    def is_child(self, mol):
        return hash(mol) == hash(self.child)

    def clone(self, parent, child):
        return Link(parent, child, parent_position=self.parent_position,
                    child_position=self.child_position,
                    parent_loss=self.parent_loss,
                    child_loss=self.child_loss,
                    id=self.id)

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

    def break_link(self, reorient_fn=None, refund=False):
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

    def reconnect(self, reorient_fn=None):
        if reorient_fn is not None:
            reorient_fn(self)
        if self.is_substituent_link():
            self.parent.substituent_links[self.parent_position] = self
        else:
            self.parent.links[self.parent_position] = self
        self.child.links[self.child_position] = self

    def refund(self):
        self.parent.composition = self.parent.composition + (self.parent_loss or default_parent_loss)
        self.child.composition = self.child.composition + (self.child_loss or default_child_loss)

    def __repr__(self):
        parent_loss_str, child_loss_str = self._glycoct_sigils()
        rep = "({parent}){parent_loss}({parent_position}+{child_position}){child_loss}({child})"
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
                parent_loss=Composition(formula="OH"),
                child_loss=Composition(formula="H"))
    return link
