import itertools
from uuid import uuid4
try:
    from collections.abc import Iterable
except ImportError:
    from collections import Iterable

from glypy.composition import Composition
from glypy.utils import uid, basestring, make_struct
from .base import SaccharideBase, SubstituentBase
from .constants import UnknownPosition, LinkageType


default_parent_loss = Composition({"O": 1, "H": 1})
default_child_loss = Composition(H=1)


linkage_configuration = make_struct("linkage_configuration", ("parent", "child", "parent_position", "child_position"))


class Link(object):
    '''
    Represents the linkage between two molecules, described as an edge in a graph with optional
    directedness semantics.

    Attributes
    ----------
    parent: :class:`~.Monosaccharide` or :class:`~.Substituent`
    child: :class:`~.Monosaccharide` or :class:`~.Substituent`
    parent_position: int
        Position of link on :attr:`parent`
    child_position: int
        position of link on :attr:`child`
    parent_loss: |Composition|
        Elemental composition lost from the :attr:`parent` when forming this bond
    child_loss: |Composition|
        Elemental composition lost from the :attr:`child` when forming this bond
    _attached: :class:`bool`
        An internal flag that keeps track of whether or not a link is attached to
        its termini. Should not be read or manipulated directly. Instead see :meth:`is_attached`.
    '''

    __slots__ = (
        "parent", "child", "parent_position", "child_position",
        "parent_loss", "child_loss", "id",
        "label", "_attached",
        "parent_linkage_type",
        "child_linkage_type"
    )

    def __init__(self, parent, child, parent_position=UnknownPosition, child_position=UnknownPosition,
                 parent_loss=None, child_loss=None, id=None, attach=True, parent_linkage_type=None,
                 child_linkage_type=None):
        '''
        Defines a bond between `parent` and `child` between the molecule positions specified.
        The bond may represent a partial loss of elemental composition from the parent and/or child
        molecules.

        Instantiating the |Link| object will automatically attach it to its endpoints, mutating them
        unless `attach=False`. If not attached on instantiation, the bond can be realized by calling
        :meth:`Link.apply` at a later time.

        Parameters
        ----------
        parent: :class:`~.Monosaccharide` or :class:`~.Substituent`
        child: :class:`~.Monosaccharide` or :class:`~.Substituent`
        parent_position: int
            The position on the parent to attach to Defaults to :const:`~.UnknownPosition`
        child_position: int
            The position on the child to attach to. Defaults to :const:`~.UnknownPosition`
        parent_loss: :class:`Composition` or str
            The elemental composition deducted from the parent when the bond is applied
        child_loss: :class:`Composition` or str
            The elemental composition deducted from the child when the bond is applied
        id: int
            A locally unique identifier within a graph. If |None|, :func:`~.glypy.utils.uid` is used to generate one. Defaults to |None|
        attach: bool
            Whether to immediately attach the |Link| object to the `parent` and `child` molecules on instantiation
            by using :meth:`Link.apply`
        '''

        if id is None:
            id = uid()

        if parent_loss is None:
            parent_loss = default_parent_loss
        elif isinstance(parent_loss, basestring):
            parent_loss = Composition(formula=parent_loss)
        if child_loss is None:
            child_loss = default_child_loss
        elif isinstance(child_loss, basestring):
            child_loss = Composition(formula=child_loss)

        self.parent = parent
        self.child = child
        self.parent_position = parent_position
        self.child_position = child_position
        self.parent_loss = parent_loss
        self.child_loss = child_loss
        self.parent_linkage_type = parent_linkage_type or LinkageType.unknown
        self.child_linkage_type = child_linkage_type or LinkageType.unknown
        self.id = id
        self.label = None
        self._attached = False
        if attach:
            self.apply()

    def has_ambiguous_linkage(self):
        """Returns whether the link position at either terminal is ambiguous, but not unknown.

        Returns
        -------
        bool
        """
        return False

    def has_ambiguous_termini(self):
        """Returns whether the link terminal is ambiguous

        Returns
        -------
        bool
        """
        return False

    def is_ambiguous(self):
        """Returns if the link position or link terminal at either end is ambiguous,
        but not unknown.

        Returns
        -------
        bool

        See Also
        --------
        :meth:`has_ambiguous_linkage`
        :meth:`has_ambiguous_termini`
        """
        return self.has_ambiguous_linkage() or self.has_ambiguous_termini()

    def apply(self):
        '''
        Actually alter :attr:`parent` and :attr:`child`'s state. Deduct their respective losses from their
        :attr:`composition` and insert a reference to :obj:`self` into their :attr:`links` or :attr:`substituent_links`
        as appropriate.

        This performs no position availability validation. If there is a possibility that either the parent or child
        position may be occupied, they should be tested using :meth:`~.Monosaccharide.is_occupied` prior to calling
        this method.

        Sets :attr:`_attached` to |True|

        See also
        --------
        :meth:`Link.break_link`
        :meth:`Link.refund`
        :meth:`Monosaccharide.is_occupied`
        '''
        # assert not self.is_attached(), ("Cannot apply an already attached link")
        self.parent.composition -= (self.parent_loss)

        self.child.composition -= (self.child_loss)
        if self.is_substituent_link():
            self.parent.substituent_links[self.parent_position] = self
        else:
            self.parent.links[self.parent_position] = self
        self.child.links[self.child_position] = self
        self.parent._degree += 1
        self.child._degree += 1
        self._attached = True

    def to(self, mol):
        '''
        Traverse the link from `mol` to its adjacent node. If `mol` is not `self.parent` or
        `self.child`, instead raises an error. :meth:`__getitem__` is an alias of `to`

        Parameters
        ----------
        mol: :class:`~.Monosaccharide` or :class:`~.Substituent`

        Returns
        -------
        :class:`~.Monosaccharide` or :class:`~.Substituent`

        Raises
        ------
        KeyError
            If `mol` is not the parent or child of this |Link|
        '''
        if mol is (self.parent):
            return self.child
        elif mol is (self.child):
            return self.parent
        else:
            raise KeyError("Could not find connection for {0}".format(mol))

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
        return (self.child.node_type is SubstituentBase.node_type) and\
            (self.parent.node_type is SaccharideBase.node_type)

    def is_bridge_link(self):
        '''
        If :attr:`parent` is a |Substituent| and :attr:`child` is a |Monosaccharide|, then `self`
        is a *bridge_link*

        Returns
        -------
        bool
        '''
        return (self.parent.node_type is SubstituentBase.node_type) and\
            (self.child.node_type is SaccharideBase.node_type)

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
                    attach=attach,
                    parent_linkage_type=self.parent_linkage_type,
                    child_linkage_type=self.child_linkage_type)

    def _get_parent_configuration(self):
        # Possible problem since Monosaccharide objects do not define __hash__
        return {self.parent}, {self.parent_position}

    def _get_child_configuration(self):
        # Possible problem since Monosaccharide objects do not define __hash__
        return {self.child}, {self.child_position}

    def __eq__(self, other):
        '''
        Performs deep equality testing between `self` and `other`
        '''
        if (other is None or not isinstance(other, Link)):
            return False
        if self.is_ambiguous() or other.is_ambiguous():
            res = self._get_parent_configuration() == other._get_parent_configuration()
            if res:
                res = self._get_child_configuration() == other._get_child_configuration()
            if res:
                res = self._flat_equality(other)
            if res:
                res = (self.parent.id == other.parent.id) and\
                      (self.parent == other.parent)
            return res
        else:
            res = self._flat_equality(other)
            if res:
                res = (self.parent.id == other.parent.id) and\
                      (self.parent == other.parent)
            return res

    def trait_equality(self, other):
        """Test whether the intrinsic attributes of another :class:`Link` object
        match this link without examining the actual termini.

        Parameters
        ----------
        other : :class:`Link`
            The object to compare against

        Returns
        -------
        bool
        """
        res = True
        if res:
            res = (self.parent_position == other.parent_position)
        if res:
            res = (self.child_position == other.child_position)
        if res:
            res = (self.parent_loss == other.parent_loss)
        if res:
            res = (self.child_loss == other.child_loss)
        return res

    def _flat_equality(self, other):
        res = self.trait_equality(other)
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
        self.parent._degree -= 1
        self.child._degree -= 1
        self._attached = False
        return (self.parent, self.child)

    def _reconnect(self, refund=False):
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

        self.child.links[self.child_position] = self

        if refund:
            self.refund()
        self.parent._degree += 1
        self.child._degree += 1
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

    def _glycoct_sigils(self):
        '''
        Helper method for determining which GlycoCT symbols and losses to present
        '''
        parent_loss_str = 'x'
        child_loss_str = 'x'
        water = Composition({"O": 1, "H": 1})

        if self.child_loss == water:
            child_loss_str = "d"
            parent_loss_str = "o"
        elif self.parent_loss == water:
            child_loss_str = 'o'
            parent_loss_str = 'd'

        if self.child_loss == Composition(
                H=1) and (self.child.node_type is SubstituentBase.node_type):
            child_loss_str = "n"
            if self.parent_loss == water:
                parent_loss_str = "d"
            else:
                parent_loss_str = "o"

        if self.child_loss is None:
            child_loss_str = 'x'
        if self.parent_loss is None:
            parent_loss_str = 'x'

        return parent_loss_str, child_loss_str

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

    def __getstate__(self):
        state = {}
        state["parent"] = self.parent
        state["child"] = self.child
        state["parent_position"] = self.parent_position
        state["child_position"] = self.child_position
        state["parent_loss"] = self.parent_loss
        state["child_loss"] = self.child_loss
        state["id"] = self.id
        state["label"] = self.label
        state["_attached"] = self._attached
        state['parent_linkage_type'] = self.parent_linkage_type
        state['child_linkage_type'] = self.child_linkage_type
        return state

    def __setstate__(self, state):
        self.parent = state["parent"]
        self.child = state["child"]
        self.parent_position = state["parent_position"]
        self.child_position = state["child_position"]
        self.parent_loss = state["parent_loss"]
        self.child_loss = state["child_loss"]
        self.id = state["id"]
        self.label = state["label"]
        self._attached = state["_attached"]
        try:
            self.parent_linkage_type = state['parent_linkage_type']
        except KeyError:
            self.parent_linkage_type = LinkageType.unknown
        try:
            self.child_linkage_type = state['child_linkage_type']
        except KeyError:
            self.child_linkage_type = LinkageType.unknown


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
    '''A :class:`Link` that may have more than one option for connection position at either
    terminal, or the terminals themselves may be chosen from a set of options.

    Attributes
    ----------
    parent_choices: list
        The set of possible choices for :attr:`parent`
    parent_position_choices: list
        The set of possible choices for :attr:`parent_position`
    child_choices: list
        The set of possible choices for :attr:`child`
    child_position_choices: list
        The set of possible choices for :attr:`child_position`
    '''

    __slots__ = (
        "parent_choices", "parent_position_choices",
        "child_choices", "child_position_choices"
    )

    def __init__(self, parent, child, parent_position=(UnknownPosition,), child_position=(UnknownPosition,),
                 parent_loss=None, child_loss=None, id=None, attach=True, parent_linkage_type=None,
                 child_linkage_type=None):
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
            parent_loss, child_loss, id, attach,
            parent_linkage_type, child_linkage_type)

    def __getstate__(self):
        state = super(AmbiguousLink, self).__getstate__()
        state['parent_choices'] = self.parent_choices
        state['parent_position_choices'] = self.parent_position_choices
        state['child_choices'] = self.child_choices
        state['child_position_choices'] = self.child_position_choices
        return state

    def __setstate__(self, state):
        super(AmbiguousLink, self).__setstate__(state)
        self.parent_choices = state['parent_choices']
        self.parent_position_choices = state['parent_position_choices']
        self.child_choices = state['child_choices']
        self.child_position_choices = state['child_position_choices']

    def has_ambiguous_linkage(self):
        return len(self.parent_position_choices) > 1 or len(self.child_position_choices) > 1

    def has_ambiguous_termini(self):
        return len(self.parent_choices) > 1 or len(self.child_choices) > 1

    def iterconfiguration(self, attach=False):
        """Iterate over possible configurations for this link, yielding them.

        If ``attach`` is :const:`True`, call :meth:`reconfigure` with each
        specification.

        Parameters
        ----------
        attach: :class:`bool`, optional
            Whether or not to apply each configuration. The default is :const:`False`

        Yields
        ------
        :class:`tuple` of :attr:`parent`, :attr:`child`, :attr:`parent_position`, :attr:`child_position`
        """
        configurations = self._configuration_crossproduct()
        for parent_o, child_o, parent_position_o, child_position_o in configurations:
            if attach:
                self.reconfigure(parent_o, child_o, parent_position_o, child_position_o)
            yield linkage_configuration(parent_o, child_o, parent_position_o, child_position_o)

    def _configuration_crossproduct(self):
        configurations = itertools.product(self.parent_choices, self.child_choices,
                                           self.parent_position_choices,
                                           self.child_position_choices)
        return configurations

    def reconfigure(self, parent, child, parent_position, child_position):
        """Apply the specified configuration, replacing the current configuration.

        This method calls :meth:`break_link` to refund the current bonds and
        then calls :meth:`apply`.

        Parameters
        ----------
        parent : :class:`~.Monosaccharide` or :class:`~.Substituent`
            The parent of the link to form
        child : :class:`~.Monosaccharide` or :class:`~.Substituent
            The child of the link to form
        parent_position : int
            The position on the parent to connect to
        child_position : int
            The position on the child to connect to
        """
        self.break_link(refund=True)
        self.parent = parent
        self.child = child
        self.parent_position = parent_position
        self.child_position = child_position
        self.apply()

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

    def __repr__(self):
        rep = super(AmbiguousLink, self).__repr__()
        return "(A)" + rep

    def _find_open_position_single(self, parent, child):
        open_parent_sites, _ = parent.open_attachment_sites()
        parent_has_unknown_position = False
        if UnknownPosition in open_parent_sites:
            open_parent_sites += list(set(self.parent_position_choices) - set(parent.occupied_attachment_sites()))
            parent_has_unknown_position = True
        parent_site_options = list((set(open_parent_sites) | {UnknownPosition}) & set(self.parent_position_choices))
        if parent_site_options:
            parent_site = parent_site_options[0]
        else:
            if parent_has_unknown_position and open_parent_sites:
                parent_site = UnknownPosition
            else:
                return None
        open_child_sites, _ = child.open_attachment_sites()
        child_has_unknown_position = False
        if UnknownPosition in open_child_sites:
            open_child_sites += list(set(self.child_position_choices) - set(child.occupied_attachment_sites()))
            child_has_unknown_position = True
        child_site_options = list((set(open_child_sites) | {UnknownPosition}) & set(self.child_position_choices))
        if child_site_options:
            child_site = child_site_options[0]
        else:
            if child_has_unknown_position and open_child_sites:
                child_site = UnknownPosition
            else:
                return None
        return parent_site, child_site

    def _find_open_position_multiple(self, parent, child):  # pragma: no cover
        '''Debugging purposes only
        '''
        open_parent_sites, _ = parent.open_attachment_sites()
        if UnknownPosition in open_parent_sites:
            open_parent_sites += self.parent_position_choices
        parent_site_options = list((set(open_parent_sites) | {UnknownPosition}) & set(self.parent_position_choices))
        open_child_sites, _ = child.open_attachment_sites()
        if UnknownPosition in open_child_sites:
            open_child_sites += self.child_position_choices
        child_site_options = list((set(open_child_sites) | {UnknownPosition}) & set(self.child_position_choices))
        return parent_site_options, child_site_options

    def _find_open_position_combinations(self):
        for parent, child in itertools.product(self.parent_choices, self.child_choices):
            positions = self._find_open_position_single(parent, child)
            if positions is not None:
                return (parent, positions[0], child, positions[1])

    def find_open_position(self, attach=True):
        """Identify a valid configuration of parent, child and respective positions
        which are all open.

        Parameters
        ----------
        attach : bool, optional
            Whether or not to actually connect to the selected configuration, which will
            break the current configuration (the default is :const:`True`).

        Raises
        ------
        ValueError:
            When no valid configuration can be found.

        """
        if attach:
            self.break_link(refund=True)
        configuration = self._find_open_position_combinations()
        if configuration is None:
            raise ValueError("Could not find a valid configurations on current parent/child pair")
        parent, parent_position, child, child_position = configuration
        self.parent = parent
        self.child = child
        self.parent_position = parent_position
        self.child_position = child_position
        if attach:
            self.apply()
