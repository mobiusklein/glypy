from ..structure import Substituent
from ..structure import Glycan
from ..structure import Monosaccharide
from ..structure import Modification

from .composition import Composition


def derivatize(saccharide, substituent):
    '''
    For each monosaccharide and viable substituent, attach the specified
    substituent.

    .. warning::
        This function mutates the input `saccharide` object, because copying objects
        is expensive. If you want to continue using the underivatized object, create
        a copy of it using the object's :meth:`.clone` method.

    See Also
    --------
    :func:`.strip_derivitization`

    '''
    if isinstance(substituent, basestring):
        substituent = Substituent(substituent)
    if isinstance(saccharide, Glycan):
        for node in saccharide:
            derivatize(node, substituent)
    elif isinstance(saccharide, Monosaccharide):
        derivatize_monosaccharide(saccharide, substituent)
    return saccharide


def derivatize_monosaccharide(monosaccharide_obj, substituent):
    open_sites, unknowns = monosaccharide_obj.open_attachment_sites()
    for site in open_sites[unknowns:]:
        s = substituent.clone()
        s._derivatize = True
        monosaccharide_obj.add_substituent(
            s, parent_loss=Composition(H=1),
            position=site, child_loss=Composition(H=1), child_position=1)
    for p, subst in monosaccharide_obj.substituents():
        if subst.is_nh_derivatizable and substituent.can_nh_derivatize:
            s = substituent.clone()
            s._derivatize = True
            subst.add_substituent(s, position=2, child_position=1)
    red_end = monosaccharide_obj.reducing_end
    if red_end is not None:
        for i in range(1, red_end.valence + 1):
            s = substituent.clone()
            s._derivatize = True
            red_end.add_substituent(
                s, parent_loss=Composition(H=1), max_occupancy=3,
                position=i, child_loss=Composition(H=1), child_position=1)

    for pos, mod in monosaccharide_obj.modifications.items():
        if mod == Modification.a:
            s = substituent.clone()
            s._derivatize = True
            monosaccharide_obj.add_substituent(
                s, position=pos, parent_loss=Composition(H=1), max_occupancy=3,
                child_loss=Composition(H=1), child_position=1)


def strip_derivitization(saccharide):
    '''
    For each monosaccharide and viable substituent, remove all substituents
    added by :func:`derivatize`.

    .. warning::
        This function, like :func:`derivatize`, will mutate the input `saccharide`.

    See Also
    --------
    :func:`.derivatize`
    '''
    if isinstance(saccharide, Glycan):
        map(strip_derivatization_monosaccharide, saccharide)
    else:
        strip_derivatization_monosaccharide(saccharide)
    return saccharide


def strip_derivatization_monosaccharide(monosaccharide_obj):
    for pos, subst_link in monosaccharide_obj.substituent_links.items():
        if hasattr(subst_link.child, "_derivatize"):
            monosaccharide_obj.drop_substituent(pos, subst_link.child)
        else:
            sub_node = subst_link.child
            for sub_pos, subst_link in sub_node.links.items():
                if hasattr(subst_link.child, "_derivatize"):
                    sub_node.drop_substituent(sub_pos, subst_link.child)
    red_end = monosaccharide_obj.reducing_end
    if red_end is not None:
        for pos, subst_link in red_end.links.items():
            if hasattr(subst_link.child, "_derivatize"):
                red_end.drop_substituent(pos, subst_link.child)




# WIP
class DerivatizeBase(object):  # pragma: no cover
    def __init__(self, substituent, white_list, black_list, *args, **kwargs):
        self.substituent = substituent
        self.white_list = white_list
        self.black_list = black_list

    def on_reducing_end(self, node):
        pass

    def on_terminal(self, node):
        pass

    def on_substituent(self, node):
        pass

