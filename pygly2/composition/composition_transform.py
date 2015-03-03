from ..structure import Substituent
from ..structure import Glycan
from ..structure import Monosaccharide
from ..structure import Modification

from .composition import Composition


def derivatize(saccharide, substituent):
    '''
    For each monosaccharide and viable substituent, attach the specified
    substituent.
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
        monosaccharide_obj.add_substituent(
            substituent.clone(), parent_loss=Composition(H=1),
            position=site, child_loss=Composition(H=1), child_position=1)
    for p, subst in monosaccharide_obj.substituents():
        if subst.is_nh_derivatizable and substituent.can_nh_derivatize:
            subst.add_substituent(substituent.clone(), position=2, child_position=1)
    reducing_end_pos = monosaccharide_obj.reducing_end
    if reducing_end_pos is not None:
        monosaccharide_obj.add_substituent(
            substituent.clone(), parent_loss=Composition(H=1), max_occupancy=3,
            position=reducing_end_pos, child_loss=Composition(H=1), child_position=1)

    for pos, mod in monosaccharide_obj.modifications.items():
        if mod == Modification.a:
            monosaccharide_obj.add_substituent(
                substituent.clone(), position=pos, parent_loss=Composition(H=1), max_occupancy=3,
                child_loss=Composition(H=1), child_position=1)
            if pos == reducing_end_pos:
                monosaccharide_obj.add_substituent(
                    substituent.clone(), position=reducing_end_pos, parent_loss=Composition(H=1), max_occupancy=4,
                    child_loss=Composition(H=1), child_position=1)


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

