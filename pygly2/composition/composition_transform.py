from ..structure import Substituent
from ..structure import Glycan
from ..structure import Monosaccharide
from .composition import Composition


def derivatize(saccharide, substituent):
    if isinstance(substituent, basestring):
        substituent = Substituent(substituent)
    if isinstance(saccharide, Glycan):
        for node in saccharide:
            derivatize(node, substituent)
    elif isinstance(saccharide, Monosaccharide):
        derivatize_monosaccharide(saccharide, substituent)

def derivatize_monosaccharide(monosaccharide_obj, substituent):
    open_sites, unknowns = monosaccharide_obj.open_attachment_sites()
    for site in open_sites:
        monosaccharide_obj.add_substituent(
            substituent.clone(), parent_loss=Composition(H=1),
            parent_position=site, child_loss=Composition(H=1), child_position=1)
    for subst in monosaccharide_obj.substituents():
        if subst.name in permethylatable_substituents:
            subst.add_substituent(substituent.clone(), parent_position=2, child_position=1)
    reducing_end_pos = monosaccharide_obj.reducing_end
    if reducing_end_pos is not None:
        monosaccharide_obj.add_substituent(
            substituent.clone(), parent_loss=Composition(H=1), max_occupancy=3,
            parent_position=reducing_end_pos, child_loss=Composition(H=1), child_position=1)



permethylatable_substituents = {
    "n_acetyl",
    "n_glycolyl",
}


class DerivatizeBase(object):

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

