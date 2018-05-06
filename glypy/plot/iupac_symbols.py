
from matplotlib.textpath import TextPath
import matplotlib.patches as patches

from glypy.io import iupac

from .common import MonosaccharidePatch
from .geometry import centroid, Affine2D
from .symbolic_nomenclature import SymbolicNomenclatureBase

default_line_color = 'lightgrey'
# default_line_color = 'black'
default_scale_factor = 1.0
zorder = 2
line_weight = 0.5


class IUPACTextSymbolicNomenclature(SymbolicNomenclatureBase):
    default_line_color = 'lightgrey'

    def draw_text(self, ax, x, y, text, scale=0.1, **kwargs):
        fs = kwargs.get("fontsize", 2) * scale * 0.13
        t_path = TextPath((x - 0.1, y), s=text, size=fs)
        center = centroid(t_path)
        dist = -1 * (center - (x, y))
        t_path_shifted = t_path.transformed(Affine2D().translate(*dist))
        patch = patches.PathPatch(t_path_shifted, facecolor="black", lw=line_weight / 20., zorder=4)
        a = ax.add_patch(patch)
        return a

    def draw(self, monosaccharide, x, y, ax, tree_node=None, scale=0.1, **kwargs):
        try:
            text = iupac.to_iupac(monosaccharide)
        except Exception:
            text = monosaccharide
        if text[0] == 'a':
            text = text.rsplit("-", 1)[1]
        elif text[0] == 'b':
            text = text.rsplit("-", 1)[1]
        node = self.draw_text(ax, x, y, text, scale=scale * 3.)
        return MonosaccharidePatch(saccharide_label=(node,))
