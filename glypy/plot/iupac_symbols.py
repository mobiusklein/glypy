from matplotlib.path import Path
from matplotlib.textpath import TextPath
import matplotlib.patches as patches
from .common import MonosaccharidePatch
from .geometry import centroid, Affine2D
from glypy.io import iupac

default_line_color = 'lightgrey'
default_line_color = 'black'
default_scale_factor = 1.0
zorder = 2
line_weight = 0.5


def draw(monosaccharide, x, y, ax, tree_node=None, scale=0.1, **kwargs):
    try:
        text = iupac.to_iupac(monosaccharide)
    except:
        text = monosaccharide
    if text[0] == 'a':
        pass
    elif text[0] == 'b':
        pass
    node = draw_text(ax, x, y, text)
    return MonosaccharidePatch(saccharide_label=(node,))


def line_to(ax, sx, sy, ex, ey, zorder=1, color=default_line_color):
    vertices = [
        (sx, sy),
        (ex, ey),
        (0, 0)
    ]
    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.STOP
    ]
    path = Path(vertices, codes)
    patch = patches.PathPatch(path, edgecolor=default_line_color, alpha=0.5)
    ax.add_patch(patch)
    return patch


def draw_text(ax, x, y, text, scale=0.1, **kwargs):
    fs = kwargs.get("fontsize", 2) * scale * 0.75
    t_path = TextPath((x - 0.1, y), s=text, size=fs)
    center = centroid(t_path)
    dist = -1 * (center - (x, y))
    t_path_shifted = t_path.transformed(Affine2D().translate(*dist))
    patch = patches.PathPatch(t_path_shifted, facecolor="black", lw=line_weight / 20., zorder=4)
    a = ax.add_patch(patch)
    return a
