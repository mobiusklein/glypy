import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib


from glypy.io import iupac

default_line_color = 'grey'
default_scale_factor = 5.0


def draw(monosaccharide, x, y, ax, tree_node=None, scale=0.1, **kwargs):
    try:
        text = iupac.to_iupac(monosaccharide)
    except:
        text = monosaccharide
    if text[0] == 'a':
        pass
    elif text[0] == 'b':
        pass
    node = ax.text(x, y-0.01, text, ha='center', fontweight=900)
    return (node,), ()


def line_to(ax, sx, sy, ex, ey, zorder=1, color='grey'):
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
    patch = patches.PathPatch(path, edgecolor=color, alpha=0.5)
    ax.add_patch(patch)
    return patch


def draw_text(ax, x, y, text, scale=0.1):
    a = ax.text(x=x, y=y, s=text, verticalalignment="center", horizontalalignment="center", fontsize=98 * scale)
    return (a,)
