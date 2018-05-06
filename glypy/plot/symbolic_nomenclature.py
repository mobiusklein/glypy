
import numpy as np
import matplotlib

from matplotlib.path import Path
from matplotlib.textpath import TextPath
import matplotlib.patches as patches

from .common import MonosaccharidePatch

from glypy.utils.enum import Enum


default_line_color = 'black'
default_scale_factor = 1.0
zorder = 2
line_weight = 0.5


class UnknownShapeException(Exception):
    pass


def draw_circle(self, ax, x, y, color, scale=0.1):
    path = Path(Path.unit_circle().vertices * scale, Path.unit_circle().codes)
    trans = matplotlib.transforms.Affine2D().translate(x, y)
    t_path = path.transformed(trans)
    patch = patches.PathPatch(
        t_path, facecolor=color.value, lw=line_weight, zorder=2)
    a = ax.add_patch(patch)
    ma = MonosaccharidePatch(saccharide_shape=(a,))
    return ma


def draw_square(self, ax, x, y, color, scale=0.1):
    square_verts = np.array([
        (0.5, 0.5),
        (0.5, -0.5),
        (-0.5, -0.5),
        (-0.5, 0.5),
        (0.5, 0.5),
        (0., 0.),
    ]) * 2
    square_codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]
    path = Path(square_verts * scale, square_codes)
    trans = matplotlib.transforms.Affine2D().translate(x, y)
    t_path = path.transformed(trans)
    patch = patches.PathPatch(
        t_path, facecolor=color.value, lw=line_weight, zorder=2)
    a = ax.add_patch(patch)
    ma = MonosaccharidePatch(saccharide_shape=(a,))
    return ma


def draw_triangle(self, ax, x, y, color, scale=0.1):
    path = Path(Path.unit_triangle().vertices * scale, Path.unit_triangle().codes)
    trans = matplotlib.transforms.Affine2D().translate(
        x, y).rotate_deg_around(x, y, -90)
    t_path = path.transformed(trans)
    patch = patches.PathPatch(
        t_path, facecolor=color.value, lw=line_weight, zorder=2)
    a = ax.add_patch(patch)
    ma = MonosaccharidePatch(saccharide_shape=(a,))
    return ma


class SymbolicNomenclatureBase(object):
    default_line_color = 'black'

    def draw(monosaccharide, x, y, ax, tree_node=None, scale=0.1, **kwargs):
        raise NotImplementedError()

    def draw_text(self, ax, x, y, text, scale=0.1, **kwargs):
        fs = kwargs.get("fontsize", 2) * scale * .25
        t_path = TextPath((x, y), s=text, size=fs)
        patch = patches.PathPatch(t_path, facecolor="black", lw=line_weight / 20., zorder=4)
        a = ax.add_patch(patch)
        return (a,)

    def format_text(self, text):
        label = ''.join([token.capitalize()[:2] for token in text.split('_', 1)])
        return label

    def line_to(self, ax, sx, sy, ex, ey, zorder=1, color=None, lw=1):
        if color is None:
            color = self.default_line_color
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
        patch = patches.PathPatch(path, color=color, lw=lw, zorder=zorder)
        ax.add_patch(patch)
        return patch
