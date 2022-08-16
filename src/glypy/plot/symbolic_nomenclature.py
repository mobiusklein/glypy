from collections import namedtuple

import numpy as np
import matplotlib

from matplotlib.path import Path
from matplotlib.textpath import TextPath
import matplotlib.patches as patches
from matplotlib.transforms import Affine2D

from .common import MonosaccharidePatch
from .geometry import centroid, bounding_box

from glypy.utils.enum import Enum


default_line_color = 'black'
default_scale_factor = 1.0
zorder = 2
line_weight = 0.5


NodeLayoutTransforms = namedtuple("NodeLayoutTransforms", (
    'layout_transform', 'substituent_transform', 'annotation_transform'))


class UnknownShapeException(Exception):
    pass


def draw_circle(ax, x, y, color, scale=0.1):
    path = Path(Path.unit_circle().vertices * scale, Path.unit_circle().codes)
    trans = matplotlib.transforms.Affine2D().translate(x, y)
    t_path = path.transformed(trans)
    patch = patches.PathPatch(
        t_path, facecolor=color.value, lw=line_weight, zorder=2)
    a = ax.add_patch(patch)
    ma = MonosaccharidePatch(saccharide_shape=(a,))
    return ma


def draw_square(ax, x, y, color, scale=0.1):
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


def draw_triangle(ax, x, y, color, scale=0.1):
    path = Path(Path.unit_triangle().vertices * scale, Path.unit_triangle().codes)
    trans = matplotlib.transforms.Affine2D().translate(
        x, y).rotate_deg_around(x, y, -90)
    t_path = path.transformed(trans)
    patch = patches.PathPatch(
        t_path, facecolor=color.value, lw=line_weight, zorder=2)
    a = ax.add_patch(patch)
    ma = MonosaccharidePatch(saccharide_shape=(a,))
    return ma


def line_to(ax, sx, sy, ex, ey, zorder=1, color=None, lw=1):
    if color is None:
        color = default_line_color
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


def line_along(ax, points, zorder=1, color=None, lw=1):
    if color is None:
        color = default_line_color
    vertices = list(points)
    vertices.append([0, 0])
    codes = [Path.MOVETO] + [Path.LINETO] * (len(points) - 1)
    codes.append(Path.STOP)
    path = Path(vertices, codes)
    patch = patches.PathPatch(path, edgecolor=color, facecolor='none', lw=lw, zorder=zorder)
    ax.add_patch(patch)
    return patch


def transform_around_coords(x, y, xorigin, yorigin, transform):
    xy2 = np.array((x, y))
    xy1 = np.array((xorigin, yorigin))
    delta = xy1 - xy2
    reloc = transform.transform(delta)
    new_pos = xy1 + reloc
    return new_pos


class SymbolicNomenclatureBase(object):
    default_line_color = 'black'

    def get_layout_transform(self, orientation, **kwargs):
        if orientation in {'h', 'horizontal'}:
            return NodeLayoutTransforms(Affine2D().rotate_deg_around(0, 0, 90),
                                        Affine2D().rotate_deg(-90).translate(0, 0.15),
                                        Affine2D())

        elif orientation in {'v', 'vertical'}:
            return NodeLayoutTransforms(Affine2D().rotate_deg_around(0, 0, 0),
                                        Affine2D().rotate_deg(180).translate(0.05, -0.05),
                                        Affine2D())
        else:
            raise ValueError("Unknown Layout %r" % (orientation,))

    def draw(self, monosaccharide, x, y, ax, tree_node=None, scale=0.1, annotation_transform=None, **kwargs):
        raise NotImplementedError()

    def _get_horizontal_alignment(self, kwargs):
        value = kwargs.get("ha")
        if value is not None:
            return value
        value = kwargs.get("horizontal_alignment")
        return value

    def _get_vertical_alignment(self, kwargs):
        value = kwargs.get("va")
        if value is not None:
            return value
        value = kwargs.get("vertical_alignment")
        return value

    def draw_text(self, ax, x, y, text, scale=0.1, center=False, **kwargs):
        kwargs.setdefault('lw', line_weight / 20.0)
        fs = kwargs.get("fontsize", 2) * scale * .25
        t_path = TextPath((x, y), s=text, size=fs)
        ha = self._get_horizontal_alignment(kwargs)
        if center or ha == 'center':
            cx, cy = centroid(t_path)
            tf = Affine2D()
            tf.translate(x - cx, y - cy)
            t_path = t_path.transformed(tf)
        va = self._get_vertical_alignment(kwargs)
        if va == 'bottom':
            bbox = bounding_box(t_path)
            ymin = bbox.y0
            cx, cy = centroid(t_path)
            tf = Affine2D()
            tf.translate(0, cy - ymin)
            t_path = t_path.transformed(tf)
        patch = patches.PathPatch(t_path, facecolor="black", lw=kwargs['lw'], zorder=4)
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

    def draw_substituents(self, ax, substituents, node_x, node_y, annotation_transform=None, **kwargs):
        subs = []
        sub_y = node_y - (0.15 * (len(substituents) - 1))
        sub_x = node_x - 0.45
        if annotation_transform is None:
            annotation_transform = Affine2D()
        for pos, subst_name in substituents:
            sub_x_tr, sub_y_tr = transform_around_coords(sub_x, sub_y, node_x, node_y, annotation_transform)
            sub_t = self.draw_text(
                ax, sub_x_tr, sub_y_tr, str(pos) + self.format_text(subst_name), center=True, fontsize=6)
            sub_y += 0.4
            subs.append(sub_t)
        return subs
