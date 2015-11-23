import logging
from functools import partial
from collections import Counter

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.textpath import TextPath
import matplotlib.patches as patches
from matplotlib.transforms import IdentityTransform
import matplotlib
from matplotlib.colors import rgb2hex

from glypy.structure import Modification, Stem, SuperClass
from glypy.utils import where
from glypy.utils.enum import Enum
from glypy.io.nomenclature import identity
from glypy.algorithms import similarity

from .common import MonosaccharidePatch
from .geometry import centroid

logger = logging.getLogger(__name__)

default_line_color = 'black'
default_scale_factor = 1.0
zorder = 2
line_weight = 0.5


class ResidueShape(Enum):
    circle = 1
    square = 2
    bisected_square = 3
    triangle = 4
    star = 5
    diamond = 6
    top_bisected_diamond = 7
    left_bisected_diamond = 8
    right_bisected_diamond = 9
    bottom_bisected_diamond = 10
    generic = 11


class UnknownShapeException(Exception):
    pass


def _pruned_substituents(residue, substituents):

    relevant_substituents = Counter(substituents)
    results = []
    for p, sub in residue.substituents():
        if relevant_substituents[sub.name] > 0:
            results.append((p, sub.name))
            relevant_substituents[sub.name] -= 1
    return results


def get_relevant_substituents(residue, shape=None):
    '''
    Given the shape for a residue, determine which of its substituents must
    be explicitly drawn.  Calls :func:`residue_shape` if `shape` is not
    provided.

    Parameters
    ----------
    residue: |Monosaccharide|
        The monosaccharide residue being rendered
    shape: ResidueShape or |None|
        The shape enum being used to represent `residue`. Defaults to None.
        If `shape` is |None|, it is calculated by :func:`residue_shape`.

    Returns
    -------
    |list| of |Substituent|s
    '''

    shape = residue_shape(residue) if shape is None else shape
    if shape != ResidueShape.generic:
        substituents = list(sub.name for p, sub in residue.substituents() if not sub._derivatize)
        if shape == ResidueShape.square:
            substituents.pop(substituents.index("n_acetyl"))
        elif shape == ResidueShape.bisected_square:
            try:
                index_amino = substituents.index("amino")
                substituents.pop(index_amino)
            except ValueError:
                pass
            try:
                index_n_sulfate = substituents.index('n_sulfate')
                # TODO:
                # Find a way to forward different substituents
                #
                # substituents.pop(index_n_sulfate)
                # substituents.append("sulfate")
            except ValueError:
                pass
            return _pruned_substituents(residue, substituents)
        elif shape == ResidueShape.diamond:
            color = residue_color(residue)
            if color == ResidueColor.neuac:
                substituents.pop(substituents.index("n_acetyl"))
            elif color == ResidueColor.neugc:
                substituents.pop(substituents.index("n_glycolyl"))
        return _pruned_substituents(residue, substituents)

    return list((p, sub.name) for p, sub in residue.substituents() if not sub._derivatize)


def residue_shape(monosaccharide):
    '''
    Determine which shape to use to represent `monosaccharide` under the CFG
    symbol nomenclature.

    Parameters
    ----------
    monosaccharide: |Monosaccharide|
        The residue to be rendered

    Returns
    -------
    ResidueShape.EnumValue

    '''
    if any(mod == Modification.a for p, mod in monosaccharide.modifications.items()):
        return resolve_acid_shape(monosaccharide)
    if "hex" in [monosaccharide.superclass]:
        if any(sub.name == 'n_acetyl' for p, sub in monosaccharide.substituents()):
            return ResidueShape.square
        elif any(sub.name == 'amino' or sub.name.startswith("n_") for p, sub in monosaccharide.substituents()):
            return ResidueShape.bisected_square
        elif any(mod == Modification.d for p, mod in monosaccharide.modifications.items()):
            return ResidueShape.triangle

        return ResidueShape.circle
    elif "pen" in [monosaccharide.superclass]:
        if 'xyl' in monosaccharide.stem:
            return ResidueShape.star
    return ResidueShape.generic


def resolve_acid_shape(monosaccharide):
    '''
    Resolve the special case in :func:`residue_shape` for acidic residues
    '''
    if ('gro' in monosaccharide.stem) and ('gal' in monosaccharide.stem):
        if any(sub.name == 'n_acetyl' for p, sub in monosaccharide.substituents()):
            return ResidueShape.diamond
        elif any(sub.name == 'n_glycolyl' for p, sub in monosaccharide.substituents()):
            return ResidueShape.diamond
        else:
            return ResidueShape.diamond
    elif 'glc' in monosaccharide.stem:
        return ResidueShape.top_bisected_diamond
    elif 'gal' in monosaccharide.stem:
        return ResidueShape.left_bisected_diamond
    elif 'man' in monosaccharide.stem:
        return ResidueShape.right_bisected_diamond
    elif 'ido' in monosaccharide.stem:
        return ResidueShape.bottom_bisected_diamond


class ResidueColor(Enum):
    gal = rgb2hex((255/255., 255/255., 0/255.))  # yellow
    glc = rgb2hex((0/255., 0/255., 255/255.))  # blue
    man = rgb2hex((0, 200/255., 50/255.))  # green
    fuc = rgb2hex((255/255., 0/255., 0/255.))  # red
    xyl = rgb2hex((250/255., 234/255., 213/255.))  # orange
    neuac = rgb2hex((200/255., 0/255., 200/255.))  # purple
    neugc = rgb2hex((233/255., 255/255., 255/255.))  # light blue
    kdn = rgb2hex((0, 200/255., 50/255.))  # green
    glca = rgb2hex((0/255., 0/255., 255/255.))  # blue
    idoa = rgb2hex((150/255., 100/255., 50/255.))  # tan
    gala = rgb2hex((255/255., 255/255., 0/255.))  # yellow
    mana = rgb2hex((0, 200/255., 50/255.))  # green
    generic = 'white'


def residue_color(monosaccharide):
    '''
    Determine which color to use to represent `monosaccharide` under the CFG
    symbol nomenclature.

    Parameters
    ----------
    monosaccharide: |Monosaccharide|
        The residue to be rendered

    Returns
    -------
    ResidueColor.EnumValue

    '''
    if any(mod == Modification.a for p, mod in monosaccharide.modifications.items()):
            return resolve_acid_color(monosaccharide)
    if "hex" in [monosaccharide.superclass]:
        if any(mod == Modification.d for p, mod in monosaccharide.modifications.items()) and\
         monosaccharide.stem == (Stem.gal,):
            return ResidueColor.fuc
    try:
        return ResidueColor[monosaccharide.stem[0]]
    except KeyError:
        return ResidueColor.generic


def resolve_acid_color(monosaccharide):
    '''
    Resolve the special case in :func:`residue_color` for acidic residues
    '''
    if ('gro' in monosaccharide.stem) and ('gal' in monosaccharide.stem):
        if any(sub.name == 'n_acetyl' for p, sub in monosaccharide.substituents()):
            return ResidueColor.neuac
        elif any(sub.name == 'n_glycolyl' for p, sub in monosaccharide.substituents()):
            return ResidueColor.neugc
        else:
            return ResidueColor.kdn
    elif 'glc' in monosaccharide.stem:
        return ResidueColor.glca
    elif 'gal' in monosaccharide.stem:
        return ResidueColor.gala
    elif 'man' in monosaccharide.stem:
        return ResidueColor.mana
    elif 'ido' in monosaccharide.stem:
        return ResidueColor.idoa


def get_symbol(monosaccharide):
    '''
    Convenience function to retrieve the shape and color of `monosaccharide`
    '''
    shp = residue_shape(monosaccharide)
    col = residue_color(monosaccharide)
    return shp, col


def resolve_generic_name(monosaccharide):
    try:
        abbrev = identity.identify(monosaccharide)
    except identity.IdentifyException:
        if monosaccharide.stem[0] == Stem.x:
            abbrev = monosaccharide.superclass.name.lower().capitalize()
        else:
            abbrev = ','.join(s.name.lower().capitalize() for s in monosaccharide.stem)
    return abbrev


def draw(monosaccharide, x, y, ax, tree_node=None, scale=0.1, **kwargs):
    '''
    Renders `monosaccharide` at the given `(x, y)` coordinates on the `matplotlib.Axis`
    `ax` provided. Determines the shape to use by :func:`residue_shape` and color by
    :func:`residue_color`. The shape value is used to select the specialized draw_* function
    '''
    abbrev = None
    shape, color = get_symbol(monosaccharide)
    if shape == ResidueShape.generic:
        abbrev = resolve_generic_name(monosaccharide)
    drawer = draw_map.get(shape)
    if drawer is None:
        raise Exception("Don't know how to draw {}".format((shape, monosaccharide)))

    res = None
    if shape == ResidueShape.generic:
        res = drawer(ax, x, y, abbrev, n_points=monosaccharide.superclass.value or 1, scale=scale)
    else:
        res = drawer(ax, x, y, color, scale=scale)
    substituents = get_relevant_substituents(monosaccharide)

    # Render substituents along the bottom of the monosaccharide
    # These layouts should be moved to be defined by the DrawTreeNode
    subs = []
    sub_y = y - (0.35 * (len(substituents) - 1))
    sub_x = x - 0.45
    for pos, subst_name in substituents:
        sub_t = draw_text(ax, sub_x, sub_y, str(pos) + format_text(subst_name), fontsize=6)
        sub_y += 0.3
        subs.append(sub_t)
    res.add_substituents(subs)
    return res

draw_map = {}


def draw_circle(ax, x, y, color, scale=0.1):
    path = Path(unit_circle.vertices * scale, unit_circle.codes)
    trans = matplotlib.transforms.Affine2D().translate(x, y)
    t_path = path.transformed(trans)
    patch = patches.PathPatch(t_path, facecolor=color.value, lw=line_weight, zorder=2)
    a = ax.add_patch(patch)
    ma = MonosaccharidePatch(saccharide_shape=(a,))
    return ma
    # return (a,)
draw_map[ResidueShape.circle] = draw_circle
unit_circle = Path.unit_circle()


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
    patch = patches.PathPatch(t_path, facecolor=color.value, lw=line_weight, zorder=2)
    a = ax.add_patch(patch)
    ma = MonosaccharidePatch(saccharide_shape=(a,))
    return ma
    # return (a,)
draw_map[ResidueShape.square] = draw_square
unit_rectangle = Path.unit_rectangle()


def draw_triangle(ax, x, y, color, scale=0.1):
    path = Path(unit_triangle.vertices * scale, unit_triangle.codes)
    trans = matplotlib.transforms.Affine2D().translate(x, y).rotate_deg_around(x, y, -90)
    t_path = path.transformed(trans)
    patch = patches.PathPatch(t_path, facecolor=color.value, lw=line_weight, zorder=2)
    a = ax.add_patch(patch)
    ma = MonosaccharidePatch(saccharide_shape=(a,))
    return ma
    # return (a,)
draw_map[ResidueShape.triangle] = draw_triangle
unit_triangle = Path.unit_regular_polygon(3)


def draw_bisected_square(ax, x, y, color, scale=0.1):
    lower_verts = (np.array([
            (0., 0.),
            (1.0, 0),
            (0, 1.0),
            (0, 0),
            (0., 0.),
            ]) - 0.5) / 5

    upper_verts = (np.array([
            (1., 1.),
            (1.0, 0),
            (0, 1.0),
            (1, 1),
            (0., 0.),
            ]) - 0.5) / 5

    codes = [Path.MOVETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.CLOSEPOLY,
             ]

    lower_path = Path(lower_verts, codes).transformed(
        matplotlib.transforms.Affine2D().translate(x, y))
    upper_path = Path(upper_verts, codes).transformed(
        matplotlib.transforms.Affine2D().translate(x, y))

    patch = patches.PathPatch(lower_path, facecolor=color.value, lw=line_weight, zorder=2)
    a = ax.add_patch(patch)
    patch = patches.PathPatch(upper_path, facecolor="white", lw=line_weight, zorder=2)
    b = ax.add_patch(patch)
    ma = MonosaccharidePatch(saccharide_shape=(a, b))
    return ma
    # return a, b
draw_map[ResidueShape.bisected_square] = draw_bisected_square


def draw_diamond(ax, x, y, color, scale=0.1):
    path = Path(unit_diamond.vertices * scale, unit_diamond.codes)
    trans = matplotlib.transforms.Affine2D().translate(x, y)
    t_path = path.transformed(trans)
    patch = patches.PathPatch(t_path, facecolor=color.value, lw=line_weight, zorder=2)
    a = (ax.add_patch(patch),)
    ma = MonosaccharidePatch(saccharide_shape=(a,))
    return ma
    # return (a,)
draw_map[ResidueShape.diamond] = draw_diamond
unit_diamond = Path.unit_regular_polygon(4)


def draw_vertical_bisected_diamond(ax, x, y, color, scale=0.1, side=None):
    lower_verts = (np.array([
            (0., 0.),
            (1.0, 0),
            (0, 1.0),
            (0, 0),
            (0., 0.),
            ]) - 0.5) / 5

    upper_verts = (np.array([
            (1., 1.),
            (1.0, 0),
            (0, 1.0),
            (1, 1),
            (0., 0.),
            ]) - 0.5) / 5

    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]

    lower_path = Path(lower_verts, codes).transformed(
        matplotlib.transforms.Affine2D().translate(x, y).rotate_deg_around(x, y, 45))
    upper_path = Path(upper_verts, codes).transformed(
        matplotlib.transforms.Affine2D().translate(x, y).rotate_deg_around(x, y, 45))

    if side == 'top':
        top_color = color.value
        bottom_color = 'white'
    elif side == 'bottom':
        top_color = 'white'
        bottom_color = color.value
    patch = patches.PathPatch(lower_path, facecolor=bottom_color, lw=line_weight, zorder=2)
    a = ax.add_patch(patch)
    patch = patches.PathPatch(upper_path, facecolor=top_color, lw=line_weight, zorder=2)
    b = ax.add_patch(patch)
    ma = MonosaccharidePatch(saccharide_shape=(a, b))
    return ma
    # return a, b
draw_map[ResidueShape.top_bisected_diamond] = partial(draw_vertical_bisected_diamond, side='top')
draw_map[ResidueShape.bottom_bisected_diamond] = partial(draw_vertical_bisected_diamond, side='bottom')


def draw_horizontal_bisected_diamond(ax, x, y, color, scale=0.1, side=None):
    left_verts = (np.array([
            (0., 0.),
            (1.0, 0),
            (0, 1.0),
            (0, 0),
            (0., 0.),
            ]) - 0.5) / 5

    right_verts = (np.array([
            (1., 1.),
            (1.0, 0),
            (0, 1.0),
            (1, 1),
            (0., 0.),
            ]) - 0.5) / 5

    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]
    if side == 'left':
        left_color = color.value
        right_color = 'white'
    elif side == 'right':
        left_color = 'white'
        right_color = color.value
    left_path = Path(left_verts, codes).transformed(
        matplotlib.transforms.Affine2D().translate(x, y).rotate_deg_around(x, y, -45))
    right_path = Path(right_verts, codes).transformed(
        matplotlib.transforms.Affine2D().translate(x, y).rotate_deg_around(x, y, -45))

    patch = patches.PathPatch(left_path, facecolor=left_color, lw=line_weight, zorder=2)
    a = ax.add_patch(patch)
    patch = patches.PathPatch(right_path, facecolor=right_color, lw=line_weight, zorder=2)
    b = ax.add_patch(patch)
    ma = MonosaccharidePatch(saccharide_shape=(a, b), saccharide_label=None)
    return ma
    # return a, b

draw_map[ResidueShape.right_bisected_diamond] = partial(draw_horizontal_bisected_diamond, side='right')
draw_map[ResidueShape.left_bisected_diamond] = partial(draw_horizontal_bisected_diamond, side='left')


def draw_star(ax, x, y, color, scale=0.1):
    path = Path(unit_star.vertices * scale, unit_star.codes)
    trans = matplotlib.transforms.Affine2D().translate(x, y)
    t_path = path.transformed(trans)
    patch = patches.PathPatch(t_path, facecolor=color.value, lw=line_weight, zorder=2)
    a = ax.add_patch(patch)
    ma = MonosaccharidePatch(saccharide_shape=(a,))
    return ma
    # return (a,)
unit_star = Path.unit_regular_star(5, 0.3)
draw_map[ResidueShape.star] = draw_star


def draw_generic(ax, x, y, name, n_points=6, scale=0.1):
    unit_polygon = Path.unit_regular_polygon(n_points)
    path = Path(unit_polygon.vertices * scale, unit_polygon.codes)
    trans = matplotlib.transforms.Affine2D().translate(x, y)
    t_path = path.transformed(trans)
    name = TextPath((x - (0.35 * scale), y), s=name, size=2 * scale * .25)
    patch = patches.PathPatch(t_path, facecolor="white", lw=line_weight, zorder=2)
    a = ax.add_patch(patch)
    patch = patches.PathPatch(name, lw=line_weight, zorder=2)
    s = ax.add_patch(patch)
    ma = MonosaccharidePatch(saccharide_shape=(a,), saccharide_label=(s,))
    return ma
    # return (a,)
draw_map[ResidueShape.generic] = draw_generic


def draw_text(ax, x, y, text, scale=0.1, **kwargs):
    fs = kwargs.get("fontsize", 2) * scale * .25
    t_path = TextPath((x, y), s=text, size=fs)
    patch = patches.PathPatch(t_path, facecolor="black", lw=line_weight / 20., zorder=4)
    a = ax.add_patch(patch)
    return (a,)


def format_text(text):
    label = ''.join([token.capitalize()[:2] for token in text.split('_', 1)])
    return label


def line_to(ax, sx, sy, ex, ey, zorder=1, color='black'):
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
    patch = patches.PathPatch(path)
    ax.add_patch(patch)
    return patch
