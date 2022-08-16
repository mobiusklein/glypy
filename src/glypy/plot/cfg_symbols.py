'''
An implementation of a subset of the following symbolic nomenclature:
    http://www.ncbi.nlm.nih.gov/books/NBK310273/table/symbolnomenclature.T.monosaccharide_symb/?report=objectonly
'''

import logging
from collections import Counter
from functools import partial

import numpy as np
from matplotlib.path import Path
from matplotlib.textpath import TextPath
import matplotlib.patches as patches
import matplotlib
from matplotlib.colors import rgb2hex
from matplotlib.transforms import Affine2D

import glypy
from glypy.structure import Modification, Stem
from glypy.utils.enum import Enum
from glypy.io.nomenclature import identity

from .common import MonosaccharidePatch
from .symbolic_nomenclature import SymbolicNomenclatureBase

logger = logging.getLogger(__name__)

default_line_color = 'black'
default_scale_factor = 1.0
zorder = 2
line_weight = 0.5


NeuAc = glypy.monosaccharides.NeuAc
NeuGc = glypy.monosaccharides.NeuGc


class CFGNomenclature(SymbolicNomenclatureBase):

    class ResidueColor(Enum):
        gal = rgb2hex((255 / 255., 255 / 255., 0 / 255.))  # yellow
        lyx = rgb2hex((255 / 255., 255 / 255., 0 / 255.))  # yellow
        glc = rgb2hex((0 / 255., 0 / 255., 255 / 255.))  # blue
        man = rgb2hex((0, 200 / 255., 50 / 255.))  # green
        fuc = rgb2hex((255 / 255., 0 / 255., 0 / 255.))  # red
        xyl = rgb2hex((250 / 255., 234 / 255., 213 / 255.))  # orange
        neuac = rgb2hex((200 / 255., 0 / 255., 200 / 255.))  # purple
        neugc = rgb2hex((233 / 255., 255 / 255., 255 / 255.))  # light blue
        kdn = rgb2hex((0, 200 / 255., 50 / 255.))  # green
        glca = rgb2hex((0 / 255., 0 / 255., 255 / 255.))  # blue
        idoa = rgb2hex((150 / 255., 100 / 255., 50 / 255.))  # tan
        gala = rgb2hex((255 / 255., 255 / 255., 0 / 255.))  # yellow
        mana = rgb2hex((0, 200 / 255., 50 / 255.))  # green
        generic = 'white'

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

    def residue_color(self, monosaccharide):
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
            return self.resolve_acid_color(monosaccharide)
        if "hex" in [monosaccharide.superclass]:
            if any(mod == Modification.d for p, mod in monosaccharide.modifications.items()) and\
                    monosaccharide.stem == (Stem.gal,):
                return self.ResidueColor.fuc
        try:
            return self.ResidueColor[monosaccharide.stem[0]]
        except KeyError:
            return self.ResidueColor.generic

    def resolve_acid_color(self, monosaccharide):
        '''
        Resolve the special case in :func:`residue_color` for acidic residues
        '''
        if ('gro' in monosaccharide.stem) and ('gal' in monosaccharide.stem):
            if any(sub.name == 'n_acetyl' for p, sub in monosaccharide.substituents()):
                return self.ResidueColor.neuac
            elif any(sub.name == 'n_glycolyl' for p, sub in monosaccharide.substituents()):
                return self.ResidueColor.neugc
            else:
                return self.ResidueColor.kdn
        elif 'glc' in monosaccharide.stem:
            return self.ResidueColor.glca
        elif 'gal' in monosaccharide.stem:
            return self.ResidueColor.gala
        elif 'man' in monosaccharide.stem:
            return self.ResidueColor.mana
        elif 'ido' in monosaccharide.stem:
            return self.ResidueColor.idoa

    def residue_shape(self, monosaccharide):
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
            return self.resolve_acid_shape(monosaccharide)
        if "hex" in [monosaccharide.superclass]:
            if any(sub.name == 'n_acetyl' for p, sub in monosaccharide.substituents()):
                return self.ResidueShape.square
            elif any(sub.name == 'amino' or sub.name.startswith("n_") for p, sub in monosaccharide.substituents()):
                return self.ResidueShape.bisected_square
            elif any(mod == Modification.d for p, mod in monosaccharide.modifications.items()):
                return self.ResidueShape.triangle

            return self.ResidueShape.circle
        elif "pen" in [monosaccharide.superclass]:
            if 'xyl' in monosaccharide.stem or 'lyx' in monosaccharide.stem:
                return self.ResidueShape.star
        return self.ResidueShape.generic

    def resolve_acid_shape(self, monosaccharide):
        '''
        Resolve the special case in :func:`residue_shape` for acidic residues
        '''
        if ('gro' in monosaccharide.stem) and ('gal' in monosaccharide.stem):
            if any(sub.name == 'n_acetyl' for p, sub in monosaccharide.substituents()):
                return self.ResidueShape.diamond
            elif any(sub.name == 'n_glycolyl' for p, sub in monosaccharide.substituents()):
                return self.ResidueShape.diamond
            else:
                return self.ResidueShape.diamond
        elif 'glc' in monosaccharide.stem:
            return self.ResidueShape.top_bisected_diamond
        elif 'gal' in monosaccharide.stem:
            return self.ResidueShape.left_bisected_diamond
        elif 'man' in monosaccharide.stem:
            return self.ResidueShape.right_bisected_diamond
        elif 'ido' in monosaccharide.stem:
            return self.ResidueShape.bottom_bisected_diamond
        else:
            return self.ResidueShape.bottom_bisected_diamond

    def draw_circle(self, ax, x, y, color, scale=0.1):
        path = Path(Path.unit_circle().vertices * scale, Path.unit_circle().codes)
        trans = Affine2D().translate(x, y)
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
        trans = Affine2D().translate(x, y)
        t_path = path.transformed(trans)
        patch = patches.PathPatch(
            t_path, facecolor=color.value, lw=line_weight, zorder=2)
        a = ax.add_patch(patch)
        ma = MonosaccharidePatch(saccharide_shape=(a,))
        return ma

    def draw_triangle(self, ax, x, y, color, scale=0.1):
        path = Path(Path.unit_regular_polygon(3).vertices * scale, Path.unit_regular_polygon(3).codes)
        trans = Affine2D().translate(
            x, y).rotate_deg_around(x, y, 0)
        t_path = path.transformed(trans)
        patch = patches.PathPatch(
            t_path, facecolor=color.value, lw=line_weight, zorder=2)
        a = ax.add_patch(patch)
        ma = MonosaccharidePatch(saccharide_shape=(a,))
        return ma

    def draw_bisected_square(self, ax, x, y, color, scale=0.1):
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
            Affine2D().translate(x, y))
        upper_path = Path(upper_verts, codes).transformed(
            Affine2D().translate(x, y))

        patch = patches.PathPatch(
            lower_path, facecolor=color.value, lw=line_weight, zorder=2)
        a = ax.add_patch(patch)
        patch = patches.PathPatch(
            upper_path, facecolor="white", lw=line_weight, zorder=2)
        b = ax.add_patch(patch)
        ma = MonosaccharidePatch(saccharide_shape=(a, b))
        return ma

    def draw_diamond(self, ax, x, y, color, scale=0.1):
        path = Path(Path.unit_regular_polygon(4).vertices * scale, Path.unit_regular_polygon(4).codes)
        trans = Affine2D().translate(x, y)
        t_path = path.transformed(trans)
        patch = patches.PathPatch(
            t_path, facecolor=color.value, lw=line_weight, zorder=2)
        a = (ax.add_patch(patch),)
        ma = MonosaccharidePatch(saccharide_shape=(a,))
        return ma

    def draw_right_bisected_diamond(self, ax, x, y, color, scale=0.1):
        return self.draw_horizontal_bisected_diamond(ax, x, y, color, scale, 'left')

    def draw_left_bisected_diamond(self, ax, x, y, color, scale=0.1):
        return self.draw_horizontal_bisected_diamond(ax, x, y, color, scale, 'right')

    def draw_top_bisected_diamond(self, ax, x, y, color, scale=0.1):
        return self.draw_vertical_bisected_diamond(ax, x, y, color, scale, 'top')

    def draw_bottom_bisected_diamond(self, ax, x, y, color, scale=0.1):
        return self.draw_vertical_bisected_diamond(ax, x, y, color, scale, 'bottom')

    def draw_vertical_bisected_diamond(self, ax, x, y, color, scale=0.1, side=None):
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
            Affine2D().translate(x, y).rotate_deg_around(x, y, 45))
        upper_path = Path(upper_verts, codes).transformed(
            Affine2D().translate(x, y).rotate_deg_around(x, y, 45))

        try:
            color = color.value
        except AttributeError:
            color = 'white'

        if side == 'top':
            top_color = color
            bottom_color = 'white'
        elif side == 'bottom':
            top_color = 'white'
            bottom_color = color
        patch = patches.PathPatch(
            lower_path, facecolor=bottom_color, lw=line_weight, zorder=2)
        a = ax.add_patch(patch)
        patch = patches.PathPatch(
            upper_path, facecolor=top_color, lw=line_weight, zorder=2)
        b = ax.add_patch(patch)
        ma = MonosaccharidePatch(saccharide_shape=(a, b))
        return ma
        # return a, b

    def draw_horizontal_bisected_diamond(self, ax, x, y, color, scale=0.1, side=None):
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

        try:
            color = color.value
        except AttributeError:
            color = 'white'

        if side == 'left':
            left_color = color
            right_color = 'white'
        elif side == 'right':
            left_color = 'white'
            right_color = color
        left_path = Path(left_verts, codes).transformed(
            Affine2D().translate(x, y).rotate_deg_around(x, y, -45))
        right_path = Path(right_verts, codes).transformed(
            Affine2D().translate(x, y).rotate_deg_around(x, y, -45))

        patch = patches.PathPatch(
            left_path, facecolor=left_color, lw=line_weight, zorder=2)
        a = ax.add_patch(patch)
        patch = patches.PathPatch(
            right_path, facecolor=right_color, lw=line_weight, zorder=2)
        b = ax.add_patch(patch)
        ma = MonosaccharidePatch(
            saccharide_shape=(a, b), saccharide_label=None)
        return ma
        # return a, b

    def draw_star(self, ax, x, y, color, scale=0.1):
        path = Path(Path.unit_regular_star(5, 0.45).vertices * scale, Path.unit_regular_star(5, 0.45).codes)
        trans = Affine2D().translate(x, y)
        t_path = path.transformed(trans)
        patch = patches.PathPatch(
            t_path, facecolor=color.value, lw=line_weight, zorder=2)
        a = ax.add_patch(patch)
        ma = MonosaccharidePatch(saccharide_shape=(a,))
        return ma

    def draw_generic(self, ax, x, y, name, n_points=6, scale=0.1):
        unit_polygon = Path.unit_regular_polygon(n_points)
        path = Path(unit_polygon.vertices * scale, unit_polygon.codes)
        trans = Affine2D().translate(x, y)
        t_path = path.transformed(trans)
        name = TextPath((x - (0.35 * scale), y), s=name, size=2 * scale * .25)
        patch = patches.PathPatch(
            t_path, facecolor="white", lw=line_weight, zorder=2)
        a = ax.add_patch(patch)
        patch = patches.PathPatch(name, lw=line_weight, zorder=2)
        s = ax.add_patch(patch)
        ma = MonosaccharidePatch(saccharide_shape=(a,), saccharide_label=(s,))
        return ma

    def get_drawer(self, shape):
        shape_name = shape.name
        drawer = getattr(self, 'draw_%s' % (shape_name,), None)
        return drawer

    def resolve_generic_name(self, monosaccharide):
        try:
            abbrev = identity.identify(monosaccharide)
        except identity.IdentifyException:
            if monosaccharide.stem[0] == Stem.x:
                abbrev = monosaccharide.superclass.name.lower().capitalize()
            else:
                abbrev = ','.join(s.name.lower().capitalize()
                                  for s in monosaccharide.stem)
        return abbrev

    def get_relevant_substituents(self, residue, shape=None):
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

        shape = self.residue_shape(residue) if shape is None else shape
        if shape != self.ResidueShape.generic:
            substituents = list(sub.name for p, sub in residue.substituents() if not sub._derivatize)
            if shape == self.ResidueShape.square:
                substituents.pop(substituents.index("n_acetyl"))
            elif shape == self.ResidueShape.bisected_square:
                try:
                    index_amino = substituents.index("amino")
                    substituents.pop(index_amino)
                except ValueError:
                    pass
                try:
                    pass
                    # index_n_sulfate = substituents.index('n_sulfate')
                    # TODO:
                    # Find a way to forward different substituents
                    #
                    # substituents.pop(index_n_sulfate)
                    # substituents.append("sulfate")
                except ValueError:
                    pass
                return self._pruned_substituents(residue, substituents)
            elif shape == self.ResidueShape.diamond:
                color = self.residue_color(residue)
                if color == self.residue_color(NeuAc):
                    substituents.pop(substituents.index("n_acetyl"))
                elif color == self.residue_color(NeuGc):
                    substituents.pop(substituents.index("n_glycolyl"))
            return self._pruned_substituents(residue, substituents)

        return list((p, sub.name) for p, sub in residue.substituents() if not sub._derivatize)

    def _pruned_substituents(self, residue, substituents):
        relevant_substituents = Counter(substituents)
        results = []
        for p, sub in residue.substituents()[::-1]:
            if relevant_substituents[sub.name] > 0:
                results.append((p, sub.name))
                relevant_substituents[sub.name] -= 1
        return results

    def get_symbol(self, monosaccharide):
        '''
        Convenience function to retrieve the shape and color of `monosaccharide`
        '''
        shp = self.residue_shape(monosaccharide)
        col = self.residue_color(monosaccharide)
        return shp, col

    def draw(self, monosaccharide, x, y, ax, tree_node=None, scale=0.1, annotation_transform=None, **kwargs):
        '''
        Renders `monosaccharide` at the given `(x, y)` coordinates on the `matplotlib.Axis`
        `ax` provided. Determines the shape to use by :func:`residue_shape` and color by
        :func:`residue_color`. The shape value is used to select the specialized draw_* function
        '''
        abbrev = None
        shape, color = self.get_symbol(monosaccharide)
        if shape == self.ResidueShape.generic:
            abbrev = self.resolve_generic_name(monosaccharide)
        drawer = self.get_drawer(shape)
        if drawer is None:
            raise Exception("Don't know how to draw {}".format(
                (shape, monosaccharide)))

        res = None
        if shape == self.ResidueShape.generic:
            res = drawer(
                ax, x, y, abbrev, n_points=monosaccharide.superclass.value or 1, scale=scale)
        else:
            res = drawer(ax, x, y, color, scale=scale)
        substituents = self.get_relevant_substituents(monosaccharide)

        # Render substituents along the bottom of the monosaccharide
        # These layouts should be moved to be defined by the DrawTreeNode
        if annotation_transform is None:
            annotation_transform = Affine2D()
        node_x, node_y = res.centroid()
        subs = self.draw_substituents(ax, substituents, node_x, node_y, annotation_transform, **kwargs)
        res.add_substituents(subs)
        return res

    def __call__(self, *args, **kwargs):
        return self.draw(*args, **kwargs)
