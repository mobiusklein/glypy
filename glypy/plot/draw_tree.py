import re
import operator
from collections import defaultdict
from uuid import uuid4
from itertools import chain
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms

from six import string_types as basestring

import numpy as np
from glypy import monosaccharides
from glypy.structure import constants
from glypy.utils import root
from glypy.io.nomenclature import identity
from glypy.algorithms.similarity import is_derivatized
from glypy.composition.composition_transform import strip_derivatization

from .buchheim import BalancedTreeLayout
from .topological_layout import TopologicalTreeLayout
from . import cfg_symbols, iupac_symbols, snfg_symbols
from .geometry import centroid, interpolate
from .fragment_annotation import BondCleavageArtist


line_to = cfg_symbols.CFGNomenclature().line_to


nomenclature_map = {
    "cfg": cfg_symbols.CFGNomenclature(),
    'iupac': iupac_symbols.IUPACTextSymbolicNomenclature(),
    "snfg": snfg_symbols.SNFGNomenclature()
}

layout_map = {
    "balanced": BalancedTreeLayout,
    "topological": TopologicalTreeLayout
}

DEFAULT_SYMBOL_SCALE_FACTOR = 0.17


#: :data:`special_cases` contains a list of names for
#: special case monosaccharides
_special_case_fucose = monosaccharides.Fuc
_special_case_fucose.configuration = None
_special_case_fucose.anomer = None
_special_case_xylose = monosaccharides.Xyl
_special_case_xylose.configuration = None
_special_case_xylose.anomer = None

special_cases = [_special_case_fucose, _special_case_xylose]


anomer_symbol_map = {
    constants.Anomer.alpha: r'\alpha',
    constants.Anomer.beta: r'\beta',
    constants.Anomer.uncyclized: r'o',
    constants.Anomer.x: r'?'
}


def sign(x):
    if x > 0:
        return 1
    else:
        return -1


def flatten(iterable):
    acc = []
    for item in iterable:
        try:
            acc.extend(flatten(item))
        except Exception:
            acc.append(item)
    return acc


def is_special_case(node):
    '''
    Check to see if `node` is a special case which requires a different layout
    scheme. See `special_cases`
    '''
    if is_derivatized(node):
        node = strip_derivatization(node.clone())
    for case in special_cases:
        if identity.is_a(node, case) and len(list(node.children())) == 0:
            return True
    return False


def box(node, ax=None):  # pragma: no cover
    if ax is None:
        ax = node.axes
    xmin, xmax, ymin, ymax = node.extrema()
    ax.plot((xmin - 0.2, xmax + 0.2), (ymin - 0.2, ymin - 0.2), c='b')
    ax.plot((xmin - 0.2, xmax + 0.2), (ymax + 0.2, ymax + 0.2), c='b')
    ax.plot((xmin - 0.2, xmin - 0.2), (ymin - 0.2, ymax + 0.2), c='b')
    ax.plot((xmax + 0.2, xmax + 0.2), (ymin - 0.2, ymax + 0.2), c='b')


def make_gid(s):  # pragma: no cover
    '''Formats a string to use ``kebab-case`` and sanitize illegal characters'''
    s = str(s)
    return re.sub(r'[\s\|\(\):]', '-', s)


def breadth_first_traversal(tree, visited=None):
    if visited is None:
        visited = set()
    if tree.id not in visited:
        visited.add(tree.id)
        yield tree
        for child in tree:
            for descend in breadth_first_traversal(child, visited=visited):
                yield descend


def _bounding_box_to_limits(dtree):
    xlim, ylim = dtree._compute_bounding_box()
    xlim = list(xlim)
    xlim[0] -= 2
    xlim[1] += 2
    ylim = list(ylim)
    ylim[0] -= 2
    ylim[1] += 2
    return np.array(xlim), np.array(ylim)


def enumerate_tree(tree, ax, labeler=None):
    """Label each node of `tree` on `ax` with the node's :attr:`id`

    Parameters
    ----------
    tree: :class:`DrawTreeNode`
        A drawn :class:`DrawTreeNode` object
    ax: matplotlib.Axes
        The axes on which `tree` has been drawn

    Returns
    -------
    tree, ax
    """
    if labeler is None:
        labeler = operator.attrgetter('id')
    nodes = list(breadth_first_traversal(tree))
    for node in nodes:
        x, y = node.coords()
        ax.text(x, y + 0.2, labeler(node.tree), color='red')
    # Update the plot
    ax.get_figure().canvas.draw()
    return tree, ax


def find_link(parent, child):  # pragma: no cover
    for p, link in parent.links.items():
        if child.parent_linkage == p:
            if link.child == child:
                return link


def get(tree, node_id):
    for node in breadth_first_traversal(tree):
        if node.tree.id == node_id:
            return node


def get_link_pair(tree, link_id):
    link = None
    for node in breadth_first_traversal(tree):
        for p, t_link in node.tree.links.items():
            if t_link.id == link_id:
                link = t_link
                break
    if link is None:
        raise Exception("Link not found")
    parent = get(tree, link.parent.id)
    child = get(tree, link.child.id)
    return parent, child


INF = float("inf")


class DrawTreeNode(object):
    def __init__(self, tree, parent=None, depth=0, number=1, parent_linkage=None, visited=None):
        if visited is None:
            visited = dict()
        self.id = tree.id
        if self.id in visited:
            raise Exception("Cycle detected")
        visited[self.id] = self
        self.x = -1.
        self.y = depth
        self.depth = depth
        self.tree = tree
        self.parent = parent
        self.parent_linkage = parent_linkage

        self._axes = None

        self.mask_special_cases = True
        # Unused. Artefact of Buccheim
        self.thread = None
        self.mod = 0
        self.ancestor = self
        self.change = self.shift = 0
        self._lmost_sibling = None

        # Number of the node in its group of siblings 1..n
        self.number = number

        if self.parent is not None:
            self.data = self.parent.data
            self.uuid = self.parent.uuid
        else:
            self.data = defaultdict(lambda: defaultdict(dict))
            self.uuid = uuid4().hex

        self.children = [self.resolve_creation_cycle(c[1], depth + 1, i + 1, c[0], visited)
                         for i, c in enumerate(sorted(tree.children(), key=lambda x: x[0]))]

        # A node owns all lines originating from it.
        self.lines = []

    def resolve_creation_cycle(self, node, depth, i, parent_linkage, visited):
        if node.id in visited:
            return visited[node.id]
        else:
            return self.__class__(node, self, depth + 1, i + 1, parent_linkage, visited)

    def __iter__(self):
        return iter(self.children)

    def traverse(self):
        return breadth_first_traversal(self)

    @property
    def children(self):
        if self.mask_special_cases:
            return [child for child in self._children if not (
                is_special_case(child.tree) and len(child.children) == 0)]
        else:
            return self._children

    @children.setter
    def children(self, value):
        self._children = value

    @property
    def axes(self):
        return self._axes

    @axes.setter
    def axes(self, value):
        self._axes = value
        for child in self:
            if child.axes != value:
                child.axes = value

    def fix_special_cases(self, offset=0.5, visited=None):
        if visited is None:
            visited = set()
        if self.id in visited:
            return
        visited.add(self.id)
        self.mask_special_cases = False
        for child in self.children:
            if is_special_case(child.tree):
                if child.parent_linkage > 3 or child.parent_linkage == -1:
                    child.y = self.y
                    child.x = self.x + offset
                elif child.parent_linkage > 1:
                    child.y = self.y
                    child.x = self.x - offset
                else:
                    raise Exception(
                        "Don't know how to handle special case child node %r, %r." % (
                            child, child.parent_linkage))
            child.fix_special_cases(offset, visited)

    def draw(self, at=(0, 0), ax=None, symbol_nomenclature="snfg", label=True, **kwargs):
        if isinstance(symbol_nomenclature, basestring):
            symbol_nomenclature = nomenclature_map[symbol_nomenclature]
        elif symbol_nomenclature is None:
            symbol_nomenclature = snfg_symbols.SNFGNomenclature()
        if ax is None:  # pragma: no cover
            ax = self.axes
            if ax is None:
                fig, ax = plt.subplots()
        self.draw_branches(at=at, ax=ax, symbol_nomenclature=symbol_nomenclature,
                           label=label, visited=set(), **kwargs)
        self.draw_nodes(at=at, ax=ax, symbol_nomenclature=symbol_nomenclature,
                        label=label, visited=set(), **kwargs)

    def coords(self, at=None):
        if at is None:
            at = (0, 0)
        if self.transform:
            at = self.transform.transform(at)
            x, y = self.transform.transform((self.x, self.y))
        else:
            x, y = self.x, self.y
        return x + at[0], y + at[1]

    def draw_branches(self, at=(0, 0), ax=None, symbol_nomenclature=None,
                      label=True, visited=None, **kwargs):
        '''
        Draw the edges linking parent nodes to child nodes. Also draws the parent outlink position
        and the child anomer symbol
        Parameters
        ----------
        at: :class:`tuple`
            The x, y coordinates at which to draw the glycan's reducing_end. Defaults to `(0, 0)`
        ax: :class:`matplotlib.axes.Axis`
            A matplotlib axis on which to draw. Useful for adding glycans to existing figures. If
            an axis is not provided, a new figure will be created and used.
        symbol_nomenclature: |str|
            A string mapping to the symbol nomenclature to use. Defaults to |None|. When `symbol_nomenclature`
            is |None|, `snfg` will be used.
        visited: |set|
            A |set| instance containing the node ids already visited while drawing. If the current
            node was already visited, its id will already be present. Used to avoid infinite cycles.
        '''
        if visited is None:
            visited = set()
        x, y = self.coords(at)
        if self.id in visited:
            return x, y
        visited.add(self.id)
        edge_weight = kwargs.get("edge_weight", 0.5)
        for child in self:
            cx, cy = child.draw_branches(at=at, ax=ax,
                                         symbol_nomenclature=symbol_nomenclature,
                                         label=label, visited=visited, **kwargs)
            patch = symbol_nomenclature.line_to(
                ax, x, y, cx, cy,
                color=kwargs.get('link_color', symbol_nomenclature.default_line_color),
                zorder=1, lw=edge_weight)
            self.data['lines'][self.id, child.id] = [patch]
            self.lines.append(patch)
        return x, y

    def draw_nodes(self, at=(0, 0), ax=None, symbol_nomenclature=None, label=True, visited=None, **kwargs):
        '''
        Draw the symbol representing the individual monosaccharide and its substituents.

        Parameters
        ----------
        at: :class:`tuple`
            The x, y coordinates at which to draw the glycan's reducing_end. Defaults to `(0, 0)`
        ax: :class:`matplotlib.axes.Axis`
            A matplotlib axis on which to draw. Useful for adding glycans to existing figures. If
            an axis is not provided, a new figure will be created and used.
        symbol_nomenclature: |str|
            A string mapping to the symbol nomenclature to use. Defaults to |None|. When `symbol_nomenclature`
            is |None|, `cfg` will be used.
        visited: |set|
            A |set| instance containing the node ids already visited while drawing. If the current
            node was already visited, its id will already be present. Used to avoid infinite cycles.
        '''
        if visited is None:
            visited = set()
        if self.id in visited:
            return
        visited.add(self.id)
        x, y = self.coords(at)
        patch_node = symbol_nomenclature.draw(self.tree, x, y, ax, tree_node=self,
                                              **kwargs)
        self.patch_node = patch_node
        residue_elements, substituent_elements, saccharide_label = patch_node
        self.data["patches"][self.tree.id] = [residue_elements]
        self.substituents_text = self.data["text"]["substituents"][self.id] = [substituent_elements]
        self.saccharide_label = self.data['text']['saccharide_label'][self.id] = [saccharide_label]
        self.data["position"][self.tree.id] = x, y

        patch_node.set_gid(self.uuid, self.id)

        for child in self:
            child.draw_nodes(at=at, ax=ax, symbol_nomenclature=symbol_nomenclature,
                             visited=visited, **kwargs)

    def _draw_linkage_label(self, child, at=(0, 0), ax=None, symbol_nomenclature=None, **kwargs):
        sx, sy = self.coords(at)
        cx, cy = child.coords(at)

        fontsize = kwargs.get("fontsize", 6)
        style = kwargs.get('label_style', 'split')

        position_x = (sx * 0.65 + cx * 0.35)
        if abs(sx - cx) < 1e-2:
            position_y = (sy + cy) / 2
        else:
            position_y = interpolate(sx, sy, cx, cy, position_x)

        position_num = -1
        for pos, link in self.tree.links.items():
            if link.is_child(child.tree):
                position_num = pos
                break

        anomer_x = (sx * 0.35 + cx * 0.65)
        if abs(sx - cx) < 1e-2:
            anomer_y = (sy + cy) / 2
        else:
            anomer_y = interpolate(sx, sy, cx, cy, anomer_x)

        if abs(anomer_x - position_x) < 1e-2 and abs(anomer_y - position_y) < 1e-2:
            anomer_x -= 0.07
            position_x += 0.05
            anomer_y -= 0.03
            position_y -= 0.03

        if style == 'split':
            position_text = symbol_nomenclature.draw_text(
                ax=ax, x=position_x, y=position_y + 0.03,
                text=str(position_num) if position_num != -1 else "?",
                fontsize=fontsize, ha='center', va='bottom')
            anomer_text = symbol_nomenclature.draw_text(
                ax, anomer_x, anomer_y + 0.03,
                r'${}$'.format(anomer_symbol_map.get(
                    child.tree.anomer, child.tree.anomer.name)),
                fontsize=fontsize, ha='center', va='bottom')
            self.data['text'][self.id, child.id]['linkage'] = [position_text, anomer_text]
        elif style == 'centered':
            symbol_nomenclature.draw_text(
                ax=ax, x=(position_x + anomer_x) / 2.0, y=(position_y + anomer_y) / 2.0,
                text=(str(position_num) if position_num != -1 else "?") + r'${}$'.format(anomer_symbol_map.get(
                    child.tree.anomer, child.tree.anomer.name)),
                fontsize=fontsize, ha='center', va='bottom')
        else:
            raise ValueError(
                "Unrecognized label_style position style %r" % (style, ))

    def draw_linkage_annotations(self, at=(0, 0), ax=None, symbol_nomenclature=None, child=None, **kwargs):
        '''
        Draw the parent outlink position and the child anomer symbol

        Parameters
        ----------
        at: :class:`tuple`
            The x, y coordinates at which to draw the glycan's reducing_end. Defaults to `(0, 0)`
        ax: :class:`matplotlib.axes.Axis`
            A matplotlib axis on which to draw. Useful for adding glycans to existing figures. If
            an axis is not provided, a new figure will be created and used.
        symbol_nomenclature: |str|
            A string mapping to the symbol nomenclature to use. Defaults to |None|. When `symbol_nomenclature`
            is |None|, `cfg` will be used.
        child: :class:`DrawTreeNode`
            The child node to draw relative to.
        '''
        if child is None:
            for child in self:
                self.draw_linkage_annotations(
                    at=at, ax=ax, symbol_nomenclature=symbol_nomenclature,
                    child=child, **kwargs)
                child.draw_linkage_annotations(
                    at=at, ax=ax, symbol_nomenclature=symbol_nomenclature, **kwargs)
            return
        self._draw_linkage_label(child, at, ax, symbol_nomenclature, **kwargs)

    # Tree Layout Helpers
    def left(self):
        return self.thread or len(self.children) and self.children[0]

    def right(self):
        return self.thread or len(self.children) and self.children[-1]

    def lbrother(self):
        n = None
        if self.parent:
            for node in self.parent.children:
                if node == self:
                    return n
                else:
                    n = node
        return n

    @property
    def lmost_sibling(self):
        if not self._lmost_sibling and self.parent and self != self.parent.children[0]:
            self._lmost_sibling = self.parent.children[0]
        return self._lmost_sibling

    def __str__(self):  # pragma: no cover
        return "%s: id=%s  x=%s mod=%s" % (self.tree, self.id, self.x, self.mod)

    def __repr__(self):  # pragma: no cover
        return self.__str__()

    def extrema(self, at=(0, 0), xmin=INF, xmax=-INF, ymin=INF, ymax=-INF, visited=None):
        '''
        Finds the most extreme points describing the area this node and its children occupy.

        Parameters
        ----------
        at: tuple of int
            The position to orient relative to in the xy-plane. A tuple of integers corresponding
            to the center x and y respectively.
        xmin: float
        xmax: float
        ymin: float
        ymax: float
            Used internally to pass along current extreme values recursively to child nodes.

        Returns
        -------
        xmin: float
        xmax: float
        ymin: float
        ymax: float

        '''
        if visited is None:
            visited = set()
        if self.id in visited:
            return xmin, xmax, ymin, ymax
        visited.add(self.id)
        x, y = self.coords(at)
        if x < xmin:
            xmin = x
        if x > xmax:
            xmax = x
        if y < ymin:
            ymin = y
        if y > ymax:
            ymax = y
        for child in sorted(self.children, key=lambda x: len(x.children)):
            n_xmin, n_xmax, n_ymin, n_ymax = child.extrema(at, xmin, xmax, ymin, ymax, visited=visited)
            xmin = min(xmin, n_xmin)
            ymin = min(ymin, n_ymin)
            xmax = max(xmax, n_xmax)
            ymax = max(ymax, n_ymax)
        return xmin, xmax, ymin, ymax

    def get(self, node_id):
        return get(self, node_id)

    def get_link_pair(self, link_id):
        return get_link_pair(self, link_id)

    def annotate_fragment(self, fragment, ax=None, color='red', label=True, **kwargs):
        if ax is None:
            ax = self.axes
        return BondCleavageArtist(self, fragment, ax, color=color, label=label, **kwargs)

    def set_transform(self, transform):
        """Apply the passed Affine2D to each graphical unit. The transform is additive
        with any previous affine transformation.

        where `M` is the value of `transform.get_matrix()`

        Parameters
        ----------
        transform : matplotlib.transforms.Affine2D
        """
        current_transform = self.transform
        if current_transform is not None:
            transform = transform + current_transform
        self.data["transform"] = transform

        # NOTE: Merge the aggregated transform of the graph
        # structure with the transformation induced by the Axes object
        # which draws each component. If this is skipped, all of the
        # geometries are wrong and the drawing will fail (or not even appear)
        transform = transform + self.axes.transData
        for i, patch_set in self.data['patches'].items():
            for entity in patch_set:
                for p in entity:
                    if isinstance(p, tuple):
                        [_.set_transform(transform) for _ in p]
                    else:
                        p.set_transform(transform)
        for i, patch_set in self.data['lines'].items():
            for p in patch_set:
                if isinstance(p, tuple):
                    [_.set_transform(transform) for _ in p]
                else:
                    p.set_transform(transform)
        for i, groups in self.data['text'].items():
            text = flatten(groups.values())
            [t.set_transform(transform) for t in text]

    def get_transform(self):
        if "transform" in self.data:
            trans = self.data["transform"]
            return trans
        else:
            return None

    transform = property(get_transform, set_transform)

    def _compute_bounding_box(self, pad_factor=0.2):
        xmin, xmax, ymin, ymax = self.extrema()
        dx = (max(xmin, xmax) - min(xmin, xmax)) * pad_factor
        dy = (max(ymin, ymax) - min(ymin, ymax)) * pad_factor
        return (xmin - dx, xmax + dx), (ymin - dy, ymax + dy)

    def update_text_position(self, degrees=0, preserve_orientation=True):
        '''
        This is a hack needed to orient text appropriately when the tree is rotated after labels
        are inscribed, as any new transform will override the correction.
        '''
        current_transform = self.transform
        if current_transform is None:
            current_transform = mtransforms.Affine2D()
        for i, groups in self.data['text'].items():
            for t in flatten(groups.values()):
                cent = (centroid(t.get_path().transformed(current_transform)))
                trans = current_transform.frozen().get_affine().translate(0, 0).rotate_deg(degrees=degrees)
                new_cent = centroid(t.get_path().transformed(trans))
                shift = cent - new_cent
                trans.translate(*shift)
                trans += self.axes.transData
                t.set_transform(trans)


class DrawTree(object):  # pragma: no cover
    def __init__(self, structure, figure=None, ax=None, layout=BalancedTreeLayout,
                 symbol_nomenclature=snfg_symbols.SNFGNomenclature(), **kwargs):
        self.structure = structure
        self.root = DrawTreeNode(root(structure))
        self.figure = figure
        self.ax = ax
        self.layout_scheme = layout
        self.symbol_nomenclature = symbol_nomenclature

    def layout(self, **kwargs):
        self.layout_scheme(self.root).layout()
        self.root.fix_special_cases()

    def draw(self, ax=None, center=True, **kwargs):
        if ax is not None:
            self.ax = ax
            self.figure = ax.get_figure()
        else:
            self.figure, self.ax = plt.subplots()

        at = (0, 0)

        self.root.draw(ax=self.ax, at=at, symbol_nomenclature=self.symbol_nomenclature)
        if center:
            xmin, xmax, ymin, ymax = self.root.extrema(at=at)
            ax = self.ax
            ax.set_xlim(sign(xmin) * (abs(xmin) + 2), (1 * abs(xmax) + 2))
            ax.set_ylim(sign(ymin) * (abs(ymin) + 2), (1 * abs(ymax) + 2))
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)

    def get(self, node_id):
        return self.root.get(node_id)

    def get_link_pair(self, link_id):
        return self.root.get_link_pair(link_id)

    @property
    def transform(self):
        return self.root.transform

    def run(self):
        self.layout()
        self.draw(ax=self.ax)
        self.figure.canvas.draw()

    def redraw(self, **kwargs):
        self.draw(**kwargs)
        self.figure.canvas.draw()


def plot(tree, at=(0, 0), ax=None, orientation='h', center=False, label=False,
         symbol_nomenclature='snfg', layout='balanced', layout_args=None, **kwargs):
    '''
    Draw the parent outlink position and the child anomer symbol

    Parameters
    ----------
    tree: Glycan or Monosaccharide
        The saccharide structure to draw.
    orientation: str
        A string corresponding to `h` or `horizontal` will draw the glycan horizontally
        from right to left. `v` or `vertical` will draw the glycan from bottom to top.
        Defaults to `h`
    at: tuple
        The x, y coordinates at which to draw the glycan's reducing_end. Defaults to `(0, 0)`
    ax: :class:`matplotlib.axes.Axis`
        A matplotlib axis on which to draw. Useful for adding glycans to existing figures. If
        an axis is not provided, a new figure will be created and used.
    center: bool
        Should the plot limits be centered around this glycan? Defaults to |False| but will be
        be set to |True| if `ax` is |None|
    label: bool
        Should the bond annotations for `tree` be drawn? Defaults to |False|
    symbol_nomenclature: str or :class:`~.SymbolicNomenclatureBase`
        Either the name of a symbol nomenclature or an instance of :class:`~.SymbolicNomenclatureBase`.
        Recognized names are "cfg", "snfg", and "iupac".
    layout: str or :class:`~.TreeLayoutBase`
        Either the name of a tree layout strategy or an instance of :class:`~.TreeLayoutBase`.
        Recognized names are "balanced" and "topological".
    layout_args: dict
        Extra parameters to pass to :class:`~.TreeLayoutBase`. It may take its own "transform" parameter
        passed here.
    scale: (float, float) or float
        Scale coefficients. Pass a pair to scale x and y dimensions respectively.
    symbol_scale: float, optional
        A scale coefficient for the size of the monosaccharide symbols, without changing the overall size
        of the layout.
    edge_weight: float, optional
        The line weight of the edges connecting monosaccharide symbols.
    label_style: str, optional
        The layout style used for drawing edge labels. May be one of "split" or "centered", where
        "split" will place the position close to the parent symbol and the anomer close to the
        child residue and "centered" places the symbols together at the centroid of the edge.
        Defaults to 'split'.
    '''

    scale = DEFAULT_SYMBOL_SCALE_FACTOR * kwargs.pop('symbol_scale', 1.0)
    transform_scale = kwargs.pop("scale", (1,))
    try:
        transform_scale = tuple(transform_scale)
    except Exception:
        transform_scale = (transform_scale, transform_scale)

    if isinstance(symbol_nomenclature, basestring):
        symbol_nomenclature = nomenclature_map.get(symbol_nomenclature)
    elif symbol_nomenclature is None:
        symbol_nomenclature = snfg_symbols.SNFGNomenclature()

    if isinstance(layout, basestring):
        layout = layout_map.get(layout)

    tree_root = root(tree)
    dtree = DrawTreeNode(tree_root)

    layout_algorithm = layout(dtree)
    layout_algorithm.layout(**(layout_args or {}))

    fig = None
    # Create a figure if no axes are provided
    if ax is None:
        fig, ax = plt.subplots()
        # at = (0, 0)
        center = True
        ax.axis('off')
    (layout_transform,
     substituent_transform,
     _) = symbol_nomenclature.get_layout_transform(orientation)
    layout_algorithm.transform(layout_transform)
    dtree.draw(at=at, ax=ax, scale=scale, label=False,
               symbol_nomenclature=symbol_nomenclature,
               annotation_transform=substituent_transform,
               **kwargs)
    dtree.axes = ax
    dtree.data['orientation'] = orientation
    dtree.set_transform(mtransforms.Affine2D())
    if label:
        dtree.draw_linkage_annotations(
            at=at, ax=ax, scale=scale,
            symbol_nomenclature=symbol_nomenclature,
            **kwargs)
    # If the figure is stand-alone, center it
    if fig is not None or center:
        xlim, ylim = _bounding_box_to_limits(dtree)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
    return (dtree, ax)
