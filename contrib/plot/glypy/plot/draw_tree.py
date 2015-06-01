import re
import operator
from collections import defaultdict
from uuid import uuid4

import matplotlib.pyplot as plt
from glypy import monosaccharides
from glypy.utils import root
from glypy.structure import Glycan
from glypy.io.nomenclature import identity
from glypy.algorithms.database import GlycanRecordBase

from .buchheim import buchheim
from . import cfg_symbols


nomenclature_map = {
    "cfg": cfg_symbols
}


DEFAULT_TREE_SCALE_FACTOR = 0.65
DEFAULT_SYMBOL_SCALE_FACTOR = 0.1


#: :data:`special_cases` contains a list of names for
#: special case monosaccharides
special_cases = [monosaccharides["Fuc"], monosaccharides["Xyl"]]


def is_special_case(node):
    '''
    Check to see if `node` is a special case which requires a different layout
    scheme. See `special_cases`
    '''
    for case in special_cases:
        if identity.is_a(node, case) and len(list(node.children())) == 0:
            return True
    return False


def sign(x, zero_dir=-1):
    '''Return 1 if x > 0, -1 if x < 0, and `zero_dir` otherwise'''
    if x > 0:
        return 1
    if x < 0:
        return -1
    return zero_dir


def make_gid(s):
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
        x, y = node.coords('h')
        ax.text(x, y + 0.2, labeler(node.tree), color='red')
    # Update the plot
    ax.get_figure().canvas.draw()
    return tree, ax


def find_link(parent, child):
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


def resolve_creation_cycle(node, parent, depth, i, parent_linkage, visited):
    if node.id in visited:
        return visited[node.id]
    else:
        return DrawTreeNode(node, parent, depth+1, i+1, parent_linkage, visited)


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
        self.tree = tree
        self.parent = parent
        self.parent_linkage = parent_linkage

        self._axes = None

        self.mask_special_cases = True
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
            self.data = defaultdict(dict)
            self.uuid = uuid4().hex

        self.children = [resolve_creation_cycle(c[1], self, depth+1, i+1, c[0], visited)
                         for i, c
                         in enumerate(tree.children())]

    def __iter__(self):
        return iter(self.children)

    @property
    def children(self):
        if self.mask_special_cases:
            return [child for child in self._children if not (is_special_case(child.tree) and len(child.children) == 0)]
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
        n = 0
        self.mask_special_cases = False
        for child in self.children:
            if is_special_case(child.tree):
                if n == 0:
                    child.y = self.y
                    child.x = self.x + offset
                elif n == 1:
                    child.y = self.y
                    child.x = self.x - offset
                else:
                    raise Exception("Don't know how to handle more than two special case child nodes.")
            child.fix_special_cases(offset, visited)

    def draw(self, orientation="h", at=(0, 0), ax=None, symbol_nomenclature="cfg", label=True, **kwargs):
        if isinstance(symbol_nomenclature, basestring):
            symbol_nomenclature = nomenclature_map[symbol_nomenclature]
        if ax is None:
            ax = self.axes
            if ax is None:
                fig, ax = plt.subplots()
        self.draw_branches(orientation, at=at, ax=ax, symbol_nomenclature=symbol_nomenclature,
                           label=label, visited=set(), **kwargs)
        self.draw_nodes(orientation, at=at, ax=ax, symbol_nomenclature=symbol_nomenclature,
                        label=label, visited=set(), **kwargs)

    def coords(self, orientation="h", at=(0, 0)):
        if orientation in {"h", "horizontal"}:
            x = -self.y + at[1]
            y = self.x + at[0]
        elif orientation in {'v', "vertical"}:
            x = self.x + at[0]
            y = self.y + at[1]
        else:
            raise Exception("Unknown Orientation: {}".format(orientation))
        return x, y

    def map_coords(self, orientation, x, y):
        if orientation in {"h", "horizontal"}:
            x = -y
            y = x
        elif orientation in {'v', "vertical"}:
            x = x
            y = y
        return x, y

    def draw_branches(self, orientation="h", at=(0, 0), ax=None, symbol_nomenclature=None,
                      label=True, visited=None, **kwargs):
        '''
        Draw the edges linking parent nodes to child nodes. Also draws the parent outlink position
        and the child anomer symbol
        Parameters
        ----------
        orientation: |str|
            A string corresponding to `h` or `horizontal` will draw the glycan horizontally
            from right to left. `v` or `vertical` will draw the glycan from bottom to top.
            Defaults to `h`
        at: :class:`tuple`
            The x, y coordinates at which to draw the glycan's reducing_end. Defaults to `(0, 0)`
        ax: :class:`matplotlib.axes.Axis`
            A matplotlib axis on which to draw. Useful for adding glycans to existing figures. If
            an axis is not provided, a new figure will be created and used.
        symbol_nomenclature: |str|
            A string mapping to the symbol nomenclature to use. Defaults to |None|. When `symbol_nomenclature`
            is |None|, `cfg` will be used.
        '''
        if visited is None:
            visited = set()
        if ax is None:
            ax = self.axes
            if ax is None:
                fig, ax = plt.subplots()
        x, y = self.coords(orientation, at)
        if self.id in visited:
            return x, y
        visited.add(self.id)
        for child in self:
            cx, cy = child.draw_branches(orientation, at, ax=ax,
                                         symbol_nomenclature=symbol_nomenclature,
                                         label=label, visited=visited, **kwargs)
            symbol_nomenclature.line_to(ax, x, y, cx, cy, color='black', zorder=1)
            if label:
                self.draw_linkage_annotations(orientation=orientation, at=at, ax=ax,
                                              symbol_nomenclature=symbol_nomenclature,
                                              child=child, visited=visited, **kwargs)
        return x, y

    def draw_nodes(self, orientation='h', at=(0, 0), ax=None, symbol_nomenclature=None,
                   label=True, visited=None, **kwargs):
        '''
        Draw the symbol representing the individual monosaccharide and its substituents.

        Parameters
        ----------
        orientation: |str|
            A string corresponding to `h` or `horizontal` will draw the glycan horizontally
            from right to left. `v` or `vertical` will draw the glycan from bottom to top.
            Defaults to `h`
        at: :class:`tuple`
            The x, y coordinates at which to draw the glycan's reducing_end. Defaults to `(0, 0)`
        ax: :class:`matplotlib.axes.Axis`
            A matplotlib axis on which to draw. Useful for adding glycans to existing figures. If
            an axis is not provided, a new figure will be created and used.
        symbol_nomenclature: |str|
            A string mapping to the symbol nomenclature to use. Defaults to |None|. When `symbol_nomenclature`
            is |None|, `cfg` will be used.
        '''
        if visited is None:
            visited = set()
        if ax is None:
            ax = self.axes
            if ax is None:
                fig, ax = plt.subplots()
        if self.id in visited:
            return
        visited.add(self.id)
        x, y = self.coords(orientation, at)
        residue_elements, substituent_elements = symbol_nomenclature.draw(self.tree, x, y, ax, tree_node=self,
                                                                          orientation=orientation, **kwargs)
        self.data["patches"][self.tree.id] = [residue_elements, substituent_elements]
        self.data["position"][self.tree.id] = x, y
        for i, res_el in enumerate(residue_elements):
            if isinstance(res_el, tuple):
                [el.set_gid(self.uuid + '-' + (str(self.id) + '-node')) for el in res_el]
            else:
                res_el.set_gid(self.uuid + '-' + str(self.id) + '-node')
        for i, sub_el in enumerate(substituent_elements):
            if isinstance(sub_el, tuple):
                [el.set_gid(self.uuid + '-' + (str(self.id) + "-subst-" + str(i))) for el in sub_el]
            else:
                sub_el.set_gid(self.uuid + '-' + (str(self.id) + "-subst-" + str(i)))
        for child in self:
            child.draw_nodes(orientation, at, ax=ax, symbol_nomenclature=symbol_nomenclature,
                             visited=visited, **kwargs)

    def draw_linkage_annotations(self, orientation='h', at=(0, 0), ax=None,
                                 symbol_nomenclature=None, child=None, **kwargs):
        '''
        Draw the parent outlink position and the child anomer symbol

        Parameters
        ----------
        orientation: |str|
            A string corresponding to `h` or `horizontal` will draw the glycan horizontally
            from right to left. `v` or `vertical` will draw the glycan from bottom to top.
            Defaults to `h`
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
        sx, sy = self.coords(orientation, at)
        cx, cy = child.coords(orientation, at)

        position_x = ((sx * 0.6 + cx * 0.4)) + (-sign(sx) if cx <= sx else sign(sx)) * 0.07
        position_y = ((sy * 0.6 + cy * 0.4))
        if sx == cx:
            position_y += 0.006
        elif (sy > cy):
            position_y += -0.11
        elif (sy < cy):
            position_y += -0.15
        else:
            position_y += -0.12
        position_num = -1
        for pos, link in self.tree.links.items():
            if link.is_child(child.tree):
                position_num = pos
                break
        symbol_nomenclature.draw_text(ax=ax, x=position_x, y=position_y, text=str(position_num))
        anomer_x = ((sx * 0.2 + cx * 0.8)) + (-sign(sx) if cx <= sx else sign(sx)) * 0.07
        anomer_y = ((sy * 0.2 + cy * 0.8))
        if sx == cx:
            anomer_y += 0.006
        elif (sy > cy):
            anomer_y += -0.11
        elif (sy < cy):
            anomer_y += -0.15
        else:
            anomer_y += -0.12

        symbol_nomenclature.draw_text(ax, anomer_x, anomer_y, r'$\{}$'.format(child.tree.anomer.name), scale=0.13)

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

    def get_lmost_sibling(self):
        if not self._lmost_sibling and self.parent and self != self.parent.children[0]:
            self._lmost_sibling = self.parent.children[0]
        return self._lmost_sibling
    lmost_sibling = property(get_lmost_sibling)

    def __str__(self):
        return "%s: id=%s  x=%s mod=%s" % (self.tree, self.id, self.x, self.mod)

    def __repr__(self):
        return self.__str__()

    def extrema(self, orientation='h', at=(0, 0), xmin=float('inf'), xmax=-float("inf"),
                ymin=float('inf'), ymax=-float("inf"), visited=None):
        '''
        Finds the most extreme points describing the area this node and its children occupy.

        Parameters
        ----------
        orientation: str
            The orientation mapping symbol. May be one of {'h', 'horizontal', 'v', 'vertical'}.
            Defaults to 'h'
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
        x, y = self.coords(orientation, at)
        if x < xmin:
            xmin = x
        if x > xmax:
            xmax = x
        if y < ymin:
            ymin = y
        if y > ymax:
            ymax = y
        for child in self:
            xmin, xmax, ymin, ymax = child.extrema(orientation, at, xmin, xmax, ymin, ymax, visited=visited)
        return xmin, xmax, ymin, ymax

    def scale(self, factor=DEFAULT_TREE_SCALE_FACTOR, visited=None):
        '''Scale all existing edge lengths by `factor`, defaults to `0.65`'''
        if visited is None:
            visited = set()
        if self.id in visited:
            return
        try:
            if len(factor) == 1:
                xfactor = yfactor = factor
            else:
                xfactor, yfactor = factor
        except TypeError:
            xfactor = yfactor = factor
        visited.add(self.id)
        self.x *= xfactor
        self.y *= yfactor
        for child in self:
            child.scale(factor, visited)

    def get(self, node_id):
        return get(self, node_id)

    def get_link_pair(self, link_id):
        return get_link_pair(self, link_id)

    def draw_cleavage(self, fragment=None, orientation="h", at=(0, 0), ax=None, scale=0.1, color='red', label=True):
        '''
        .. warning::
            Here be magical numbers

        '''
        if ax is None:
            ax = self.axes
        if ax is None:
            raise Exception("`ax` is required")
        if fragment is None:
            raise Exception("`fragment` is required")
        scale *= 2

        break_targets = fragment.link_ids
        crossring_targets = fragment.crossring_cleavages

        for link_break in break_targets:
            parent, child = self.get_link_pair(link_break)
            px, py = parent.coords(orientation, at)
            cx, cy = child.coords(orientation, at)

            branch_point = False
            if py == cy and px != cx:
                center = (max(px, cx) - min(px, cx)) / 2. + min(px, cx)
                lx1, ly1 = center, py + scale
                lx2, ly2 = center, py - scale
            elif py != cy and px == cx:
                center = (max(py, cy) - min(py, cy)) / 2. + min(py, cy)
                lx1, ly1 = px + scale, center
                lx2, ly2 = px - scale, center
            else:
                branch_point = True
                xcenter = (max(px, cx) - min(px, cx)) / 2. + min(px, cx)
                ycenter = (max(py, cy) - min(py, cy)) / 2. + min(py, cy)
                if py < cy:
                    lx2, ly2 = xcenter + scale, ycenter + scale
                    lx1, ly1 = xcenter - scale, ycenter - scale
                else:
                    lx1, ly1 = xcenter - scale, ycenter + scale
                    lx2, ly2 = xcenter + scale, ycenter - scale

            line = ax.plot((lx1, lx2), (ly1, ly2), color=color, zorder=3)
            self.data['patches'][fragment.name] = line[0]
            self.data['position'][fragment.name] = (lx1, lx2), (ly1, ly2)
            line[0].set_gid(self.uuid + '-' + fragment.name)
            if fragment.kind[-1] in {"B", "C"}:
                label_x = (lx2, lx2 - scale)
                if branch_point:
                    label_x = [x - 2 * scale for x in label_x]

                line = ax.plot(label_x, (ly1, ly1), color=color, zorder=3)
                line[0].set_gid(self.uuid + '-' + fragment.name + "-direction")
                self.data['patches'][fragment.name + "_direction"] = line[0]
                self.data['position'][fragment.name + "_direction"] = label_x, (ly1, ly1)
                if label:
                    text = ax.text(label_x[0] - .4, ly1 + 0.05, fragment.fname)
                    self.data['patches'][fragment.name + "_text"] = text
                    self.data['position'][fragment.name + "_text"] = label_x[0] - .4, ly1 + 0.05
            else:
                line = ax.plot((lx2, lx2 + scale), (ly2, ly2), color=color, zorder=3)
                line[0].set_gid(self.uuid + '-' + fragment.name + "-direction")
                self.data['patches'][fragment.name + "_direction"] = line[0]
                self.data['position'][fragment.name + "_direction"] = (lx2, lx2 + scale), (ly2, ly2)

                if label:
                    self.data['patches'][fragment.name + "_text"] = ax.text(lx2 + 0.05, ly2 - 0.15, fragment.fname)
                    self.data['position'][fragment.name + "_text"] = lx2 + 0.05, ly2 - 0.15

        for crossring in crossring_targets:
            target = self.get(crossring)
            cx, cy = target.coords(orientation, at)
            line = ax.plot((cx - scale, cx + scale), (cy + scale, cy - scale), color=color, zorder=3)
            self.data['patches'][fragment.name] = line[0]
            self.data['position'][fragment.name] = (cx - scale, cx + scale), (cy + scale, cy - scale)
            line[0].set_gid(self.uuid + '-' + fragment.name.replace(",", '_'))
            annotation_name = re.sub(r'\d,\d', '', fragment.fname)
            if fragment.kind[-1] == "X":
                line = ax.plot((cx + scale, cx + 2 * scale), (cy - scale, cy - scale), color=color, zorder=3)
                line[0].set_gid(self.uuid + '-' + fragment.name + "-direction")
                self.data['patches'][fragment.name + "_direction"] = line[0]
                self.data['position'][fragment.name + "_direction"] =\
                    (cx + scale, cx + 2 * scale), (cy - scale, cy - scale)

                if label:
                    ax.text((cx + scale) - 0.4, (cy - scale) - .15, annotation_name)
            else:
                line = ax.plot((cx - scale, cx - scale * 2), (cy + scale, cy + scale), color=color, zorder=3)
                line[0].set_gid(self.uuid + '-' + fragment.name.replace(",", '_') + "-direction")
                self.data['patches'][fragment.name + "_direction"] = line[0]
                self.data['position'][fragment.name + "_direction"] =\
                    (cx - scale, cx - scale * 2), (cy + scale, cy + scale)
                if label:
                    ax.text((cx - scale) - 0.32, (cy + scale) + .035, annotation_name)


class DrawTree(object):
    def __init__(self, structure, figure=None, ax=None, **kwargs):
        self.structure = structure
        self.root = DrawTreeNode(root(structure))
        self.figure = figure
        self.axes = ax
        self.layout = kwargs.get("layout", buchheim)
        self.symbol_nomenclature = kwargs.get("symbol_nomenclature", cfg_symbols)


def plot(tree, orientation='h', at=(1, 1), ax=None, center=False, label=False, **kwargs):
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
    scale: (float, float) or float
        Node scale coefficients. Pass a pair to scale x and y dimensions respectively.
    stretch: (float, float) or float
        Edge length coefficients. Pass a pair to scale x and y dimensions respectively. Order is
        affected by `orientation`
    '''

    scale = kwargs.get("scale", DEFAULT_SYMBOL_SCALE_FACTOR)

    # if ax is not None and "scale" not in kwargs:
    #     x_max = ax.get_xlim()[1]
    #     y_max = ax.get_ylim()[1]

    tree_root = root(tree)
    dtree = DrawTreeNode(tree_root)
    buchheim(dtree)
    dtree.fix_special_cases()
    dtree.scale(kwargs.get("squeeze", kwargs.get("stretch", DEFAULT_TREE_SCALE_FACTOR)))
    fig = None
    # Create a figure if no axes are provided
    if ax is None:
        fig, ax = plt.subplots()
        at = (0, 0)
        center = True
    dtree.draw(orientation, at=at, ax=ax, scale=scale, label=label)
    # If the figure is stand-alone, center it
    if fig is not None or center:
        xmin, xmax, ymin, ymax = dtree.extrema(orientation=orientation, at=at)
        if orientation in {'h', 'horizontal'}:
            ax.set_xlim(-1 * (abs(xmin) + 2) if xmin < 1 else xmin - 2, (1 * abs(xmax) + 2))
            ax.set_ylim(-1 * (abs(ymin) + 2) if ymin < 1 else ymin - 2, (1 * abs(ymax) + 2))
        elif orientation in {'v', 'vertical'}:
            ax.set_xlim(xmin - sign(xmin, 1) * 2, xmax + 2)
            ax.set_ylim(ymin - sign(ymin, 1) * 2, ymax + 2)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
    dtree.axes = ax
    return (dtree, ax)
