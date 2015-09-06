import re
import operator
from collections import defaultdict
from uuid import uuid4
from itertools import chain
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
from glypy import monosaccharides
from glypy.utils import root
from glypy.io.nomenclature import identity

from .buchheim import buchheim
from .topological_layout import layout as topological
from . import cfg_symbols, iupac_symbols
from .geometry import centroid

nomenclature_map = {
    "cfg": cfg_symbols,
    'iupac': iupac_symbols
}

layout_map = {
    "balanced": buchheim,
    "topological": topological
}

DEFAULT_SYMBOL_SCALE_FACTOR = 0.17


#: :data:`special_cases` contains a list of names for
#: special case monosaccharides
special_cases = [monosaccharides["Fuc"], monosaccharides["Xyl"]]


def flatten(iterable):
    acc = []
    for item in iterable:
        try:
            acc.extend(flatten(item))
        except:
            acc.append(item)
    return acc


def is_special_case(node):
    '''
    Check to see if `node` is a special case which requires a different layout
    scheme. See `special_cases`
    '''
    for case in special_cases:
        if identity.is_a(node, case) and len(list(node.children())) == 0:
            return True
    return False


def box(node, ax=None):
    if ax is None:
        ax = node.axes
    xmin, xmax, ymin, ymax = node.extrema()
    ax.plot((xmin-0.2, xmax+0.2), (ymin-0.2, ymin-0.2), c='b')
    ax.plot((xmin-0.2, xmax+0.2), (ymax+0.2, ymax+0.2), c='b')
    ax.plot((xmin-0.2, xmin-0.2), (ymin-0.2, ymax+0.2), c='b')
    ax.plot((xmax+0.2, xmax+0.2), (ymin-0.2, ymax+0.2), c='b')


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
        x, y = node.coords()
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

        self.children = [resolve_creation_cycle(c[1], self, depth+1, i+1, c[0], visited)
                         for i, c
                         in enumerate(tree.children())]

        self.patches = []
        # A node owns all lines originating from it.
        self.lines = []

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
        self.mask_special_cases = False
        for child in self.children:
            if is_special_case(child.tree):
                if child.parent_linkage > 3:
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

    def draw(self, at=(0, 0), ax=None, symbol_nomenclature="cfg", label=True, **kwargs):
        if isinstance(symbol_nomenclature, basestring):
            symbol_nomenclature = nomenclature_map[symbol_nomenclature]
        if ax is None:
            ax = self.axes
            if ax is None:
                fig, ax = plt.subplots()
        self.draw_branches(at=at, ax=ax, symbol_nomenclature=symbol_nomenclature,
                           label=label, visited=set(), **kwargs)
        self.draw_nodes(at=at, ax=ax, symbol_nomenclature=symbol_nomenclature,
                        label=label, visited=set(), **kwargs)

    def coords(self, *args, **kwargs):
        return self.x, self.y

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
            is |None|, `cfg` will be used.
        '''
        if visited is None:
            visited = set()
        if ax is None:
            ax = self.axes
            if ax is None:
                fig, ax = plt.subplots()
        x, y = self.coords()
        if self.id in visited:
            return x, y
        visited.add(self.id)
        for child in self:
            cx, cy = child.draw_branches(at, ax=ax,
                                         symbol_nomenclature=symbol_nomenclature,
                                         label=label, visited=visited, **kwargs)
            patch = symbol_nomenclature.line_to(
                ax, x, y, cx, cy,
                color=kwargs.get('link_color', symbol_nomenclature.default_line_color),
                zorder=1)
            self.data['lines'][self.id, child.id] = [patch]
            self.lines.append(patch)
            # if label:
            #     self.draw_linkage_annotations(at=at, ax=ax,
            #                                   symbol_nomenclature=symbol_nomenclature,
            #                                   child=child, visited=visited, **kwargs)
        return x, y

    def draw_nodes(self, at=(0, 0), ax=None, symbol_nomenclature=None,
                   label=True, visited=None, **kwargs):
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
        x, y = self.coords(at)
        residue_elements, substituent_elements = symbol_nomenclature.draw(self.tree, x, y, ax, tree_node=self,
                                                                          **kwargs)
        self.patches = self.data["patches"][self.tree.id] = [residue_elements]
        self.substituents_text = self.data["text"]["substituents"][self.id] = [substituent_elements]
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
            child.draw_nodes(at, ax=ax, symbol_nomenclature=symbol_nomenclature,
                             visited=visited, **kwargs)

    def draw_linkage_annotations(self, at=(0, 0), ax=None,
                                 symbol_nomenclature=None, child=None, **kwargs):
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
        sx, sy = self.coords(at)
        if child is None:
            for child in self:
                self.draw_linkage_annotations(
                    at=at, ax=ax, symbol_nomenclature=symbol_nomenclature,
                    child=child, **kwargs)
                child.draw_linkage_annotations(
                    at=at, ax=ax, symbol_nomenclature=symbol_nomenclature, **kwargs)
            return
        cx, cy = child.coords(at)

        fontsize = kwargs.get("fontsize", 6)

        position_x = (sx * 0.5 + cx * 0.5)
        position_y = (sy * 0.7 + cy * 0.3)
        # if sx == cx:
        #     position_y += -0.01
        if sx > cx:
            position_x = (sx * 0.7) + (cx * 0.3) - 0.2
            position_y = (sy * 0.8 + cy * 0.2)
        elif sx < cx:
            position_x = (sx * 0.7) + (cx * 0.3) + 0.1
            position_y = (sy * 0.8 + cy * 0.2)
        else:
            position_x += 0.05

        position_num = -1
        for pos, link in self.tree.links.items():
            if link.is_child(child.tree):
                position_num = pos
                break

        anomer_x = ((sx * 0.4 + cx * 0.6))# + (-sign(sx) if cx <= sx else sign(sx)) * 0.11
        anomer_y = ((sy * 0.4 + cy * 0.6))
        if sx > cx:
            anomer_x = (sx * 0.3) + (cx * 0.7) - 0.2
        elif sx < cx:
            anomer_x = (sx * 0.3) + (cx * 0.7) + 0.05
        else:
            anomer_x += 0.05
        if sy == cy:
            _x = (position_x * 0.9) + (anomer_x * 0.1)
            position_x = anomer_x = _x
            anomer_y = ((sy * 0.5 + cy * 0.5)) + 0.10
            position_y = ((sy * 0.5 + cy * 0.5)) - 0.16

        overlapping = flatten(
            (chain.from_iterable(self.data['patches'][self.id]),
             chain.from_iterable(self.data['patches'][child.id]),
             self.data['lines'][self.id, child.id]))

        position_text = symbol_nomenclature.draw_text(
            ax=ax, x=position_x, y=position_y, text=str(position_num),
            fontsize=fontsize, minimize_overlap=overlapping)
        anomer_text = symbol_nomenclature.draw_text(
            ax, anomer_x, anomer_y, r'$\{}$'.format(child.tree.anomer.name), fontsize=fontsize + 2,
            minimize_overlap=overlapping)
        self.data['text'][self.id, child.id]['linkage'] = [position_text, anomer_text]

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

    def get_lmost_sibling(self):
        if not self._lmost_sibling and self.parent and self != self.parent.children[0]:
            self._lmost_sibling = self.parent.children[0]
        return self._lmost_sibling
    lmost_sibling = property(get_lmost_sibling)

    def __str__(self):
        return "%s: id=%s  x=%s mod=%s" % (self.tree, self.id, self.x, self.mod)

    def __repr__(self):
        return self.__str__()

    def extrema(self, at=(0, 0), xmin=float('inf'), xmax=-float("inf"),
                ymin=float('inf'), ymax=-float("inf"), visited=None):
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

    def draw_cleavage(self, fragment=None, at=(0, 0), ax=None, scale=0.1, color='red', label=True):
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
            px, py = parent.coords(at)
            cx, cy = child.coords(at)

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
            cx, cy = target.coords(at)
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

    def transform(self, transform):
        base_transform = transform
        current_transform = self.data.get("transform")
        if current_transform is not None:
            transform = transform + current_transform
        self.data["transform"] = transform
        transform = transform + self.axes.transData
        for i, patches in self.data['patches'].items():
            for entity in patches:
                for p in entity:
                    if isinstance(p, tuple):
                        [_.set_transform(transform) for _ in p]
                    else:
                        p.set_transform(transform)
        for i, patches in self.data['lines'].items():
            for p in patches:
                if isinstance(p, tuple):
                    [_.set_transform(transform) for _ in p]
                else:
                    p.set_transform(transform)
        for i, groups in self.data['text'].items():
            text = flatten( groups.values())
            [t.set_transform(transform) for t in text]

        mat = base_transform.get_matrix()
        for node in breadth_first_traversal(self):
            x, y, w = mat.dot(np.array((node.x, node.y, 1)))
            x /= w
            y /= w
            node.x = x
            node.y = y

        if self.data['orientation'] == 'h':
            self._rotate_text(-90)

    def _rotate_text(self, degrees):
        '''
        This is a hack needed to orient text appropriately when the tree is rotated after labels
        are inscribed, as any new transform will override the correction.

        Two solutions exist:
            Generate all text placements, then apply all transforms, and then inscribe the actual
            text paths once they are 

        '''
        current_transform = self.data.get("transform")
        for i, groups in self.data['text'].items():
            for t in flatten(groups.values()):
                cent = (centroid(t.get_path().transformed(current_transform)))
                trans = current_transform.frozen().get_affine().translate(0, 0).rotate_deg(degrees=degrees)
                new_cent = centroid(t.get_path().transformed(trans))
                shift = cent - new_cent
                trans.translate(*shift)
                trans += self.axes.transData
                t.set_transform(trans)


class DrawTree(object):
    def __init__(self, structure, figure=None, ax=None, layout=buchheim, symbol_nomenclature=cfg_symbols, **kwargs):
        self.structure = structure
        self.root = DrawTreeNode(root(structure))
        self.figure = figure
        self.ax = ax
        self.layout_scheme = layout
        self.symbol_nomenclature = symbol_nomenclature

    def layout(self, **kwargs):
        self.layout_scheme.layout(self.root)
        self.root.fix_special_cases()

    def draw(self, ax=None, at=(1, 1), center=True, **kwargs):
        if ax is not None:
            self.ax = ax
            self.figure = ax.get_figure()
        else:
            self.figure, self.ax = plt.subplots()

        self.root.draw(ax=self.ax, at=at, symbol_nomenclature=self.symbol_nomenclature)
        if center:
            xmin, xmax, ymin, ymax = self.root.extrema(at=at)
            ax = self.ax
            ax.set_xlim(sign(xmin) * (abs(xmin) + 2), (1 * abs(xmax) + 2))
            ax.set_ylim(sign(ymin) * (abs(ymin) + 2), (1 * abs(ymax) + 2))
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)

    def transform(self, *args, **kwargs):
        self.root.transform(*args, **kwargs)

    def run(self, at=(1, 1)):
        self.layout()
        self.draw(ax=self.ax, at=at)
        self.figure.canvas.draw()

    def redraw(self, **kwargs):
        self.draw(**kwargs)
        self.figure.canvas.draw()


def plot(tree, at=(1, 1), ax=None, orientation='h', center=False, label=False,
         symbol_nomenclature='cfg', layout='balanced', **kwargs):
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

    if isinstance(symbol_nomenclature, basestring):
        symbol_nomenclature = nomenclature_map.get(symbol_nomenclature)

    if isinstance(layout, basestring):
        layout = layout_map.get(layout)

    tree_root = root(tree)
    dtree = DrawTreeNode(tree_root)
    if layout == topological:
        for node in breadth_first_traversal(dtree):
            node.mask_special_cases = False

    layout(dtree)
    if layout != topological:
        dtree.fix_special_cases()

    fig = None
    # Create a figure if no axes are provided
    if ax is None:
        fig, ax = plt.subplots()
        at = (0, 0)
        center = True
        ax.axis('off')
    dtree.draw(at=at, ax=ax, scale=scale, label=False,
               symbol_nomenclature=symbol_nomenclature, **kwargs)
    dtree.axes = ax
    dtree.data['orientation'] = orientation
    if label:
        dtree.draw_linkage_annotations(
            at=at, ax=ax, scale=scale,
            symbol_nomenclature=symbol_nomenclature, **kwargs)
    if orientation in {"h", "horizontal"}:
        dtree.transform(mtransforms.Affine2D().rotate_deg(90))
        dtree._rotate_text(-90)
    # If the figure is stand-alone, center it
    if fig is not None or center:
        xmin, xmax, ymin, ymax = dtree.extrema(at=at)
        ax.set_xlim(xmin - 2, xmax + 2)
        ax.set_ylim(ymin - 2, ymax + 2)
        ax.autoscale_view()
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
    return (dtree, ax)
