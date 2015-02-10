import matplotlib.pyplot as plt

from pygly2.structure import Glycan

from .buchheim import buchheim
import cfg_symbols


nomenclature_map = {
    "cfg": cfg_symbols
}


class DrawTree(object):
    def __init__(self, tree, parent=None, depth=0, number=1):
        self.x = -1.
        self.y = depth
        self.tree = tree
        self.children = [DrawTree(c, self, depth+1, i+1)
                         for i, c
                         in enumerate(ch for p, ch in tree.children())]
        self.parent = parent
        self.thread = None
        self.mod = 0
        self.ancestor = self
        self.change = self.shift = 0
        self._lmost_sibling = None
        # Number of the node in its group of siblings 1..n
        self.number = number

    def __iter__(self):
        return iter(self.children)

    def check(self):
        print(self.x, self.y)
        for child in self:
            child.check()

    def draw(self, orientation="h", at=(0, 0), ax=None, symbol_nomenclature="cfg", **kwargs):
        if isinstance(symbol_nomenclature, basestring):
            symbol_nomenclature = nomenclature_map[symbol_nomenclature]
        if ax is None:
            fig, ax = plt.subplots()
        self.draw_branches(orientation, at=at, ax=ax, symbol_nomenclature=symbol_nomenclature, **kwargs)
        self.draw_nodes(orientation, at=at, ax=ax, symbol_nomenclature=symbol_nomenclature, **kwargs)

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

    def draw_branches(self, orientation="h", at=(0, 0), ax=None, symbol_nomenclature=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots()
        x, y = self.coords(orientation, at)
        for child in self:
            cx, cy = child.draw_branches(orientation, at, ax=ax, **kwargs)
            ax.plot((x, cx), (y, cy), color='black', zorder=1)
        return x, y

    def draw_nodes(self, orientation='h', at=(0, 0), ax=None, symbol_nomenclature=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots()
        x, y = self.coords(orientation, at)
        # ax.scatter(x, y, s=150, zorder=2)
        symbol_nomenclature.draw(self.tree, x, y, ax, **kwargs)
        for child in self:
            child.draw_nodes(orientation, at, ax=ax, symbol_nomenclature=symbol_nomenclature, **kwargs)

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
        return "%s: x=%s mod=%s" % (self.tree, self.x, self.mod)

    def __repr__(self):
        return self.__str__()

    def extrema(self, orientation='h', at=(0, 0), xmin=float('inf'), xmax=-float("inf"),
                ymin=float('inf'), ymax=-float("inf")):
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
            xmin, xmax, ymin, ymax = child.extrema(orientation, at, xmin, xmax, ymin, ymax)
        return xmin, xmax, ymin, ymax

    def scale(self, factor=0.5):
        self.x *= factor
        self.y *= factor
        for child in self:
            child.scale(factor)


def sign(x, zero_dir=-1):
    if x > 0:
        return 1
    if x < 0:
        return -1
    return zero_dir


def plot(tree, orientation='h', at=(1, 1), ax=None):
    root = tree
    if isinstance(tree, Glycan):
        root = tree.root
    dtree = DrawTree(root)
    buchheim(dtree)
    dtree.scale()
    fig = None
    if ax is None:
        fig, ax = plt.subplots()
        at = (0, 0)
    dtree.draw(orientation, at=at, ax=ax)
    if fig is not None:
        xmin, xmax, ymin, ymax = dtree.extrema(orientation=orientation, at=at)
        print(xmin, xmax, ymin, ymax)
        if orientation in {'h', 'horizontal'}:
            ax.set_xlim(-1 * (abs(xmin) + 2) if xmin < 1 else xmin - 2, (1 * abs(xmax) + 2))
            ax.set_ylim(-1 * (abs(ymin) + 2) if ymin < 1 else ymin - 2, (1 * abs(ymax) + 2))
        elif orientation in {'v', 'vertical'}:
            ax.set_xlim(xmin - sign(xmin, 1) * 2, xmax + 2)
            ax.set_ylim(ymin - sign(ymin, 1) * 2, ymax + 2)
        fig.tight_layout()
    return (dtree, ax)
