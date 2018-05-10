# Adapted from https://github.com/llimllib/pymag-trees/blob/master/buchheim.py

from .geometry import TreeLayoutBase


class BalancedTreeLayout(TreeLayoutBase):

    def layout_tree(self, **kwargs):
        buchheim(self.root)


def buchheim(tree, visited=None):
    if visited is None:
        visited = set()
    dt = first_walk(tree, visited=visited)
    visited = set()
    min = second_walk(dt, visited=visited)
    if min < 0:
        third_walk(dt, -min)
    return dt


layout = buchheim


def third_walk(tree, n, visited=None):  # pragma: no cover
    if visited is None:
        visited = set()

    tree.x += n

    if tree.id in visited:
        return

    visited.add(tree.id)

    for c in tree.children:
        third_walk(c, n, visited=visited)


def first_walk(v, distance=1., visited=None):
    if visited is None:  # pragma: no cover
        visited = set()
    if v.id in visited:  # pragma: no cover
        return v
    visited.add(v.id)
    if len(v.children) == 0:
        if v.lmost_sibling:
            v.x = v.lbrother().x + distance
        else:
            v.x = 0.
    else:
        default_ancestor = v.children[0]
        for w in v.children:
            first_walk(w, distance, visited)
            default_ancestor = apportion(w, default_ancestor, distance)
        execute_shifts(v)

        midpoint = (v.children[0].x + v.children[-1].x) / 2

        # ell = v.children[0]
        # arr = v.children[-1]
        w = v.lbrother()
        if w:
            v.x = w.x + distance
            v.mod = v.x - midpoint
        else:
            v.x = midpoint
    return v


def apportion(v, default_ancestor, distance):
    w = v.lbrother()
    if w is not None:
        # in buchheim notation:
        # i == inner; o == outer; r == right; l == left; r = +; l = -
        vir = vor = v
        vil = w
        vol = v.lmost_sibling
        sir = sor = v.mod
        sil = vil.mod
        sol = vol.mod
        while vil.right() and vir.left():
            vil = vil.right()
            vir = vir.left()
            vol = vol.left()
            vor = vor.right()
            vor.ancestor = v
            shift = (vil.x + sil) - (vir.x + sir) + distance
            if shift > 0:  # pragma: no cover
                move_subtree(ancestor(vil, v, default_ancestor), v, shift)
                sir = sir + shift
                sor = sor + shift
            sil += vil.mod
            sir += vir.mod
            sol += vol.mod
            sor += vor.mod
        if vil.right() and not vor.right():
            vor.thread = vil.right()
            vor.mod += sil - sor
        else:
            if vir.left() and not vol.left():
                vol.thread = vir.left()
                vol.mod += sir - sol
            default_ancestor = v
    return default_ancestor


def move_subtree(wl, wr, shift):  # pragma: no cover
    subtrees = wr.number - wl.number
    wr.change -= shift / subtrees
    wr.shift += shift
    wl.change += shift / subtrees
    wr.x += shift
    wr.mod += shift


def execute_shifts(v):
    shift = change = 0
    for w in v.children[::-1]:
        w.x += shift
        w.mod += shift
        change += w.change
        shift += w.shift + change


def ancestor(vil, v, default_ancestor):  # pragma: no cover
    # the relevant text is at the bottom of page 7 of
    # "Improving Walker's Algorithm to Run in Linear Time" by Buchheim et al, (2002)
    # http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.16.8757&rep=rep1&type=pdf
    if vil.ancestor in v.parent.children:
        return vil.ancestor
    else:
        return default_ancestor


def second_walk(v, m=0, depth=0, min=None, visited=None):
    if visited is None:
        visited = set()
    v.x += m
    v.y = depth

    if v.id in visited:
        return v.x

    visited.add(v.id)

    if min is None or v.x < min:
        min = v.x

    for w in v.children:
        min = second_walk(w, m + v.mod, depth + 1, min, visited=visited)

    return min
