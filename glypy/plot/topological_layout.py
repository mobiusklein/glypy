from .geometry import breadth_first_traversal, centroid, make_path
import numpy as np

#: clock-face coordinates assume left-horizontal growth
index_position_shift = {
    0: (-0.0, -0.0),
    2: (-0.0, 1),
    3: (-0.5, 1),
    4: (-0.8, 0),
    5: (-0.5, -1),
    6: (-0.6, -1),
    8: (-0, -1)
}
index_position_shift = {k: np.array(v[::-1], dtype=np.float64) for k, v in index_position_shift.items()}


def _test_draw_unit_perimeter(ax):  # pragma: no cover
    ax.text(0, 0, 'origin')
    for i in index_position_shift:
        ax.text(*square_perimeter(0, 0, i), s=i)
    return ax


def square_perimeter(cx, cy, index_position, scale=1.0):
    sx, sy = index_position_shift[index_position]
    sx *= scale
    sy *= scale
    return cx - sx, cy - sy


def layout_node(node, cx, cy, parent_position=3, visited=None):
    if visited is None:
        visited = set()
    node.x, node.y = square_perimeter(cx, cy, parent_position)
    node.mask_special_cases = False
    for child_node, pos_child in zip(node.children, node.tree.children()):
        pos = pos_child[0]
        layout_node(child_node, node.x, node.y, pos)
    return node


def layout(tree, visited=None):
    if visited is None:
        visited = set()
    layout_node(tree, 0, 0, 0, visited)
    spread(tree)
    return tree


# def make_path(node):
#     pathing = []
#     last = [node.x, node.y]
#     for node in breadth_first_traversal(node):
#         pathing.append([node.x, node.y])
#         pathing.append(last)
#         last = [node.x, node.y]
#     codes = [mpath.Path.MOVETO] + [mpath.Path.LINETO for i in pathing[1:]]
#     return mpath.Path(pathing, codes)


# def centroid(path):
#     point = np.zeros_like(path.vertices[0])
#     c = 0.
#     for p in path.vertices:
#         point += p
#         c += 1
#     return point / c


def spread(tree):
    next_layer = tree.children
    if len(next_layer) == 0:
        return
    for node in next_layer:
        spread(node)

    if len(next_layer) < 2:
        return

    last = next_layer[0]
    shifted = [last]
    for node in next_layer[1:]:
        lpath = make_path(last)
        npath = make_path(node)
        if lpath.intersects_path(npath):
            lcentroid = centroid(lpath)
            ncentroid = centroid(npath)
            lxmin_lymin, lxmax_lymax = lpath.get_extents().get_points()
            nxmin_nymin, nxmax_nymax = npath.get_extents().get_points()
            if last.x == tree.x:
                delta = (nxmax_nymax[0] - ncentroid[0])
                shift_subtree(node, delta)
            elif ncentroid[0] > lcentroid[0]:
                delta = (nxmin_nymin[0] - ncentroid[0]) * 1.1
                for prev in shifted:
                    if prev.x == tree.x:
                        continue
                    pcentroid = centroid(make_path(prev))
                    if pcentroid[0] > tree.x:
                        shift_subtree(prev, delta)
            elif ncentroid[0] < lcentroid[0]:
                delta = (lxmax_lymax - lcentroid[0]) * -1.1
                for prev in shifted:
                    if prev.x == tree.x:
                        continue
                    pcentroid = centroid(make_path(prev))
                    if pcentroid[0] < tree.x:
                        shift_subtree(prev, delta)
        shifted.append(node)
        last = node


def shift_subtree(tree, dx):
    for n in breadth_first_traversal(tree):
        n.x += dx
    return tree
