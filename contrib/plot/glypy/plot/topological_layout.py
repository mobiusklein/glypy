import numpy as np

#: clock-face coordinates assume left-horizontal growth
index_position_shift = {
    0: (0., 0.),
    2: (0.0, -1),
    3: (0.5, -1),
    4: (0.8, 0),
    5: (0.5, 1),
    6: (0.5, 1),
    8: (0, 1)
}
index_position_shift = {k: np.array(v, dtype=np.float64) for k, v in index_position_shift.items()}


def find_link(parent, child):
    for p, link in parent.links.items():
        if child.parent_linkage == p:
            if link.child == child:
                return link


def test_draw_unit_perimeter(ax):
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
    return tree
