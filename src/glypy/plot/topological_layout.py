from .geometry import breadth_first_traversal, centroid, make_path, TreeLayoutBase
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


class TopologicalTreeLayout(TreeLayoutBase):

    def __init__(self, root, **kwargs):
        super(TopologicalTreeLayout, self).__init__(root, **kwargs)
        self.linkage_index = dict()

    def before_layout(self, **kwargs):
        for node in self.traverse():
            node.mask_special_cases = False

    def after_layout(self, **kwargs):
        """Apply any cleanup after applying layout
        """
        pass

    def layout_tree(self, **kwargs):
        visited = set()
        iterations = kwargs.get("iterations", 25)
        self.layout_node(self.root, 0, 0, 0, visited)
        for i in range(int(iterations)):
            delta = self.adjust_subtrees(self.root)
            if delta == 0:
                break

    def position_relative_to_parent(self, cx, cy, index_position, scale=1.0):
        sx, sy = index_position_shift[index_position]
        sx *= scale
        sy *= scale
        return cx - sx, cy - sy

    def layout_node(self, node, cx, cy, parent_position=3, visited=None):
        if visited is None:
            visited = set()
        if node.id in visited:
            return node
        self.linkage_index[node.id] = parent_position
        node.x, node.y = self.position_relative_to_parent(cx, cy, parent_position)
        node.mask_special_cases = False
        # default positions for unknown linkages
        default_positions = ([3, 5, 4, 6, 2, 8])
        # when there is only one child node, and that node is unknown,
        # prefer a straight line instead of a slanted line.
        if len(node.tree.children()) == 1:
            default_positions = [4]
        deferred_children = []
        for child_node, pos_child in zip(node.children, sorted(node.tree.children())):
            pos = pos_child[0]
            if pos != -1:
                if pos in default_positions:
                    default_positions.remove(pos)
            else:
                deferred_children.append(child_node)
                continue
            self.layout_node(child_node, node.x, node.y, pos)
        if len(deferred_children) > len(default_positions):
            raise ValueError("Out of child positions")
        for child_node, pos in zip(deferred_children, default_positions):
            self.layout_node(child_node, node.x, node.y, pos)
        return node

    def adjust_subtrees(self, tree):
        next_layer = tree.children
        if len(next_layer) == 0:
            return 0

        total_shift = 0

        # layout later layers first in a depth-first traversal
        for node in next_layer:
            total_shift += self.adjust_subtrees(node)

        if len(next_layer) < 2:
            return total_shift

        next_layer = sorted(next_layer, key=lambda x: x.x, reverse=False)

        prior_nodes = [next_layer[0]]
        shifted = [next_layer[0]]
        for node in next_layer[1:]:
            for last in prior_nodes:
                # This creates a circle around each node and reports on overlaps
                # of nodes, but doesn't care about edges. The path extraction
                # function should try to cover edges too.
                lpath = make_path(last, 0.3)
                npath = make_path(node, 0.3)
                delta = 0
                # if the covering paths intersect
                if lpath.intersects_path(npath):
                    lcentroid = centroid(lpath)
                    ncentroid = centroid(npath)
                    lxmin_lymin, lxmax_lymax = lpath.get_extents().get_points()
                    nxmin_nymin, nxmax_nymax = npath.get_extents().get_points()
                    # the node is in a special position which should not be shifted
                    # according to examples
                    # if last.x == tree.x or self.linkage_index[last.id] == 4:
                    if self.linkage_index[node.id] == 4:
                        delta = (nxmax_nymax[0] - ncentroid[0]) / -2
                        self.shift_subtree(last, dx=delta)
                    elif ncentroid[0] > lcentroid[0]:
                        delta = (nxmin_nymin[0] - ncentroid[0]) / -2
                        self.shift_subtree(node, dx=delta)
                    elif ncentroid[0] < lcentroid[0]:
                        delta = (lxmax_lymax[0] - lcentroid[0]) / -2
                        self.shift_subtree(last, dx=delta)
                    total_shift += abs(delta)
            shifted.append(node)
            prior_nodes.append(node)
        return total_shift

    def shift_subtree(self, tree, dx=0, dy=0):
        for n in breadth_first_traversal(tree):
            n.x += dx
            n.y += dy
        return tree


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
    for child_node, pos_child in zip(node.children, sorted(node.tree.children())):
        pos = pos_child[0]
        layout_node(child_node, node.x, node.y, pos)
    return node


def layout(tree, visited=None):
    if visited is None:
        visited = set()
    layout_node(tree, 0, 0, 0, visited)
    for i in range(10):
        delta = spread(tree)
        if delta == 0:
            break
    return tree


def spread(tree):
    next_layer = tree.children
    if len(next_layer) == 0:
        return 0

    total_shift = 0

    # layout later layers first in a depth-first traversal
    for node in next_layer:
        total_shift += spread(node)

    if len(next_layer) < 2:
        return total_shift

    # layout this layer in x-ascending order
    next_layer = sorted(next_layer, key=lambda x: x.x)

    prior_nodes = [next_layer[0]]
    shifted = [next_layer[0]]
    for node in next_layer[1:]:
        for last in prior_nodes:
            lpath = make_path(last)
            npath = make_path(node)
            delta = 0
            if lpath.intersects_path(npath):
                lcentroid = centroid(lpath)
                ncentroid = centroid(npath)
                lxmin_lymin, lxmax_lymax = lpath.get_extents().get_points()
                nxmin_nymin, nxmax_nymax = npath.get_extents().get_points()
                if last.x == tree.x:
                    delta = (nxmax_nymax[0] - ncentroid[0]) / 2
                    shift_subtree(node, delta)
                elif ncentroid[0] > lcentroid[0]:
                    delta = (nxmin_nymin[0] - ncentroid[0]) / -2
                    shift_subtree(node, delta)
                elif ncentroid[0] < lcentroid[0]:
                    delta = (lxmax_lymax - lcentroid[0]) / -2
                    shift_subtree(last, delta)
                total_shift += abs(delta)
        shifted.append(node)
        prior_nodes.append(node)
    return total_shift


def shift_subtree(tree, dx):
    for n in breadth_first_traversal(tree):
        n.x += dx
    return tree
