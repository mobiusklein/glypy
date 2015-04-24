import numpy as np

from pygly2.structure.constants import Anomer


index_position_shift = {
    1: (-1, 0),
    2: (-0.5, -1),
    3: (0.5, -1),
    4: (1, 0),
    5: (0.5, 0.5),
    6: (0.5, 1)
}


def square_perimeter(cx, cy, index_position, scale=0.1):
    sx, sy = index_position_shift[index_position]
    sx *= scale
    sy *= scale
    return cx + sx, cy + sy


def layout_node(node, cx, cy, anomer=Anomer.alpha, parent_position=3, child_position=1):
    pass
