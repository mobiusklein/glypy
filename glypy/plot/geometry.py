from matplotlib import patches as mpatches
from matplotlib import path as mpath
from matplotlib import transforms as mtransform
import numpy as np

from glypy.utils import groupby


Affine2D = mtransform.Affine2D


def breadth_first_traversal(tree, visited=None):
    if visited is None:
        visited = set()
    if tree.id not in visited:
        visited.add(tree.id)
        yield tree
        for child in tree:
            for descend in breadth_first_traversal(child, visited=visited):
                yield descend


def _build_path_with_radial_bounding_box(node, implicit_radius=0.3):
    pathing = []
    last = [node.x, node.y]
    # build perimeter
    for node in breadth_first_traversal(node):
        pathing.append([node.x, node.y])
        pathing.append(last)
        last = [node.x, node.y]
    codes = [mpath.Path.MOVETO] + [mpath.Path.LINETO for i in pathing[1:]]
    path = mpath.Path(pathing, codes)
    for node in breadth_first_traversal(node):
        center = (node.x, node.y)
        path = path.make_compound_path(path, mpath.Path.circle(center, implicit_radius))
    # for node in breadth_first_traversal(node):
    #     for part in node.patch_node.shapes():
    #         path = path.make_compound_path(path, part.get_path())
    return path


def _build_path_with_centroids(node, **kwargs):
    pathing = []
    last = [node.x, node.y]
    for node in breadth_first_traversal(node):
        pathing.append([node.x, node.y])
        pathing.append(last)
        last = [node.x, node.y]
    codes = [mpath.Path.MOVETO] + [mpath.Path.LINETO for i in pathing[1:]]
    return mpath.Path(pathing, codes)


def _build_path_with_circles(node, implicit_radius=0.3):
    points = []
    for node in breadth_first_traversal(node):
        points.append(mpath.Path.circle((node.x, node.y), implicit_radius))
    return mpath.Path.make_compound_path(*points)


def make_path(node):
    return _build_path_with_circles(node)


def centroid(path):
    point = np.zeros_like(path.vertices[0])
    c = 0.
    for p in path.vertices:
        point += p
        c += 1
    return point / c


# def shift_path_2d(path, hits):
#     centers = map(centroid, hits)
#     stacked = np.vstack(centers)
#     min_x = stacked[:, 0].min()
#     max_x = stacked[:, 0].max()
#     min_y = stacked[:, 1].min()
#     max_y = stacked[:, 1].max()
#     new_center = np.array([(min_x + max_x)/2., (min_y + max_y)/2.]) * 0.2
#     current_center = centroid(path)
#     distance = current_center - new_center
#     path = path.transformed(Affine2D().translate(*distance))
#     return path

# def shift_path_away(path, hit):
#     center = centroid(path)
#     distance = ((centroid(hit) - center) * 0.2)
#     return path.transformed(Affine2D().translate(*distance))



# def find_overlaps(path, objs):
#     hits = []
#     for obj in objs:
#         if path.intersects_path(obj):
#             hits.append(obj)
#     return hits

# def minimize_overlap(path, objs, max_distance=1.0):
#     hits = find_overlaps(path, objs)
#     # We hit two or more objects
#     if len(hits) > 1:
#         path = shift_path_2d(path, hits)
#     elif len(hits) == 1:
#         path = shift_path_away(path, hits[0])
#     again_hits = find_overlaps(path, objs)
#     return path
