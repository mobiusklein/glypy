from matplotlib import patches as mpatches
from matplotlib import path as mpath
from matplotlib import transforms as mtransform
import numpy as np

from glypy.utils import groupby


Affine2D = mtransform.Affine2D


class TreeLayoutBase(object):
    def __init__(self, root, **kwargs):
        self.root = root

    def layout(self, **kwargs):
        self.before_layout(**kwargs)
        self.layout_tree(**kwargs)
        self.after_layout(**kwargs)

    def traverse(self):
        return self.root.traverse()

    def layout_tree(self, **kwargs):
        raise NotImplementedError()

    def before_layout(self, **kwargs):
        pass

    def after_layout(self, **kwargs):
        self.root.fix_special_cases()

    def transform(self, transform):
        for node in self.traverse():
            node.x, node.y = transform.transform([node.x, node.y])


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


def interpolate(x0, y0, x1, y1, x):
    y = y0 + (x - x0) * ((y1 - y0) / (x1 - x0))
    return y


def make_path(node, implicit_radius=0.3):
    return _build_path_with_circles(node, implicit_radius)


def centroid(path):
    point = np.zeros_like(path.vertices[0])
    c = 0.
    for p in path.vertices:
        point += p
        c += 1
    return point / c


def bounding_box(path, padding=0):
    if not padding:
        padding = 0
    xmin = float('inf')
    xmax = -float('inf')
    ymin = float('inf')
    ymax = -float('inf')
    for p in path.vertices:
        xmin = min(xmin, p[0])
        xmax = max(xmax, p[0])
        ymin = min(ymin, p[1])
        ymax = max(ymax, p[1])
    xmin -= padding
    xmax += padding
    ymin -= padding
    ymax += padding

    return mtransform.Bbox.from_extents(xmin, ymin, xmax, ymax)
    # return xmin, xmax, ymin, ymax


class TransformProxy(object):
    def __init__(self, base_transform, owner_update):
        self.base_transform = base_transform
        self.owner_update = owner_update

    def _copy_transform(self):
        return self.base_transform.frozen()

    def scale(self, *args, **kwargs):
        transform = self._copy_transform()
        transform.scale(*args, **kwargs)
        self.owner_update(transform)
        self.base_transform = transform
        return self

    def rotate(self, *args, **kwargs):
        transform = self._copy_transform()
        transform.rotate(*args, **kwargs)
        self.owner_update(transform)
        self.base_transform = transform
        return self

    def rotate_deg(self, *args, **kwargs):
        transform = self._copy_transform()
        transform.rotate_deg(*args, **kwargs)
        self.owner_update(transform)
        self.base_transform = transform
        return self

    def rotate_around(self, *args, **kwargs):
        transform = self._copy_transform()
        transform.rotate_around(*args, **kwargs)
        self.owner_update(transform)
        self.base_transform = transform
        return self

    def rotate_deg_around(self, *args, **kwargs):
        transform = self._copy_transform()
        transform.rotate_deg_around(*args, **kwargs)
        self.owner_update(transform)
        self.base_transform = transform
        return self

    def __getattr__(self, name):
        if name == 'base_transform':
            raise AttributeError("Accessing Base Transform Too Early")
        return getattr(self.base_transform, name)


def l2_distance(a, b):
    return sum([(ai - bi) ** 2 for ai, bi in zip(a, b)])


def l1_distance(a, b):
    return sum([abs(ai - bi) for ai, bi in zip(a, b)])


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
