'''
A work in progress to abstract over the messy shared dictionary scheme used in
draw_tree.DrawTreeNode
'''

import weakref
from itertools import chain

import numpy as np
from matplotlib.path import Path

try:
    from matplotlib.path import get_paths_extents
except ImportError:
    # TODO: Check why this was deprecated, beyond redundancy
    def get_paths_extents(paths):
        from matplotlib.path import _path
        from matplotlib.transforms import BBox, Affine2D

        transforms = []
        if len(paths) == 0:
            raise ValueError("No paths provided")
        return Bbox.from_extents(*_path.get_path_collection_extents(
            Affine2D(), paths, transforms, [], Affine2D()))


chain_iterable = chain.from_iterable


def ensure_iterable(obj):
    try:
        return [x for x in obj if x is not None]
    except Exception:
        if obj is None:
            return []
        return [obj]


def flatten(iterable):
    return list(chain_iterable(map(ensure_iterable, iterable)))


class GraphicalPatch(object):  # pragma: no cover
    def __init__(self, shape_patches, text_patches, parent=None, **kwargs):
        self.shape_patches = shape_patches
        self.text_patches = text_patches
        self.parent = weakref.ref(parent) if parent is not None else None
        self.options = kwargs

    def transform_shapes(self, transform):
        for shape in self.shapes():
            shape.set_transform(transform)

    def transform_text(self, transform):
        for text in self.text():
            text.set_transform(transform)

    def shapes(self):
        for shape in self.shape_patches:
            if shape is not None:
                yield shape

    def text(self):
        for t in self.text_patches:
            if t is not None:
                yield t

    def set_gid(self, tree_id, node_id):
        for i, shape in enumerate(self.shapes()):
            shape.set_gid("%s-%s-%d-shape" % (tree_id, node_id, i))
        for i, text in enumerate(self.text()):
            text.set_gid("%s-%s-%d-text" % (tree_id, node_id, i))

    def centroid(self):
        point = np.zeros(2)
        c = 0
        for shape in flatten(self.shapes()):
            path_ = shape.get_path()
            for p in path_.vertices:
                point += p
                c += 1
        return point / c

    def get_paths(self):
        result = []
        for shape in flatten(self.shapes()):
            result.append(shape.get_path())
        return result

    def get_bbox(self):
        return get_paths_extents(self.get_paths())


class MonosaccharidePatch(GraphicalPatch):
    def __init__(self, saccharide_shape=(None,), substituent_text=(None,), saccharide_label=(None,),
                 parent=None, **kwargs):
        saccharide_shape = ensure_iterable(saccharide_shape)
        substituent_text = ensure_iterable(substituent_text)
        saccharide_label = ensure_iterable(saccharide_label)
        self.saccharide_shape = flatten(saccharide_shape)
        self.substituent_text = flatten(substituent_text)
        self.saccharide_label = saccharide_label
        super(MonosaccharidePatch, self).__init__(
            (saccharide_shape),
            (substituent_text + saccharide_label),
            parent, **kwargs)

    def set_gid(self, tree_id, node_id):
        node_id = str(node_id)
        for el in self.saccharide_shape:
            el.set_gid("%s-%s-node" % (tree_id, node_id))
        for i, el in enumerate(self.substituent_text):
            el.set_gid("%s-%s-subst-%d" % (tree_id, node_id, i))
        for el in self.saccharide_label:
            el.set_gid("%s-%s-label" % (tree_id, node_id))

    def add_substituents(self, substituents):
        substituents = flatten(ensure_iterable(substituents))
        self.substituent_text.extend(substituents)
        self.text_patches.extend(substituents)

    def __iter__(self):
        yield self.saccharide_shape
        yield self.substituent_text
        yield self.saccharide_label


class EdgePatch(GraphicalPatch):  # pragma: no cover
    def __init__(self, line, parent_label=None, child_label=None, parent=None, **kwargs):
        self.line_shape = ensure_iterable(line)
        self.parent_label = ensure_iterable(parent_label)
        self.child_label = ensure_iterable(child_label)
        super(EdgePatch, self).__init__(
            line, (parent_label, child_label), parent=parent, **kwargs)


class GraphicalAggregate(list):
    def __init__(self, elements):
        list.__init__(self, elements)
