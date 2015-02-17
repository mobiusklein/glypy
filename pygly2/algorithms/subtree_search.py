import logging
from pygly2.utils import makestruct
logger = logging.getLogger(__name__)

SubtreeRecord = makestruct("SubtreeRecord", ("subtree", "include", "link_ids"))


class SubtreeMatchRecord(object):
    def __init__(self, subtree, include=None, link_ids=None):
        if include is None:
            include = []
        if link_ids is None:
            link_ids = []
        self.subtree = subtree
        self.include = map(frozenset, include)
        self.link_ids = link_ids

    def size(self):
        return len(self.include[0])

    __len__ = size

    def __repr__(self):
        return "SubtreeMatchRecord({size})".format(size=len(self))


def max_order(a):
    return max([node.order() for node in a])


def exhaustive_subtrees(a, b, min_size=2):
    '''
    A naive implementation for finding the maximum common subtree between two |Glycan|s
    '''
    index_a = {}
    max_size = min_size

    for record in exhaustive_subtree_comparison(a, b, min_size):
        if len(record) > max_size:
            index_a = {}
            max_size = len(record)
        match = index_a.get(record.include[0])
        if match is not None:
            if len(match.link_ids[0]) >= len(record.link_ids[0]) and len(match.link_ids[1]) >= len(record.link_ids[1]):
                index_a[record.include[0]] = record
        else:
            index_a[record.include[0]] = record
    return sorted(index_a.values(), key=len, reverse=True)


def exhaustive_subtree_comparison(a, b, min_size=2):
    subtrees_look_up = list()
    max_observed = min_size

    for parent_tree, parent_include, child_tree, child_include, link_ids in b.fragments(
            max_cleavages=max_order(b), structures=True):
        if len(parent_include) > min_size:
            subtrees_look_up.append(SubtreeRecord(parent_tree, parent_include, link_ids))
        if len(child_include) > min_size:
            subtrees_look_up.append(SubtreeRecord(child_tree, child_include, link_ids))
    subtrees_look_up.append(SubtreeRecord(b.clone(), list(n.id for n in b), []))

    a_index = list(a.index)
    for tree in subtrees_look_up:
        if len(a_index) == len(tree.include) and a == tree.subtree:
            yield SubtreeMatchRecord(a.clone(), [list(a_index), tree.include], [[], tree.link_ids])
            max_observed = len(a_index)
    for tree in a.fragments(max_cleavages=max_order(a), structures=True):
        for subtree in subtrees_look_up:
            if len(tree.parent_include_nodes) >= max_observed:
                if len(tree.parent_include_nodes) == len(subtree.include) and tree.parent_tree == subtree.subtree:
                    max_observed = len(subtree.include)
                    yield SubtreeMatchRecord(tree.parent_tree,
                                             [tree.parent_include_nodes, subtree.include],
                                             [tree.link_ids, subtree.link_ids])

            if len(tree.child_include_nodes) >= max_observed:
                if len(tree.child_include_nodes) == len(subtree.include) and tree.child_tree == subtree.subtree:
                    max_observed = len(subtree.include)
                    yield SubtreeMatchRecord(tree.child_tree,
                                             [tree.child_include_nodes, subtree.include],
                                             [tree.link_ids, subtree.link_ids])
