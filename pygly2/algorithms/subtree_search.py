import logging
from pygly2.utils import make_struct
logger = logging.getLogger(__name__)

SubtreeRecord = make_struct("SubtreeRecord", ("subtree", "include", "link_ids"))


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
            max_cleavages=b.order(), structures=True):
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
    for tree in a.fragments(max_cleavages=a.order(), structures=True):
        if len(tree.parent_include_nodes) >= max_observed:
            for subtree in subtrees_look_up:
                if len(tree.parent_include_nodes) == len(subtree.include) and tree.parent_tree == subtree.subtree:
                    max_observed = len(subtree.include)
                    yield SubtreeMatchRecord(tree.parent_tree,
                                             [tree.parent_include_nodes, subtree.include],
                                             [tree.link_ids, subtree.link_ids])

        if len(tree.child_include_nodes) >= max_observed:
            for subtree in subtrees_look_up:
                if len(tree.child_include_nodes) == len(subtree.include) and tree.child_tree == subtree.subtree:
                    max_observed = len(subtree.include)
                    yield SubtreeMatchRecord(tree.child_tree,
                                             [tree.child_include_nodes, subtree.include],
                                             [tree.link_ids, subtree.link_ids])


def substituent_inclusion(self, other):
    taken_b = set()
    b_substituents = list(other.substituents())
    cntr = 0
    for a_pos, a_substituent in self.substituents():
        matched = False
        cntr += 1
        for b_pos, b_substituent in b_substituents:
            if b_pos in taken_b:
                continue
            if b_substituent == a_substituent:
                matched = True
                taken_b.add(b_pos)
                break
        if not matched and cntr > 0:
            return False
    return True


def topological_inclusion(self, other, substituents=True):
    '''
    Performs equality testing between two monosaccharides where
    the exact ordering of child links does not have match between
    the input |Monosaccharide|s, so long as an `a` is included in `b`

    Returns
    -------
    |bool|
    '''
    if self._flat_equality(other, lengths=False) and (not substituents or substituent_inclusion(self, other)):
        taken_b = set()
        for a_pos, a_child in self.children():
            matched = False
            for b_pos, b_child in other.children():
                if (b_pos, b_child.id) in taken_b:
                    continue
                if topological_inclusion(a_child, b_child, substituents=substituents):
                    matched = True
                    taken_b.add((b_pos, b_child.id))
                    break
            if not matched and len(list(self.children())) > 0:
                return False
        return True
    return False


def subtree_of(subtree, tree):
    for node in tree:
        if topological_inclusion(subtree.root, node):
            return node.id
    return None


class BalancedSequence(object):
    def __init__(self, tree):
        self.sequence = []
