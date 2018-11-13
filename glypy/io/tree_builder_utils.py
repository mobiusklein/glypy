from glypy.utils import enum, root as rootp, uid
from glypy.structure import monosaccharide


class StructurePrecisionEnum(enum.Enum):
    unknown = -1
    ranging = 1
    exact = 2


class AbstractGraphEntryEnum(enum.Enum):
    parent = 1
    child = 2
    internal = 3


def try_int(v):
    try:
        return int(v)
    except Exception:
        return None


def decorate_tree(tree, decorate_value):
    '''
    Transform each node in the tree's id value
    into a tuple of `(decorate_value, node.id)`
    to use common ids across identical graphs but
    discriminable between attached copies.

    Parameters
    ----------
    tree: Glycan
    decorate_value: int

    Returns
    -------
    Glycan
    '''
    for node in tuple(tree):
        node.id = (decorate_value, node.id)
        for i, substituent in node.substituents():
            substituent.id = (decorate_value, substituent.id)
    return tree


def decorated_value(tree):  # pragma: no cover
    '''
    Get the decorating value from a tree's root node.id.

    Returns `None` if the tree is undecorated

    Parameters
    ----------
    tree: Glycan

    Returns
    -------
    int or None
    '''
    try:
        return next(iter(tree.root.id))
    except Exception:
        return None


def get_decorated(tree, id):  # pragma: no cover
    '''
    As :meth:`Glycan.get`, but with awareness of decorated
    node.id attributes. Will fall back to the normal get method
    if the tree is undecorated.

    Parameters
    ----------
    tree: Glycan
    id: int

    Returns
    -------
    Monosaccharide
    '''
    d = decorated_value(tree)
    if d is None:
        return tree.get(id)
    else:
        return tree.get((d, id))


def undecorate_tree(tree):
    '''
    Remove decoration from a tree and reindex its nodes.

    Depends upon  :attr:`Glycan.index` to find nodes

    Parameters
    ----------
    tree: Glycan

    Returns
    -------
    Glycan
    '''
    i = 1
    for node in list(tree):
        node.id = i
        i += 1
        for j, substituent in node.substituents():
            substituent.id = i
            i += 1

    # tree.reindex()
    return tree


def find_root(tree):
    root = rootp(tree)

    while True:
        parents = root.parents()
        if not parents:
            break
        root = parents[0][1]
    return root


class NodeCollection(object):

    @classmethod
    def from_node(cls, root):
        nodes = []
        for node in monosaccharide.traverse(root):
            nodes.append(node)
            try:
                for _, subst in node.substituents():
                    nodes.append(subst)
            except AttributeError:
                pass
        return cls(nodes)

    def __init__(self, nodes):
        self.nodes = list(nodes)
        self.root = self.nodes[0]

    def __iter__(self):
        return monosaccharide.traverse(self.root)

    def deindex(self):
        base = uid()
        for i, node in enumerate(self.nodes):
            node.id += base
            node.id *= -1
        return self

    def reindex(self):
        for i, node in enumerate(self.nodes):
            node.id = i + 1
        return self

    def clone(self, *args, **kwargs):
        root = monosaccharide.graph_clone(self.root)
        return self.from_node(root)

    def get(self, key):
        for node in self.nodes:
            if node.id == key:
                return node
        raise IndexError(key)

    def __getitem__(self, i):
        return self.get(i)

    def __root__(self):
        return self.root
