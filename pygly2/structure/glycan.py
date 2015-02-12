import logging
import itertools
logger = logging.getLogger("Glycan")
from functools import partial
from collections import deque, defaultdict, Callable

from .base import SaccharideBase
from .monosaccharide import Monosaccharide

from ..utils import make_counter, identity, StringIO
from ..composition import Composition


fragment_shift = {
    "B": Composition(O=1, H=2),
    "Y": Composition(),
    "C": Composition(),
    "Z": Composition(H=2, O=1),
}


class Glycan(SaccharideBase):
    '''
    Represents a full graph of connected |Monosaccharide| objects and their connecting bonds.

    Attributes
    ----------
    root: |Monosaccharide|
        The first monosaccharide unit of the glycan, and the reducing end if present.
    index: :class:`list`
        A list of the |Monosaccharide| instances in `self` in the order they are encountered
        by traversal by `traversal_methods[index_method]`
    link_index: :class:`list`
        A list of the |Link| connecting the |Monosaccharide| instances in `self` in the order they
        are encountered by traversal by `traversal_methods[index_method]`
    reducing_end: |int| or |None|
        The index of the reducing end on :attr:`root`.

    '''
    traversal_methods = {}

    def __init__(self, graph_root=None, reducing_end=None, index_method='dfs'):
        '''
        Constructs a new Glycan from the set of |Monosaccharide|
        rooted at `graph_root`. 

        If `reducing_end` is not |None|, it is set as the reducing end.

        If index_method is not |None|, the graph is indexed by the default search method
        given by `traversal_methods[index_method]`
        '''
        logger.debug("Creating new glycan with root: %s", graph_root)
        if graph_root is None:
            graph_root = Monosaccharide()
        self.root = graph_root
        self.index = []
        self.link_index = []
        logger.debug("Reindexing")
        if index_method is not None:
            self.reindex(index_method)
        if reducing_end is not None:
            self.root.reducing_end = reducing_end

    def reindex(self, method='dfs'):
        '''
        Traverse the graph using the function specified by ``method``. The order of 
        traversal defines the new :attr:`id` value for each |Monosaccharide|
        and |Link|.

        The order of traversal also defines the ordering of the |Monosaccharide|
        in :attr:`index` and |Link| in :attr:`link_index`.

        '''
        traversal = self._get_traversal_method(method)
        index = []
        i = 1
        for node in traversal():
            index.append(node)
        for node in index:
            node.id = i
            i += 1
        link_index = []
        for pos, link in self.iterlinks(method=method):
            link_index.append(link)

        i = 1
        for link in link_index:
            link.id = i
            i += 1

        self.index = index
        self.link_index = link_index
        return self

    def deindex(self):
        for node in self.index:
            node.id *= -1
        for link in self.link_index:
            link.id *= -1
        return self

    def __getitem__(self, ix, method='dfs'):
        '''
        Alias for :attr:`index.__getitem__`
        '''
        if not isinstance(ix, int):
            raise TypeError("Residue indices must be integers")
        if self.index is None:
            self.reindex(method=method)
        return self.index[ix]

    @property
    def root(self):
        return self._root

    @root.setter
    def root(self, value):
        self._root = value

    @property
    def reducing_end(self):
        '''
        An alias for :attr:`Monosaccharide.reducing_end` for :attr:`root`
        '''
        return self.root.reducing_end

    @reducing_end.setter
    def reducing_end(self, value):
        self.root.reducing_end = value

    def depth_first_traversal(self, from_node=None, apply_fn=identity, visited=None):
        node_stack = list([self.root if from_node is None else from_node])
        visited = set() if visited is None else visited
        while len(node_stack) > 0:
            node = node_stack.pop()
            if node.id in visited:
                continue
            visited.add(node.id)
            yield apply_fn(node)
            node_stack.extend(sorted((terminal for pos, link in node.links.items()
                              for terminal in link if terminal.id not in visited), key=Monosaccharide.order))

    # Convenience aliases and the set up the traversal_methods entry
    dfs = depth_first_traversal
    traversal_methods['dfs'] = "dfs"
    traversal_methods['depth_first_traversal'] = "dfs"

    def breadth_first_traversal(self, from_node=None, apply_fn=identity, visited=None):
        node_queue = deque([self.root if from_node is None else from_node])
        visited = set() if visited is None else visited
        while len(node_queue) > 0:
            node = node_queue.popleft()
            if node.id in visited:
                continue
            visited.add(node.id)
            yield apply_fn(node)
            node_queue.extend(sorted((terminal for pos, link in node.links.items()
                              for terminal in link if terminal.id not in visited), key=Monosaccharide.order))

    # Convenience aliases and the set up the traversal_methods entry
    bfs = breadth_first_traversal
    traversal_methods['bfs'] = "bfs"
    traversal_methods['breadth_first_traversal'] = "bfs"

    def _get_traversal_method(self, method):
        if isinstance(method, Callable):
            return partial(method, self)
        if method == 'dfs':
            return self.dfs
        elif method == 'bfs':
            return self.bfs
        traversal = self.traversal_methods.get(method, None)
        if traversal is None:
            raise AttributeError("Unknown traversal method: {}".format(method))
        traversal = getattr(self, traversal)
        return traversal

    def __iter__(self):
        return self.depth_first_traversal()

    def iternodes(self, from_node=None, method='dfs', visited=None):
        traversal = self._get_traversal_method(method)
        return traversal(from_node=from_node, visited=visited)

    def iterlinks(self, substituents=False, method='dfs', visited=None):
        traversal = self._get_traversal_method(method)
        links_visited = set()

        def links(obj):
            if substituents:
                for pos, link in obj.substituent_links.items():
                    yield (pos, link)
            for pos, link in obj.links.items():
                if link.id in links_visited:
                    continue
                links_visited.add(link.id)
                yield (pos, link)

        return itertools.chain.from_iterable(traversal(apply_fn=links, visited=visited))

    def leaves(self, method='dfs', visited=None):
        traversal = self._get_traversal_method(method)

        def no_leaves(obj):
            if len(list(obj.children())) == 0:
                yield obj

        return itertools.chain.from_iterable(traversal(apply_fn=no_leaves, visited=visited))

    def to_glycoct(self, buffer=None, close=False):
        '''
        Serialize the |Glycan| graph object into condensed GlycoCT, using
        `buffer` to store the result. If `buffer` is |None|, then the
        function will operate on a newly created :class:`~pygly2.utils.StringIO` object.

        Parameters
        ----------
        buffer: `file-like` or |None|
            The stream to write the serialized structure to. If |None|, uses an instance
            of `StringIO`
        close: |bool|
            Whether or not to close the stream in `buffer` after writing is done

        Returns
        -------
        file-like or str if ``buffer`` is :const:`None`

        '''
        is_stringio = False
        if buffer is None:
            buffer = StringIO()
            is_stringio = True

        buffer.write("RES\n")

        res_counter = make_counter()
        lin_counter = make_counter()

        # Look-ups for mapping RES nodes to objects by section index and id, respectively
        index_to_residue = {}
        residue_to_index = {}

        # Accumulator for linkage indices and mapping linkage indices to dependent RES indices
        lin_accumulator = []
        dependencies = defaultdict(dict)

        for node in (self):
            res, lin, index = node.to_glycoct(res_counter, lin_counter, complete=False)

            lin_accumulator.append((index, lin))
            residue_to_index[node.id] = index
            index_to_residue[index] = node

            for pos, link in node.links.items():
                if link.is_child(node):
                    continue
                dependencies[link.child.id][node.id] = ((lin_counter(), link))
            for line in res:
                buffer.write(line + '\n')

        buffer.write("LIN\n")
        for res_ix, links in lin_accumulator:
            for line in links:
                buffer.write(line + '\n')
            residue = index_to_residue[res_ix]
            for pos, link in residue.links.items():
                if link.is_child(residue):
                    continue
                child_res = link.child
                ix, link = dependencies[child_res.id][residue.id]
                buffer.write(link.to_glycoct(ix, res_ix, residue_to_index[child_res.id]) + "\n")

        if is_stringio:
            return buffer.getvalue()
        else:
            if close:
                buffer.close()
            return buffer

    __repr__ = to_glycoct

    def mass(self, average=False, charge=0, mass_data=None):
        '''
        Calculates the total mass of the intact graph by querying each
        node for its mass.

        Parameters
        ----------
        average: |bool|, optional, defaults to |False|
            Whether or not to use the average isotopic composition when calculating masses.
            When ``average == False``, masses are calculated using monoisotopic mass.
        charge: |int|, optional, defaults to 0
            If charge is non-zero, m/z is calculated, where m is the theoretical mass, and z is `charge`
        mass_data: |dict|, optional, defaults to |None|
            If mass_data is None, standard NIST mass and isotopic abundance data are used. Otherwise the 
            contents of mass_data are assumed to contain elemental mass and isotopic abundance information.

        See also
        --------
        :func:`pygly2.composition.composition.calculate_mass`
        '''
        return sum(node.mass(average=average, charge=charge, mass_data=mass_data) for node in self)


    def clone(self):
        clone_root = self.root.clone()
        node_stack = [(clone_root, self.root)]
        visited = set()
        while(len(node_stack) > 0):
            clone, ref = node_stack.pop()
            if ref.id in visited:
                continue
            visited.add(ref.id)
            for pos, link in ref.links.items():
                child = link.to(ref)
                if child.id in visited:
                    continue            
                clone_child = child.clone()
                clone_child.id = child.id
                link.clone(clone, clone_child)
                node_stack.append((clone_child, child))
        return Glycan(clone_root)

    def __eq__(self, other):
        if (other is None):
            return False
        traversal = zip(self.dfs(), other.dfs())
        for a, b in traversal:
            if a != b:
                return False
        return True

    def __ne__(self, other):
        return not self == other
        
    def fragments(self, kind=('B','Y'), max_cleavages=1, average=False, charge=0, mass_data=None):
        '''
        Generate carbohydrate backbone fragments from this glycan by examining the disjoint subtrees 
        created by removing one or more monosaccharide-monosaccharide bond.

        .. note::
            Does not copy the structure, though it does mutate the structure but reverses the effect so there should be no side-effects. 
            If an error should occur during iteration, there is a chance that the glycan structure may become corrupted. Efforts have been
            taken to prevent such corruption from happening, but as the library develops, new errors may not yet be anticipated.
        
        Does not produce cross-ring cleavages. 

        Parameters
        ----------

        kind: `sequence`
            Any `iterable` or `sequence` of characters corresponding to glycan fragment nomenclature
            as published by :title-reference:`Domon and Costello`
        max_cleavages: |int|
            The maximum number of bonds to break per fragment
        average: |bool|, optional, defaults to |False|
            Whether or not to use the average isotopic composition when calculating masses.
            When ``average == False``, masses are calculated using monoisotopic mass.
        charge: |int|, optional, defaults to 0
            If charge is non-zero, m/z is calculated, where m is the theoretical mass, and z is `charge`
        mass_data: |dict|, optional, defaults to |None|
            If mass_data is |None|, standard NIST mass and isotopic abundance data are used. Otherwise the 
            contents of `mass_data` are assumed to contain elemental mass and isotopic abundance information.

        See also
        --------
        :func:`pygly2.composition.composition.calculate_mass`
        '''
        for i in range(1, max_cleavages + 1):
            link_group_gen = itertools.combinations((link for ind, link in self.iterlinks()), i)
            for links in link_group_gen:
                for frag in self.break_links(links, kind, average, charge, mass_data):
                    yield frag



    def break_links(self, link_list, kind=('B','Y'), average=False, charge=0, mass_data=None):
        if len(link_list) == 0:
            raise StopIteration
        try:
            link = link_list[0]
            parent, child = link.break_link()
            link.refund()
            break_id = link.id
            parent_tree = Glycan(graph_root=parent, index_method=None)
            child_tree = Glycan(graph_root=child, index_method=None)

            if len(link_list) > 1:
                parent_frags = list(parent_tree.break_links(link_list[1:], kind=kind, average=average, charge=charge, mass_data=mass_data))
                child_frags = list(child_tree.break_links(link_list[1:], kind=kind, average=average, charge=charge, mass_data=mass_data))
                
                if 'Y' in kind:
                    offset = fragment_shift['Y'].calc_mass(average=average, charge=charge, mass_data=mass_data)
                    for ion_type, link_ids, include, mass in parent_frags:
                        yield ion_type + 'Y',  link_ids + [break_id], include, mass - offset
                if 'Z' in kind:
                    offset = fragment_shift['Z'].calc_mass(average=average, charge=charge, mass_data=mass_data)
                    for ion_type, link_ids, include, mass in parent_frags:
                        yield ion_type + 'Z', link_ids + [break_id], include, mass - offset
                if 'B' in kind:
                    offset = fragment_shift['B'].calc_mass(average=average, charge=charge, mass_data=mass_data)
                    for ion_type, link_ids, include, mass in child_frags:
                        yield ion_type + 'B', link_ids + [break_id], include, mass - offset
                if 'C' in kind:
                    offset = fragment_shift['C'].calc_mass(average=average, charge=charge, mass_data=mass_data)
                    for ion_type, link_ids, include, mass in child_frags:
                        yield ion_type + 'C', include, link_ids + [break_id], mass - offset
            else:
                parent_include  = [n.id for n in parent_tree]
                child_include  = [n.id for n in child_tree]

                if 'Y' in kind:
                    offset = fragment_shift['Y'].calc_mass(
                            average=average, charge=charge, mass_data=mass_data)                 
                    y_mass = parent_tree.mass(
                                average=average, charge=charge, mass_data=mass_data) - offset
                    yield ('Y', [break_id], parent_include, y_mass)
                if 'B' in kind:
                    offset = fragment_shift['B'].calc_mass(
                            average=average, charge=charge, mass_data=mass_data)
                    b_mass = child_tree.mass(
                        average=average, charge=charge, mass_data=mass_data) - offset
                    yield 'B', [break_id], child_include, b_mass
                if "C" in kind:
                    offset = fragment_shift['C'].calc_mass(
                            average=average, charge=charge, mass_data=mass_data)
                    c_mass = child_tree.mass(
                        average=average, charge=charge, mass_data=mass_data) - offset
                    yield 'C', [break_id], child_include, c_mass
                if 'Z' in kind:
                    offset = fragment_shift['Z'].calc_mass(
                            average=average, charge=charge, mass_data=mass_data)
                    z_mass = parent_tree.mass(
                            average=average, charge=charge, mass_data=mass_data) - offset
                    yield ('Z', [break_id], parent_include, z_mass)
        except Exception, e:
            print(e)
        finally:
            link.apply()

