import operator
import random
import logging
import itertools
from functools import partial
from collections import deque, defaultdict, namedtuple, Callable

from .base import SaccharideBase
from .monosaccharide import Monosaccharide
from ..utils import make_counter, identity, StringIO, chrinc
from ..composition import Composition

logger = logging.getLogger("Glycan")

fragment_shift = {
    "B": Composition(O=1, H=2),
    "Y": Composition(),
    "C": Composition(),
    "Z": Composition(H=2, O=1),
}

fragment_direction = {
    "B": -1,
    "C": -1,
    "Y": 1,
    "Z": 1
}

MAIN_BRANCH_SYM = '-'

Fragment = namedtuple("Fragment", ("kind", "link_ids", "included_nodes", "mass"))
DisjointTrees = namedtuple("DisjointTrees", ("parent_tree", "parent_include_nodes",
                                             "child_tree", "child_include_nodes", "link_ids"))


class Glycan(SaccharideBase):
    '''
    Represents a full graph of connected |Monosaccharide| objects and their connecting bonds.

    Attributes
    ----------
    root: |Monosaccharide|
        The first monosaccharide unit of the glycan, and the reducing end if present.
    index: |list|
        A list of the |Monosaccharide| instances in `self` in the order they are encountered
        by traversal by `traversal_methods[index_method]`
    link_index: |list|
        A list of the |Link| connecting the |Monosaccharide| instances in `self` in the order they
        are encountered by traversal by `traversal_methods[index_method]`
    reducing_end: |int| or |None|
        The index of the reducing end on :attr:`root`.
    branch_lengths: |dict|
        A dictionary mapping branch symbols to their lengths
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
        if graph_root is None:
            graph_root = Monosaccharide()
        self.root = graph_root
        self.index = []
        self.link_index = []
        self.branch_lengths = {}
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

        self.label_branches()

        return self

    def deindex(self):
        '''
        When combining two Glycan structures, very often their component ids will
        overlap, making it impossible to differentiate between a cycle and the new
        graph. This function mangles all of the node and link ids so that they are
        distinct from the pre-existing nodes.
        '''
        base = random.randint(1, 100000)
        for node in self.index:
            node.id += base
            node.id *= -1
        for link in self.link_index:
            link.id += base
            link.id *= -1
        return self

    def reroot(self):
        self.root = sorted(self, key=operator.attrgetter('id'))[0]
        self.reindex()
        return self

    def __getitem__(self, ix, method='dfs'):
        '''
        Alias for :attr:`index.__getitem__`
        '''
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
        '''
        Make a depth-first traversal of the glycan graph. Children are explored in descending bond-order.

        This is the default traversal method for all |Glycan|s. :meth:`~.dfs` is an alias of this method.
        Both names can be used to specify this strategy to :meth:`~._get_traversal_method`.

        Parameters
        ----------
        from_node: |None| or |Monosaccharide|
            If `from_node` is |None|, then traversal starts from the root node. Otherwise it begins
            from the given node.
        apply_fn: `function`
            A function applied to each node on arrival. If this function returns a non-None value,
            the result is yielded from the generator, otherwise it is ignored. Defaults to :func:`.identity`
        visited: :class:`set` or |None|
            A :class:`set` of node ID values to ignore. If |None|, defaults to the empty `set`

        Returns
        -------
        generator

        See also
        --------
        :meth:`~.breadth_first_traversal`
        '''
        node_stack = list([self.root if from_node is None else from_node])
        visited = set() if visited is None else visited
        while len(node_stack) > 0:
            node = node_stack.pop()
            visited.add(node.id)
            res = apply_fn(node)
            if res is not None:
                yield res
            node_stack.extend(sorted((terminal for pos, link in node.links.items()
                                      for terminal in link if terminal.id not in visited), key=Monosaccharide.order))

    # Convenience aliases and the set up the traversal_methods entry
    dfs = depth_first_traversal
    traversal_methods['dfs'] = "dfs"
    traversal_methods['depth_first_traversal'] = "dfs"

    def breadth_first_traversal(self, from_node=None, apply_fn=identity, visited=None):
        '''
        Make a breadth-first traversal of the glycan graph. Children are explored in descending bond-order.

        :meth:`~.bfs` is an alias of this method.
        Both names can be used to specify this strategy to :meth:`~._get_traversal_method`.

        Parameters
        ----------
        from_node: |None| or |Monosaccharide|
            If `from_node` is |None|, then traversal starts from the root node. Otherwise it begins
            from the given node.
        apply_fn: `function`
            A function applied to each node on arrival. If this function returns a non-None value,
            the result is yielded from the generator, otherwise it is ignored. Defaults to :func:`.identity`
        visited: :class:`set` or |None|
            A :class:`set` of node ID values to ignore. If |None|, defaults to the empty `set`

        Returns
        -------
        generator

        See also
        --------
        :meth:`~.depth_first_traversal`
        '''
        node_queue = deque([self.root if from_node is None else from_node])
        visited = set() if visited is None else visited
        while len(node_queue) > 0:
            node = node_queue.popleft()
            visited.add(node.id)
            res = apply_fn(node)
            if res is not None:
                yield res
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
        return self.iternodes()

    def iternodes(self, from_node=None, apply_fn=identity, method='dfs', visited=None):
        '''
        Generic iterator over nodes. :meth:`Glycan.__iter__` is an alias of this method

        Parameters
        ----------
        from_node: |None| or |Monosaccharide|
            If `from_node` is |None|, then traversal starts from the root node. Otherwise it begins
            from the given node.
        apply_fn: `function`
            A function applied to each node on arrival. If this function returns a non-None value,
            the result is yielded from the generator, otherwise it is ignored. Defaults to :func:`.identity`
        method: |str| or `function`
            Traversal method to use. See :meth:`._get_traversal_method`
        visited: :class:`set` or |None|
            A :class:`set` of node ID values to ignore. If |None|, defaults to the empty `set`

        Returns
        -------
        generator

        See also
        --------
        :meth:`~.depth_first_traversal`
        :meth:`~.breadth_first_traversal`
        :meth:`~._get_traversal_method`
        '''
        traversal = self._get_traversal_method(method)
        return traversal(from_node=from_node, apply_fn=apply_fn, visited=visited)

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

    def leaves(self, bidirectional=False, method='dfs', visited=None):
        traversal = self._get_traversal_method(method)
        if bidirectional:
            def is_leaf(obj):
                if len(obj.links) == 1:
                    yield obj
        else:
            def is_leaf(obj):
                if len(list(obj.children())) == 0:
                    yield obj

        return itertools.chain.from_iterable(traversal(apply_fn=is_leaf, visited=visited))

    def label_branches(self):
        '''
        Labels each branch point with an alphabetical symbol. Also computes and stores each branch's
        length and stores it in :attr:`branch_lengths`
        '''
        last_branch_label = MAIN_BRANCH_SYM
        self.branch_lengths = {MAIN_BRANCH_SYM: 0}

        def get_parent_link(node):
            try:
                return node.links[node.parents().next()[0]][0].label[0]
            except StopIteration:
                return MAIN_BRANCH_SYM

        for node in self:
            links = []
            for pos, link in node.links.items():
                if link.is_child(node):
                    continue
                links.append(link)
            if len(links) == 1:
                label_key = get_parent_link(node)
                self.branch_lengths[label_key] += 1
                label = "{}{}".format(label_key, self.branch_lengths[label_key])
                links[0].label = label
            else:
                last_label_key = label_key = get_parent_link(node)
                count = self.branch_lengths[last_label_key]
                for link in links:
                    last_branch_label = chrinc(last_branch_label) if last_branch_label != MAIN_BRANCH_SYM else 'a'
                    new_label_key = last_branch_label
                    self.branch_lengths[new_label_key] = count + 1
                    label = "{}{}".format(new_label_key, self.branch_lengths[new_label_key])
                    link.label = label
        self.branch_lengths["-"] = max(self.branch_lengths.values())

    def count_branches(self):
        '''
        Count the number of branches in the Glycan tree

        Returns
        -------
        |int|
        '''
        count = 0
        for node in self:
            if len(node.links) > 2:
                count += 2 if count == 0 else 1
        return count

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

        # Look-ups for mapping RES nodes to objects by section index and id,
        # respectively
        index_to_residue = {}
        residue_to_index = {}

        # Accumulator for linkage indices and mapping linkage indices to
        # dependent RES indices
        lin_accumulator = []
        dependencies = defaultdict(dict)

        for node in (self):
            res, lin, index = node.to_glycoct(
                res_counter, lin_counter, complete=False)

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
                buffer.write(
                    link.to_glycoct(ix, res_ix, residue_to_index[child_res.id]) + "\n")

        if is_stringio:
            return buffer.getvalue()
        else:  # pragma: no cover
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

    def clone(self, index_method='dfs'):
        '''
        Create a copy of `self`, indexed using `index_method`, a *traversal method*  or |None|.

        Returns
        -------
        |Glycan|
        '''
        clone_root = self.root.clone()
        clone_root.id = self.root.id
        node_stack = [(clone_root, self.root)]
        visited = set()
        while(len(node_stack) > 0):
            clone, ref = node_stack.pop()
            if ref.id in visited:
                continue
            visited.add(ref.id)
            links = sorted([link for pos, link in ref.links.items()], key=lambda x: x[ref].order())
            for link in links:
                terminal = link.to(ref)
                if terminal.id in visited:
                    continue
                clone_terminal = terminal.clone()
                clone_terminal.id = terminal.id
                if link.is_child(terminal):
                    link.clone(clone, clone_terminal)
                else:
                    link.clone(clone_terminal, clone)

                node_stack.append((clone_terminal, terminal))
        return Glycan(clone_root, index_method=index_method)

    def __eq__(self, other):
        '''
        Two glycans are considered equal if they are topologically equal.

        Parameters
        ----------
        self, other: Glycan

        Returns
        -------
        |bool|

        See also
        --------
        :func:`~pygly2.structure.Monosaccharide.topological_equality`
        '''
        if other is None:
            return False
        elif not isinstance(other, Glycan):
            return False
        return self.root.topological_equality(other.root)

    def __ne__(self, other):
        return not self == other

    def fragments(self, kind=('B', 'Y'), max_cleavages=1, average=False, charge=0, mass_data=None,
                  structures=False, min_cleavages=1):
        '''
        Generate carbohydrate backbone fragments from this glycan by examining the disjoint subtrees
        created by removing one or more monosaccharide-monosaccharide bond.

        .. note::
            While generating fragments, the glycan structure is being permuted. All of the
            changes being made are reversed during the generation process, and the glycan is
            returned to the same state it was in when :meth:`~.fragments` was called by the end
            of the generator. Do not attempt to use the glycan structure for other things while
            fragmenting it. If you must, copy it first with :meth:`~.clone`.

        Does not produce cross-ring cleavages.

        Parameters
        ----------

        kind: `sequence`
            Any `iterable` or `sequence` of characters corresponding to B/C/Y/Z
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
        structures: |bool| optional, defaults to |False|
            If `structures` is |True|, then instead of yielding fragment masses, yield disjoint subtrees formed by
            breaking a glycosidic bond.
        See also
        --------
        :func:`pygly2.composition.composition.calculate_mass`
        '''
        results_container = DisjointTrees if structures else Fragment
        for i in range(min_cleavages, max_cleavages + 1):
            for frag in self.break_links(i, kind, average, charge, mass_data, structures):
                    yield results_container(*frag)

    def break_links(self, n_links=0, kind=('B', 'Y'), average=False, charge=0, mass_data=None, structures=False):
        if n_links < 0:  # pragma: no cover
            raise ValueError("Cannot break a negative number of Links")
        n_links -= 1
        for pos, link in self.iterlinks():
            try:
                parent, child = link.break_link(refund=True)
                break_id = link.id
                logger.debug("Breaking %d", break_id)
                parent_tree = Glycan(graph_root=parent, index_method=None)
                child_tree = Glycan(graph_root=child, index_method=None)
                if structures:
                    if n_links > 0:
                        logger.debug("%d more links to break", n_links)

                        logger.debug("Breaking parent tree which contains %r", [
                                     l.id for p, l in parent_tree.iterlinks()])
                        for p_parent, parent_include, p_child, child_include, p_link_ids in parent_tree.break_links(
                                n_links, kind=kind, structures=structures):
                            logger.debug(
                                "Received parent tree from %r", p_link_ids)
                            yield p_parent, parent_include, p_child, child_include, p_link_ids + [break_id]

                        logger.debug("Breaking child tree which contains %r", [
                                     l.id for p, l in child_tree.iterlinks()])
                        for c_parent, parent_include, c_child, child_include, c_link_ids in child_tree.break_links(
                                n_links, kind=kind, structures=structures):
                            logger.debug(
                                "Received child tree from %r", c_link_ids)
                            yield c_parent, parent_include, c_child, child_include, c_link_ids + [break_id]
                    else:
                        parent_include = [n.id for n in parent_tree.dfs()]
                        child_include = [n.id for n in child_tree.dfs()]

                        # Copy the trees so that when the Link object is unmasked the copies returned
                        # to the user are not suddenly joined again.
                        parent_clone = parent_tree.clone(index_method=None).reroot()
                        child_clone = child_tree.clone(index_method=None).reroot()
                        yield (parent_clone, parent_include, child_clone, child_include, [break_id])

                elif n_links > 0:
                    parent_frags = list(parent_tree.break_links(
                        n_links, kind=kind, average=average, charge=charge, mass_data=mass_data))
                    child_frags = list(child_tree.break_links(
                        n_links, kind=kind, average=average, charge=charge, mass_data=mass_data))

                    if 'Y' in kind:
                        offset = fragment_shift['Y'].calc_mass(
                            average=average, charge=charge, mass_data=mass_data)
                        for ion_type, link_ids, include, mass in parent_frags:
                            yield ion_type + 'Y',  link_ids + [break_id], include, mass - offset
                    if 'Z' in kind:
                        offset = fragment_shift['Z'].calc_mass(
                            average=average, charge=charge, mass_data=mass_data)
                        for ion_type, link_ids, include, mass in parent_frags:
                            yield ion_type + 'Z', link_ids + [break_id], include, mass - offset
                    if 'B' in kind:
                        offset = fragment_shift['B'].calc_mass(
                            average=average, charge=charge, mass_data=mass_data)
                        for ion_type, link_ids, include, mass in child_frags:
                            yield ion_type + 'B', link_ids + [break_id], include, mass - offset
                    if 'C' in kind:
                        offset = fragment_shift['C'].calc_mass(
                            average=average, charge=charge, mass_data=mass_data)
                        for ion_type, link_ids, include, mass in child_frags:
                            yield ion_type + 'C', include, link_ids + [break_id], mass - offset

                else:
                    parent_include = [n.id for n in parent_tree]
                    child_include = [n.id for n in child_tree]

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

            except Exception, e:  # pragma: no cover
                link.apply()
                raise e
            finally:
                # Unmask the Link and apply its composition shifts
                logger.debug("Reapplying %d", break_id)
                link.apply()

    def name_fragment(self, fragment):
        '''
        Attempt to assign a full name to a fragment based on the branch and position relative to
        the reducing end along side B/C/Y/Z, according to :title-reference:`Domon and Costello`
        '''
        if fragment.kind not in fragment_direction:
            raise ValueError("Cannot determine fragment orientation, {}".format(fragment))
        if fragment_direction[fragment.kind] > 0:
            link = self.link_index[fragment.link_ids[0] - 1]
            label = link.label
            fragment_name = "{}{}".format(fragment.kind, label.replace(MAIN_BRANCH_SYM, ""))
        else:
            link = self.link_index[fragment.link_ids[0] - 1]
            label = link.label
            label_key = label[0]
            distance = int(label[1:])
            inverted_distance = self.branch_lengths[label_key] - (distance - 1)
            fragment_name = "{}{}{}".format(fragment.kind, label_key.replace(MAIN_BRANCH_SYM, ""), inverted_distance)
        return fragment_name
