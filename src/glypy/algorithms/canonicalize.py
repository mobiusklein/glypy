from functools import cmp_to_key
from collections import Counter, deque


def all_node_depth(node, visited=None):
    if visited is None:
        visited = set()
    if node.id in visited:  # pragma: no cover
        return 0
    visited.add(node.id)
    depth_count = 1
    children = list(node.children())
    try:
        children += list(node.substituents())
    except AttributeError:
        pass
    if children:
        depth_count += max(all_node_depth(ch, visited) for p, ch in children)
    return depth_count


class CanonicalizerBase(object):
    def __init__(self, structure, reverse=False):
        self.structure = structure
        self.reverse = reverse
        if not self.structure.has_index():
            self.structure.reindex()

    def sort_links(self, links, reverse=False):
        raise NotImplementedError()

    def resort(self):
        nodes = deque(self.structure.leaves())
        seen = set()
        while nodes:
            node = nodes.popleft()
            if node.id in seen:
                continue
            seen.add(node.id)
            link_map = node.links.__class__()
            for link in self.sort_links(node.links.values(), reverse=self.reverse):
                if link.is_parent(node):
                    link_map[link.parent_position] = link
                else:
                    link_map[link.child_position] = link
            node.links = link_map
            for pos, parent in node.parents():
                nodes.append(parent)
        for link in self.structure.link_index:
            link.label = None
        self.structure.label_branches()

    @classmethod
    def canonicalize(cls, structure, **kwargs):
        sorter = cls(structure, **kwargs)
        sorter.resort()
        return structure


class GlycoCTCanonicalizer(CanonicalizerBase):
    def __init__(self, structure, reverse=False):
        super(GlycoCTCanonicalizer, self).__init__(structure, reverse=reverse)
        self.branch_to_terminal_count = self.build_branch_to_terminal_count()

    def get_branch_from_link_label(self, link):
        return link.label[0]

    def build_branch_to_terminal_count(self):
        counter = Counter()
        try:
            for key in sorted(self.structure.branch_parent_map.keys(), reverse=True):
                parent = self.structure.branch_parent_map[key]
                counter[parent] += counter[key] + 1
        except AttributeError:
            pass
        return counter

    def _compare_residue_ordering(self, res_a, res_b):
        n_child_residues_a = all_node_depth(res_a)
        n_child_residues_b = all_node_depth(res_b)
        diff_child_res = n_child_residues_a - n_child_residues_b

        if diff_child_res != 0:
            if diff_child_res < 0:
                return -1
            else:
                return 1

        try:
            branch_length_a = max((all_node_depth(cr) for p, cr in res_a.children()))
        except ValueError:
            branch_length_a = 0
        try:
            branch_length_b = max((all_node_depth(cr) for p, cr in res_b.children()))
        except ValueError:
            branch_length_b = 0

        diff_longest_branch = branch_length_a - branch_length_b

        if diff_longest_branch != 0:
            if diff_longest_branch < 0:
                return -1
            else:
                return 1

        n_branches_from_a = 0
        n_branches_from_b = 0
        for link in res_a.links.values():
            if link.is_parent(res_a):
                branch_label = self.get_branch_from_link_label(link)
                n_branches_from_a = max(n_branches_from_a, self.branch_to_terminal_count[branch_label])

        for link in res_b.links.values():
            if link.is_parent(res_b):
                branch_label = self.get_branch_from_link_label(link)
                n_branches_from_b = max(n_branches_from_b, self.branch_to_terminal_count[branch_label])
        diff_n_branches_from = n_branches_from_a - n_branches_from_b

        if diff_n_branches_from != 0:
            if diff_n_branches_from < 0:
                return -1
            else:
                return 1

        if res_a == res_b:
            return 0

        subtree_a = str(self.structure.subtree_from(self.structure, res_a))
        subtree_b = str(self.structure.subtree_from(self.structure, res_b))
        return (subtree_b > subtree_a) - (subtree_b < subtree_a)

    def compare_residue_ordering(self, res_a, res_b):
        ordered = self._compare_residue_ordering(res_a, res_b)
        return ordered

    def _compare_link_ordering(self, link_a, link_b):
        # Ignoring # of links for now since it is difficult
        # to compute
        parent_pos_a = link_a.parent_position
        parent_pos_b = link_b.parent_position
        try:
            diff_parent = parent_pos_a - parent_pos_b
        except TypeError as e:
            print(parent_pos_a, parent_pos_b, link_a, link_b)
            raise e

        if diff_parent != 0:
            if diff_parent < 0:
                return -1
            else:
                return 1

        child_pos_a = link_a.child_position
        child_pos_b = link_b.child_position
        diff_child = child_pos_a - child_pos_b

        if diff_child != 0:
            if diff_child < 0:
                return -1
            else:
                return 1

        sigils_a = link_a._glycoct_sigils()
        sigils_b = link_b._glycoct_sigils()

        if sigils_a[0] != sigils_b[0]:
            diff_sig0 = ord(sigils_a[0]) - ord(sigils_b[0])
            if diff_sig0 < 0:
                return -1
            else:
                return 1

        if sigils_a[1] != sigils_b[1]:
            diff_sig1 = ord(sigils_a[1]) - ord(sigils_b[1])
            if diff_sig1 < 0:
                return -1
            else:
                return 1

        child_a = link_a.child
        child_b = link_b.child
        ordered = self.compare_residue_ordering(child_a, child_b)
        return ordered

    def compare_link_ordering(self, link_a, link_b):
        ordered = self._compare_link_ordering(link_a, link_b)
        return ordered

    def sort_links(self, links, reverse=False):
        return sorted(links, key=cmp_to_key(self.compare_link_ordering),
                      reverse=reverse)

    def sort_residues(self, residues, reverse=False):
        return sorted(residues, key=cmp_to_key(self.compare_residue_ordering),
                      reverse=reverse)


def canonicalize(structure, canonicalizer=GlycoCTCanonicalizer, **kwargs):
    if canonicalizer is None:
        canonicalizer = GlycoCTCanonicalizer
    return canonicalizer.canonicalize(structure, **kwargs)
