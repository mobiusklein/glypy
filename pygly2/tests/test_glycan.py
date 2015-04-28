import unittest
from common import load, glycoct, glycan, multimap, pickle, named_structures, monosaccharides

Glycan = glycan.Glycan


class GlycanTests(unittest.TestCase):
    _file_path = "./test_data/glycoct.txt"

    def test_from_glycoct(self):
        for structure in glycoct.read(self._file_path):
            self.assertEqual(
                structure, glycoct.loads(structure.to_glycoct()).next())

    def test_fragments_preserve(self):
        structure = load("branchy_glycan")
        dup = structure.clone()
        self.assertEqual(structure, dup)
        list(dup.fragments('AXZCBY', 2))

        self.assertEqual(structure, dup)

    def test_branch_counts(self):
        structure = load("branchy_glycan")
        self.assertEqual(structure.count_branches(), 3)

    def test_fragments_mass(self):
        structure = load("common_glycan")
        frags = list(structure.fragments('ZCBY', 1))
        import json
        frags_ref = json.load(open('test_data/fragments-example.json'))
        container = multimap.MultiMap()
        for frag in frags_ref:
            container[frag['kind']] = frag['mass']

        def almost_equal(a, b):
            e = 0.0001
            return (a - e) <= b <= (a + e)

        for frag in frags:
            structure.name_fragment(frag)
            kind = frag.kind
            mass = frag.mass
            candidates = container[kind] + container[kind[::-1]]
            if len(candidates) == 0:
                raise AssertionError("No candidates found for {}".format(frag))
            res = (any(almost_equal(mass, x) for x in candidates))
            if not res:
                raise AssertionError(
                    "{} found no matches in {}".format(frag, candidates))

    def test_disjoint_subtrees(self):
        structure = load("common_glycan")
        f = open('test_data/test_disjoint_subtrees.pkl', 'rb')
        refs = {}
        while True:
            try:
                for ref in pickle.load(f):
                    refs[tuple(ref.link_ids)] = ref
            except EOFError:
                break
        for tree in structure.substructures(3):
            ref = refs[tuple(tree.link_ids)]
            link_id_eq = tree.link_ids == ref.link_ids
            parent_tree_eq = tree.parent_tree == ref.parent_tree
            if not parent_tree_eq or not link_id_eq:
                raise AssertionError("DisjointTrees did not match")
                self.assertEqual(list(tree.parent_tree), list(ref.parent_tree))

    def test_reducing_end(self):
        structure = load("common_glycan")
        self.assertEqual(structure.reducing_end, None)
        structure.reducing_end = 1
        self.assertEqual(structure.reducing_end, 1)

    def test_clone(self):
        structure = load("common_glycan")
        ref = structure.clone()
        structure.reducing_end = 1
        self.assertTrue(structure != ref)

    def test_indexing(self):
        structure = load("common_glycan")
        ref = structure.clone()
        for i, node in enumerate(structure.index):
            self.assertEqual(node.id, ref[i].id)
        structure.deindex()
        for i, node in enumerate(structure.index):
            self.assertNotEqual(node.id, ref[i].id)
        structure.index = None
        self.assertRaises(IndexError, lambda: structure[0])

    def test_traversal(self):
        structure = load("common_glycan")
        structure[-
                  1].add_monosaccharide(named_structures.monosaccharides['NeuGc'])
        structure.reindex(method='dfs')
        ref = structure.clone()
        self.assertEqual(structure[-1], ref[-1])
        structure.reindex(method='bfs')
        self.assertNotEqual(structure[-1], ref[-1])

    def test_traversal_by_name(self):
        structure = load("common_glycan")
        structure[-
                  1].add_monosaccharide(named_structures.monosaccharides['NeuGc'])
        structure.reindex(method='dfs')
        ref = structure.clone()
        self.assertEqual(structure, ref)
        structure.reindex(method='depth_first_traversal')
        self.assertEqual(structure, ref)
        self.assertRaises(
            AttributeError, lambda: structure.reindex(method='not_real_traversal'))

    def test_leaves(self):
        structure = load("common_glycan")
        leaves = list(structure.leaves())
        for node in leaves:
            self.assertTrue(len(list(node.children())) == 0)
        leaves = list(structure.leaves(bidirectional=True))
        for node in leaves:
            self.assertTrue(
                len(list(node.children())) == 0 or node == structure.root)

    def test_custom_traversal_method(self):
        def rev_sort_dfs(self, visited=None, from_node=None, *args, **kwargs):
            node_stack = list([self.root])
            visited = set()
            while len(node_stack) > 0:
                node = node_stack.pop()
                if node.id in visited:
                    continue
                visited.add(node.id)
                yield (node)
                node_stack.extend(reversed(list(terminal for pos, link in node.links.items()
                                                for terminal in link if terminal.id not in visited and
                                                len(link.child.substituent_links) < 1)))
        structure = load("common_glycan")
        structure[-3].add_monosaccharide(
            named_structures.monosaccharides['Hex'], 4).add_substituent('methyl', 5)
        structure.reindex()
        ref = structure.clone()
        self.assertEqual(structure[-1], ref[-1])
        structure.reindex(method=rev_sort_dfs)
        self.assertNotEqual(structure[-1], ref[-1])

    def test_topological_equality(self):
        base = load("branchy_glycan")
        a = base.clone()
        b = base.clone()
        c = base.clone()
        self.assertEqual(base, b)
        self.assertEqual(a, b)
        list(a.leaves())[0].add_monosaccharide(monosaccharides["NeuGc"])
        list(b.leaves())[1].add_monosaccharide(monosaccharides["NeuGc"])
        list(c.leaves())[2].add_monosaccharide(monosaccharides["NeuGc"])
        d = base.clone()
        d_children = list(d.leaves())
        d_children[0].add_monosaccharide(monosaccharides["NeuGc"])
        d_children[1].add_monosaccharide(monosaccharides["NeuGc"])
        self.assertTrue(a.topological_equality(b))
        self.assertFalse(a.topological_equality(c))
        self.assertFalse(a.topological_equality(base))
        self.assertFalse(a.topological_equality(d))
        self.assertFalse(b.topological_equality(d))

    def test_substructures_does_not_mutate(self):
        structure = load("broad_n_glycan")
        ref = structure.clone()
        for substructure in structure.substructures(max_cleavages=2):
            pass
        self.assertEqual(structure, ref)

    def test_crossring(self):
        structure = load("branchy_glycan")
        frag_data = {
            '0,2A': 485.174,
            '0,2X': 1317.470,
            '0,3A': 455.163,
            '0,3X': 1347.481,
            '0,4A': 425.153,
            '0,4X': 1377.491,
            '1,3A': 425.153,
            '1,3X': 1377.491,
            '1,4A': 455.163,
            '1,4X': 1347.481,
            '1,5A': 864.322,
            '1,5X': 938.322,
            '2,4X': 1742.623,
            '2,5A': 469.179,
            '2,5X': 1333.465,
            '3,5A': 439.168,
            '3,5X': 1363.4769
        }
        for fragment in structure.fragments("AX"):
            if fragment.link_ids == [3]:
                self.assertAlmostEqual(frag_data[fragment.kind], fragment.mass, 2)

    def test_subtree_from(self):
        structure = load("branchy_glycan")
        child = structure.root.children().next()[1]
        subtree = Glycan.subtree_from(structure, child)
        temp = structure.clone()
        temproot = temp.root.children().next()[1]
        for link in temp.root.links.values():
            link.break_link(refund=True)
        temp.root = temproot
        self.assertEqual(temp, subtree)
        self.assertEqual(Glycan.subtree_from(structure, 1), temp)

    def test_cyclic_clone(self):
        structure = load("cyclical_glycan")
        self.assertEqual(structure, structure.clone())

if __name__ == '__main__':
    unittest.main()
