import unittest
from .common import load, glycoct, glycan, multimap, pickle, named_structures, monosaccharides

from glypy import Substituent, tree
from glypy.structure.fragment import Fragment

Glycan = glycan.Glycan


class GlycanTests(unittest.TestCase):
    _file_path = "./test_data/glycoct.txt"

    def test_from_glycoct(self):
        with open(self._file_path) as stream:
            for structure in glycoct.read(stream):
                self.assertEqual(
                    structure, glycoct.loads(structure.serialize("glycoct")))

    def test_fragments_preserve(self):
        structure = load("branchy_glycan")
        dup = structure.clone()
        self.assertEqual(structure, dup)
        list(dup.fragments('ABY', 2))

        self.assertEqual(structure, dup)

    def test_branch_counts(self):
        structure = load("branchy_glycan")
        self.assertEqual(structure.count_branches(), 3)

    def test_fragments_mass(self):
        structure = load("common_glycan")
        frags = list(structure.fragments('ZCBY', 1))
        import json
        with open('test_data/fragments-example.json') as f:
            frags_ref = json.load(f)
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
        structure[-1].add_monosaccharide(
            named_structures.monosaccharides['NeuGc'])
        structure.reindex(method='dfs')
        ref = structure.clone()
        self.assertEqual(structure, ref)
        structure.reindex(method='depth_first_traversal')
        self.assertEqual(structure, ref)
        self.assertRaises(
            KeyError, lambda: structure.reindex(
                method='not_real_traversal'))

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
        def rev_sort_dfs(self, apply_fn=lambda x: x, visited=None, from_node=None, *args, **kwargs):
            node_stack = list([self.root])
            visited = set()
            while len(node_stack) > 0:
                node = node_stack.pop()
                if node.id in visited:
                    continue
                visited.add(node.id)
                yield apply_fn(node)
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
        for fragment in structure.fragments("AXB", max_cleavages=2):
            if fragment.link_ids == [3]:
                self.assertAlmostEqual(frag_data[fragment.kind], fragment.mass, 2)

    def test_subtree_from(self):
        structure = load("branchy_glycan")
        child = structure.root.children()[0][1]
        subtree = Glycan.subtree_from(structure, child)
        temp = structure.clone()
        temproot = temp.root.children()[0][1]
        for link in temp.root.links.values():
            link.break_link(refund=True)
        temp.root = temproot
        self.assertEqual(temp, subtree)
        self.assertEqual(Glycan.subtree_from(structure, 1), temp)

    def test_cyclic_clone(self):
        import warnings
        warnings.simplefilter("error")
        with self.assertRaises(UserWarning):
            structure = load("cyclical_glycan")
            self.assertEqual(structure, structure.clone())

    def test_subtree_from_fragment(self):
        structure = load("branchy_glycan")
        fragment = Fragment(mass=1463.5284494, kind="1,5A0,3X", included_nodes=set([2, 3, 4, 5, 6, 7, 8, 9, 10]),
                            link_ids={}, name="0,3Xc4-1,5A4",
                            crossring_cleavages={3: ('1,5', 'A'), 10: ('0,3', 'X')})
        subtree = glycan.fragment_to_substructure(fragment, structure)
        self.assertAlmostEqual(subtree.mass(), fragment.mass)

    def test_get(self):
        structure = load("branchy_glycan")
        self.assertEqual(structure.get(1), structure.root)
        self.assertEqual(structure.get_link(1), next(structure.iterlinks())[1])

    def test_iternodes(self):
        structure = load("branchy_glycan")
        for a, b in zip(iter(structure), structure.iternodes()):
            self.assertEqual(a, b)

    def test_iterlinks_substituents(self):
        had_substituents = False
        structure = load("branchy_glycan")
        for p, link in structure.iterlinks(substituents=True):
            if isinstance(link.child, Substituent):
                had_substituents = True
        self.assertTrue(had_substituents)

    def test_subtrees(self):
        structure = load("branchy_glycan")
        g = structure.subtrees()
        s = next(g)
        self.assertAlmostEqual(tree(s).mass(), 221.0899, 2)
        s = next(g)
        self.assertAlmostEqual(tree(s).mass(), 1599.565, 2)

    def test_fragment_properties(self):
        structure = load("branchy_glycan")
        frags = {f.name: f for f in structure.fragments("Y")}
        self.assertEqual(frags["Ya2"].fname, r"Y$_\alpha$2")
        Ya2 = frags["Ya2"]
        self.assertTrue(Ya2.is_reducing())
        self.assertFalse(Ya2.is_non_reducing())
        self.assertFalse(Ya2.is_internal())

    def test_iterconfigurations(self):
        G81339YK_iupac = (
            '?-D-GalpNAc-(1->3/4)[?-L-Fucp-(1->3/4)]-beta-D-GlcpNAc-(1->3)[?-D-GalpNAc-(1->3/4)[?-L-Fucp-(1->3/4)]'
            '-beta-D-GlcpNAc-(1->6)]-?-D-GalpNAc(1->')
        G81339YK_glycoct = '''RES
1b:x-dgal-HEX-1:5
2s:n-acetyl
3b:b-dglc-HEX-1:5
4b:x-dgal-HEX-1:5
5s:n-acetyl
6b:x-lgal-HEX-1:5|6:d
7s:n-acetyl
8b:b-dglc-HEX-1:5
9b:x-dgal-HEX-1:5
10s:n-acetyl
11b:x-lgal-HEX-1:5|6:d
12s:n-acetyl
LIN
1:1d(2+1)2n
2:1o(3+1)3d
3:3o(3|4+1)4d
4:4d(2+1)5n
5:3o(3|4+1)6d
6:3d(2+1)7n
7:1o(6+1)8d
8:8o(3|4+1)9d
9:9d(2+1)10n
10:8o(3|4+1)11d
11:8d(2+1)12n
'''
        from glypy.io import iupac, glycoct
        structure_iupac = iupac.loads(G81339YK_iupac)
        structure_glycoct = glycoct.loads(G81339YK_glycoct)

        structure = structure_iupac
        configurations = []
        for config_list in structure.iterconfiguration():
            instance = structure.clone()
            for link, conf in config_list:
                link = instance.get_link(link.id)
                parent = instance.get(conf[0].id)
                child = instance.get(conf[1].id)
                link.reconfigure(parent, child, conf[2], conf[3])
            configurations.append(instance)
        configurations_iupac = configurations
        assert len(configurations) == 4

        structure = structure_glycoct
        configurations = []
        for config_list in structure.iterconfiguration():
            instance = structure.clone()
            for link, conf in config_list:
                link = instance.get_link(link.id)
                parent = instance.get(conf[0].id)
                child = instance.get(conf[1].id)
                link.reconfigure(parent, child, conf[2], conf[3])
            configurations.append(instance)
        configurations_glycoct = configurations
        assert len(configurations) == 4

        pairs = zip(sorted(configurations_glycoct, key=str), sorted(configurations_iupac, key=str))
        for a, b in pairs:
            assert a.canonicalize() == b.canonicalize()

    def test_pickles(self):
        ambig = '''
        RES
        1b:x-dglc-HEX-x:x
        2s:n-acetyl
        3b:b-dglc-HEX-1:5
        4s:n-acetyl
        5b:b-dman-HEX-1:5
        6b:a-dman-HEX-1:5
        7b:b-dglc-HEX-1:5
        8s:n-acetyl
        9b:b-dgal-HEX-1:5
        10b:a-dman-HEX-1:5
        11b:a-dman-HEX-1:5
        LIN
        1:1d(2+1)2n
        2:1o(4+1)3d
        3:3d(2+1)4n
        4:3o(4+1)5d
        5:5o(3|6+1)6d
        6:6o(2+1)7d
        7:7d(2+1)8n
        8:7o(4+1)9d
        9:5o(3|6+1)10d
        10:10o(3|6+1)11d
        '''

        x = glycoct.loads(ambig)

        result = pickle.loads(pickle.dumps(x, -1))
        assert x == result

        ambig = """
        RES
        1b:x-dglc-HEX-1:5
        2b:x-dglc-HEX-1:5
        3b:x-dman-HEX-1:5
        4b:a-dman-HEX-1:5
        5b:b-dglc-HEX-1:5
        6s:n-acetyl
        7b:a-dman-HEX-1:5
        8b:b-dglc-HEX-1:5
        9s:n-acetyl
        10b:b-dglc-HEX-1:5
        11b:x-dgal-HEX-1:5
        12b:x-dgro-dgal-NON-2:6|1:a|2:keto|3:d
        13s:n-acetyl
        14s:n-acetyl
        15s:n-acetyl
        16s:n-acetyl
        17b:a-lgal-HEX-1:5|6:d
        LIN
        1:1o(-1+1)2d
        2:2o(-1+1)3d
        3:3o(3+1)4d
        4:4o(2|4+1)5d
        5:5d(2+1)6n
        6:3o(6+1)7d
        7:7o(2+1)8d
        8:8d(2+1)9n
        9:7o(6+1)10d
        10:10o(-1+1)11d
        11:11o(-1+2)12d
        12:12d(5+1)13n
        13:10d(2+1)14n
        14:2d(2+1)15n
        15:1d(2+1)16n
        16:1o(6+1)17d
        """
        x = glycoct.loads(ambig)
        result = pickle.loads(pickle.dumps(x, -1))
        assert x == result

        ambig = """
        RES
        1b:b-dglc-HEX-1:5
        2s:n-acetyl
        3b:b-dglc-HEX-1:5
        4s:n-acetyl
        5b:b-dman-HEX-1:5
        6b:a-dman-HEX-1:5
        7b:b-dglc-HEX-1:5
        8s:n-acetyl
        9b:b-dgal-HEX-1:5
        10b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
        11s:n-acetyl
        12b:b-dglc-HEX-1:5
        13s:n-acetyl
        14b:b-dgal-HEX-1:5
        15b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
        16s:n-acetyl
        17b:a-dman-HEX-1:5
        18b:b-dglc-HEX-1:5
        19s:n-acetyl
        20b:b-dgal-HEX-1:5
        21b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
        22s:n-acetyl
        23b:a-lgal-HEX-1:5|6:d
        LIN
        1:1d(2+1)2n
        2:1o(4+1)3d
        3:3d(2+1)4n
        4:3o(4+1)5d
        5:5o(3|6+1)6d
        6:6o(2+1)7d
        7:7d(2+1)8n
        8:7o(4+1)9d
        9:9o(3|6+2)10d
        10:10d(5+1)11n
        11:6o(4|6+1)12d
        12:12d(2+1)13n
        13:12o(4+1)14d
        14:14o(3|6+2)15d
        15:15d(5+1)16n
        16:5o(3|6+1)17d
        17:17o(2+1)18d
        18:18d(2+1)19n
        19:18o(4+1)20d
        20:20o(3|6+2)21d
        21:21d(5+1)22n
        22:1o(6+1)23d
        """

        x = glycoct.loads(ambig)
        result = pickle.loads(pickle.dumps(x, -1))
        assert x == result

        ambig = """
        RES
        1b:x-dglc-HEX-1:5
        2s:n-acetyl
        3b:b-dglc-HEX-1:5
        4s:n-acetyl
        5b:b-dman-HEX-1:5
        6b:a-dman-HEX-1:5
        7b:b-dglc-HEX-1:5
        8s:n-acetyl
        9b:b-dgal-HEX-1:5
        10b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
        11s:n-acetyl
        12b:b-dglc-HEX-1:5
        13s:n-acetyl
        14b:b-dgal-HEX-1:5
        15b:a-dman-HEX-1:5
        16b:b-dglc-HEX-1:5
        17s:n-acetyl
        18b:b-dgal-HEX-1:5
        19b:b-dglc-HEX-1:5
        20s:n-acetyl
        21b:b-dgal-HEX-1:5
        22b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
        23s:n-acetyl
        24b:b-dglc-HEX-1:5
        25s:n-acetyl
        26b:b-dgal-HEX-1:5
        27b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
        28s:n-acetyl
        29b:a-lgal-HEX-1:5|6:d
        LIN
        1:1d(2+1)2n
        2:1o(4+1)3d
        3:3d(2+1)4n
        4:3o(4+1)5d
        5:5o(3+1)6d
        6:6o(2+1)7d
        7:7d(2+1)8n
        8:7o(4+1)9d
        9:9o(3|6+2)10d
        10:10d(5+1)11n
        11:6o(4+1)12d
        12:12d(2+1)13n
        13:12o(4+1)14d
        14:5o(6+1)15d
        15:15o(2+1)16d
        16:16d(2+1)17n
        17:16o(4+1)18d
        18:18o(3+1)19d
        19:19d(2+1)20n
        20:19o(4+1)21d
        21:21o(3|6+2)22d
        22:22d(5+1)23n
        23:15o(6+1)24d
        24:24d(2+1)25n
        25:24o(4+1)26d
        26:26o(3|6+2)27d
        27:27d(5+1)28n
        28:1o(6+1)29d
        """
        x = glycoct.loads(ambig)
        result = pickle.loads(pickle.dumps(x, -1))
        assert x == result


if __name__ == '__main__':
    unittest.main()
