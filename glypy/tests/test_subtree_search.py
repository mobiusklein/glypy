import unittest

from glypy.io import glycoct
from glypy.structure.named_structures import glycans, monosaccharides
from glypy.algorithms import subtree_search

from .common import load


bokan_bao_gly_a = """RES
1b:b-dman-HEX-1:5
2b:a-dman-HEX-1:5
3b:a-dman-HEX-1:5
4b:b-dglc-HEX-1:5
5b:b-dglc-HEX-1:5
6s:n-acetyl
7b:b-dgal-HEX-1:5
8s:n-acetyl
9b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
10s:n-acetyl
LIN
1:1o(3+1)2d
2:1o(6+1)3d
3:2o(4+1)4d
4:3o(6+1)5d
5:4d(2+1)6n
6:4o(4+1)7d
7:5d(2+1)8n
8:7o(3+2)9d
9:9d(5+1)10n"""

bokan_bao_gly_b = """RES
1b:b-dman-HEX-1:5
2b:a-dman-HEX-1:5
3b:a-dman-HEX-1:5
4b:b-dglc-HEX-1:5
5b:b-dglc-HEX-1:5
6s:n-acetyl
7s:n-acetyl
8b:b-dgal-HEX-1:5
9b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
10s:n-acetyl
LIN
1:1o(3+1)2d
2:1o(6+1)3d
3:2o(4+1)4d
4:3o(6+1)5d
5:4d(2+1)6n
6:5d(2+1)7n
7:5o(4+1)8d
8:8o(3+2)9d
9:9d(5+1)10n"""


class SubtreeSearchTests(unittest.TestCase):

    def test_topological_inclusion(self):
        branchy = load("branchy_glycan")
        broad = load("broad_n_glycan")
        match = subtree_search.topological_inclusion(branchy.root, branchy.root)
        assert match == 10
        match = subtree_search.topological_inclusion(branchy.root, broad.get(3))
        assert match
        match = subtree_search.topological_inclusion(branchy.root, broad.root)
        assert not match

        gly_a = glycoct.loads(bokan_bao_gly_a)
        gly_b = glycoct.loads(bokan_bao_gly_b)
        assert subtree_search.topological_inclusion(gly_a.root, gly_b.root)
        assert subtree_search.topological_inclusion(gly_b.root, gly_a.root)

    def test_subtree_of(self):
        branchy = load("branchy_glycan")
        broad = load("broad_n_glycan")
        match = subtree_search.subtree_of(branchy, broad)
        assert match == 3
        core = glycans['N-Linked Core']
        self.assertTrue(subtree_search.subtree_of(core, broad))
        self.assertTrue(subtree_search.subtree_of(broad, core) is None)

    def test_maximum_common_subtree(self):
        core = glycans['N-Linked Core']
        branchy = load("branchy_glycan")
        res = subtree_search.maximum_common_subgraph(core, branchy)
        self.assertEqual(res.score, 6.0)

    def test_find_matching_subtree_roots(self):
        branchy = load("branchy_glycan")
        broad = load("broad_n_glycan")

        roots = subtree_search.find_matching_subtree_roots(branchy, broad)
        assert broad.get(3) in roots

    def test_is_n_glycan(self):
        core = glycans['N-Linked Core']
        tree = load("broad_n_glycan")
        result = (subtree_search.subtree_of(core, tree))
        self.assertTrue(result == 1)
        tree = load("complex_glycan")
        result = (subtree_search.subtree_of(core, tree, exact=False))
        self.assertTrue(result == 1)
        result = (subtree_search.subtree_of(core, tree, exact=True))
        self.assertTrue(result == 1)
        tree = load("branchy_glycan")
        result = (subtree_search.subtree_of(core, tree, exact=False))
        self.assertTrue(result is None)

    def test_disaccharide_similarity(self):
        core = glycans['N-Linked Core']
        self.assertEqual(subtree_search.n_saccharide_similarity(core, core), 1.0)
        copy = core.clone()
        copy.root.add_monosaccharide(monosaccharides.Fucose, 3)
        self.assertEqual(subtree_search.n_saccharide_similarity(core, copy), 0.7)

    def test_treelet_iterator(self):
        complex_glycan = load("complex_glycan")
        treelets = list(subtree_search.treelets(complex_glycan, 3))
        assert len(treelets) == 27
