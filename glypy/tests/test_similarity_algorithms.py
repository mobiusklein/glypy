import unittest
import operator

import glypy
from glypy.composition import composition_transform
from glypy.algorithms import similarity

from .common import load


class SimilarityTests(unittest.TestCase):

    def test_deep_similarity(self):
        branchy = load("branchy_glycan")
        broad = load("broad_n_glycan")
        ref = broad.clone()
        self.assertEqual(similarity.monosaccharide_similarity(branchy.root, branchy.root), (7, 7))
        self.assertEqual(
            similarity.monosaccharide_similarity(branchy.root, branchy.root, include_children=True),
            (64, 64))
        self.assertEqual(similarity.monosaccharide_similarity(branchy.root, broad.root), (4, 7))
        self.assertEqual(
            similarity.monosaccharide_similarity(branchy.root, broad.root, include_children=True),
            (9, 14))
        self.assertEqual(
            similarity.monosaccharide_similarity(broad.root, branchy.root, include_children=True),
            (17, 20))
        self.assertEqual(similarity.monosaccharide_similarity(broad.root, broad.root, include_children=True), (91, 91))
        self.assertEqual(ref, broad)

    def test_build_index_pairs(self):
        nsc = similarity.NodeSimilarityComparator()
        pairs = {(18, 15): (9, 9), (15, 15): (9, 9), (15, 18): (9, 9), (18, 18): (9, 9)}
        expected = set([((18, 15), (15, 18)), ((18, 18), (15, 15))])
        result = nsc.build_unique_index_pairs(pairs)
        self.assertEqual(set(result), expected)

    def test_optimal_assignment(self):
        nsc = similarity.NodeSimilarityComparator()
        pairs = {(18, 15): (9, 9), (15, 15): (9, 9), (15, 18): (9, 9), (18, 18): (9, 9)}
        expected = set(((18, 15), (15, 18)))
        result = nsc.optimal_assignment(pairs)
        self.assertEqual(set(result), expected)

    def test_partial_similarity(self):
        broad = load("broad_n_glycan")
        expected = [
            (1, (91, 91)),
            (3, (84, 84)),
            (5, (77, 77)),
            (6, (39, 39)),
            (7, (20, 20)),
            (9, (6, 6)),
            (10, (7, 7)),
            (11, (13, 13)),
            (13, (6, 6)),
            (14, (32, 32)),
            (15, (13, 13)),
            (17, (6, 6)),
            (18, (13, 13)),
            (20, (6, 6))
        ]
        result = list(map(lambda x: (
            x[0].id, similarity.monosaccharide_similarity(
                x[0], x[1], include_children=1)), zip(broad, broad)))
        self.assertEqual(result, expected)

    def test_is_aminated(self):
        broad = load("broad_n_glycan")
        self.assertTrue(similarity.is_aminated(broad.root))
        self.assertTrue(similarity.has_substituent(broad.root, glypy.Substituent("n-acetyl")))
        self.assertFalse(similarity.is_aminated(broad[3]))

    def test_is_reduced(self):
        broad = load("broad_n_glycan")
        broad.reducing_end = True
        self.assertTrue(similarity.is_reduced(broad))
        self.assertFalse(similarity.is_reduced(None))

    def test_ignore_reduction(self):
        broad = load("broad_n_glycan")
        broad.reducing_end = True
        self.assertTrue(similarity.is_reduced(broad))
        r, q = similarity.monosaccharide_similarity(
            broad.root, broad.root)
        r2, q2 = similarity.monosaccharide_similarity(
            broad.root, broad.root, ignore_reduction=True)
        self.assertTrue((r - 1) == r2)

    def test_has_monosaccharide(self):
        broad = load("broad_n_glycan")
        self.assertTrue(similarity.has_fucose(broad))
        self.assertTrue(similarity.has_modification(broad[6], 'd'))
        self.assertFalse(similarity.has_modification(broad[5], 'd'))

    def test_similar_substituents(self):
        self.assertTrue(
            similarity.monosaccharide_similarity(
                glypy.Substituent('n-acetyl'),
                glypy.Substituent('n-acetyl')) == (1, 1))

    def test_is_derivatized(self):
        broad = load("broad_n_glycan")
        self.assertFalse(similarity.is_derivatized(broad.root))
        composition_transform.derivatize(broad, 'methyl')
        self.assertTrue(similarity.is_derivatized(broad.root))


if __name__ == '__main__':
    unittest.main()
