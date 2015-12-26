import unittest

import glypy
from glypy.algorithms import similarity

from common import load


class SimilarityTests(unittest.TestCase):

    def test_deep_similarity(self):
        branchy = load("branchy_glycan")
        broad = load("broad_n_glycan")
        ref = broad.clone()
        self.assertEqual(similarity.monosaccharide_similarity(branchy.root, branchy.root), (5, 5))
        self.assertEqual(
            similarity.monosaccharide_similarity(branchy.root, branchy.root, include_children=True),
            (26, 26))
        self.assertEqual(similarity.monosaccharide_similarity(branchy.root, broad.root), (4, 5))
        self.assertEqual(
            similarity.monosaccharide_similarity(branchy.root, broad.root, include_children=True),
            (7, 10))
        self.assertEqual(
            similarity.monosaccharide_similarity(broad.root, branchy.root, include_children=True),
            (11, 14))
        self.assertEqual(similarity.monosaccharide_similarity(broad.root, broad.root, include_children=True), (54, 54))
        self.assertEqual(ref, broad)

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

    def test_has_monosaccharide(self):
        broad = load("broad_n_glycan")
        self.assertTrue(similarity.has_fucose(broad))
        self.assertTrue(similarity.has_modification(broad[6], 'd'))
        self.assertFalse(similarity.has_modification(broad[5], 'd'))


if __name__ == '__main__':
    unittest.main()
