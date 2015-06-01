import unittest

from glypy.io import glycoct
from common import load, glycan


class GlycoCTParserTests(unittest.TestCase):
    _file_path = "./test_data/glycoct.txt"

    def test_parse_file(self):
        for g in glycoct.read(self._file_path):
            self.assertTrue(isinstance(g, glycan.Glycan))

    def test_parse_cyclical(self):
        structure = load("cyclical_glycan")
        self.assertAlmostEqual(structure.mass(), 810.2641170925)

    def test_parse_repeating(self):
        structure = load("repeating_glycan")
        self.assertAlmostEqual(structure.mass(), 662.190558, 3)
