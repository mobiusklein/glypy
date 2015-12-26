import unittest

from glypy.io import glycoct
from common import load, glycan


example_multiplicity = '''
RES
1r:r1
REP
REP1:3o(3+1)2d=1-4
RES
2b:a-dman-HEX-1:5
3b:a-lman-HEX-1:5|6:d
4b:b-lxyl-PEN-1:5
LIN
1:2o(4+1)3d
2:3o(2+1)4d
'''


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

        structure = load("big_repeating_glycan")
        self.assertAlmostEqual(structure.mass(), 5248.87193, 3)

        structure = glycoct.loads(example_multiplicity)
        self.assertEqual(len(structure), 4 * 3)
