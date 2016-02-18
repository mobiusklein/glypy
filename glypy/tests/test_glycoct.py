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

osvaldo_repeat_1 = '''
RES
1r:r1
REP
REP1:5o(4+1)2d=-1--1
RES
2b:b-dglc-HEX-1:5
3b:b-dgal-HEX-1:5
4b:b-dglc-HEX-1:5
5b:a-dglc-HEX-1:5
6b:a-dglc-HEX-1:5
7s:n-acetyl
8b:b-dgal-HEX-1:5
9b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
10s:n-acetyl
LIN
1:2o(4+1)3d
2:3o(3+1)4d
3:3o(4+1)5d
4:5o(6+1)6d
5:6d(2+1)7n
6:6o(4+1)8d
7:8o(3+2)9d
8:9d(5+1)10n
'''

osvaldo_repeat_2 = '''
RES
1r:r1
REP
REP1:6o(3+1)2n=-1--1
RES
2r:r2
3r:r3
4r:r4
5r:r5
6b:b-dgal-HEX-1:5
LIN
1:2n(4+1)3n
2:3n(3+1)4n
3:4n(3+1)5n
4:5n(3+1)6d
REP2:7o(4+1)7d=-1--1
RES
7b:b-dxyl-PEN-1:5
REP3:8o(3+1)8d=-1--1
RES
8b:b-dgal-HEX-1:5
9b:b-dgal-HEX-1:5
LIN
5:8o(6+1)9d
REP4:10o(3+1)10d=-1--1
RES
10b:b-dgal-HEX-1:5
11b:b-dgal-HEX-1:5
12b:b-dgal-HEX-1:5
13b:b-dgal-HEX-1:5
LIN
6:10o(6+1)11d
7:11o(3+1)12d
8:12o(3+1)13d
REP5:14o(3+1)14d=-1--1
RES
14b:b-dgal-HEX-1:5
15b:b-dgal-HEX-1:5
16b:a-lara-PEN-1:4
17b:a-lara-PEN-1:4
LIN
9:14o(6+1)15d
10:15o(6+1)16d
11:16o(3+1)17d
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

        structure = glycoct.loads(osvaldo_repeat_1)
        self.assertAlmostEqual(structure.mass(), 1322.449, 3)
        self.assertEqual(len(structure), 7)

        structure = glycoct.loads(osvaldo_repeat_2)
        self.assertAlmostEqual(structure.mass(), 1872.612751, 3)
        self.assertEqual(len(structure), 12)


if __name__ == '__main__':
    unittest.main()
