import unittest

import glypy
from glypy.io import glycoct
from glypy.tests.common import load, glycan, structures as raw_structures


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


test_multiple_repeating_buffer = '''
RES
1b:b-dglc-HEX-1:5
2s:n-acetyl
3b:b-dglc-HEX-1:5
4s:n-acetyl
5b:b-dman-HEX-1:5
6b:a-dman-HEX-1:5
7b:b-dglc-HEX-1:5
8s:n-acetyl
9r:r1
10b:b-dgal-HEX-1:5
11b:a-dman-HEX-1:5
12b:b-dglc-HEX-1:5
13s:n-acetyl
14r:r2
15b:b-dgal-HEX-1:5
16b:a-lgal-HEX-1:5|6:d
LIN
1:1d(2+1)2n
2:1o(4+1)3d
3:3d(2+1)4n
4:3o(4+1)5d
5:5o(3+1)6d
6:6o(2+1)7d
7:7d(2+1)8n
8:7o(4+1)9n
9:9n(4+1)10d
10:5o(6+1)11d
11:11o(2+1)12d
12:12d(2+1)13n
13:12o(4+1)14n
14:14n(4+1)15d
15:1o(6+1)16d
REP
REP1:18o(4+1)17d=-1--1
RES
17b:b-dgal-HEX-1:5
18b:b-dglc-HEX-1:5
19s:n-acetyl
LIN
16:17o(3+1)18d
17:18d(2+1)19n
REP2:21o(4+1)20d=-1--1
RES
20b:b-dgal-HEX-1:5
21b:b-dglc-HEX-1:5
22s:n-acetyl
LIN
18:20o(3+1)21d
19:21d(2+1)22n


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
10b:b-dglc-HEX-1:5
11s:n-acetyl
12b:b-dgal-HEX-1:5
13b:a-dman-HEX-1:5
14b:b-dglc-HEX-1:5
15s:n-acetyl
16b:b-dgal-HEX-1:5
17b:b-dglc-HEX-1:5
18s:n-acetyl
19b:b-dgal-HEX-1:5
20b:a-lgal-HEX-1:5|6:d
LIN
1:1d(2+1)2n
2:1o(4+1)3d
3:3d(2+1)4n
4:3o(4+1)5d
5:5o(3+1)6d
6:6o(2+1)7d
7:7d(2+1)8n
8:7o(4+1)9d
9:6o(4+1)10d
10:10d(2+1)11n
11:10o(4+1)12d
12:5o(6+1)13d
13:13o(2+1)14d
14:14d(2+1)15n
15:14o(4+1)16d
16:13o(6+1)17d
17:17d(2+1)18n
18:17o(4+1)19d
19:1o(6+1)20d
'''


und_case1 = '''
RES
1b:a-dman-HEX-1:5
2r:r1
3b:a-dgal-HEX-1:5
LIN
1:1o(2+1)2n
2:2n(2+1)3d
REP
REP1:4o(2+1)4d=-1--1
RES
4b:a-dglc-HEX-1:5
UND
UND1:100.0:100.0
ParentIDs:3|1
SubtreeLinkageID1:o(2+4)d
RES
5b:b-dgro-dgal-NON-0:0|1:a|2:keto|3:d
6s:n-acetyl
LIN
3:5d(5+2)6n
'''


und_case2 = '''
RES
1b:a-dman-HEX-1:5
2r:r1
3b:a-dgal-HEX-1:5
LIN
1:1o(2+1)2n
2:2n(2+1)3d
REP
REP1:4o(2+1)4d=1-5
RES
4b:a-dglc-HEX1:5
UND
UND1:60.0:80.0
ParentIDs:4
SubtreeLinkageID1:o(4+1)n
RES
5s:sulfate

'''


class GlycoCTParserTests(unittest.TestCase):
    _file_path = "./test_data/glycoct.txt"

    def test_parse_file(self):
        with open(self._file_path) as stream:
            for g in glycoct.read(stream):
                self.assertTrue(isinstance(g, glycan.Glycan))

    def test_parse_cyclical(self):
        structure = load("cyclical_glycan")
        self.assertAlmostEqual(structure.mass(), 810.2641170925)

    def test_parse_multiple_structures_in_buffer(self):
        structures = glycoct.loads(test_multiple_repeating_buffer)
        self.assertEqual(len(structures), 2)

        self.assertAlmostEqual(structures[0].mass(), structures[1].mass(), 3)
        self.assertNotEqual(*structures)

    def test_parse_repeating(self):
        structure = load("repeating_glycan")
        self.assertAlmostEqual(structure.mass(), 666.22185, 3)

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

    def test_parse_undetermined(self):
        structure = glycoct.loads(und_case1)
        for link in structure.link_index:
            if link.is_ambiguous():
                self.assertTrue(len(link.child.stem) > 1)

        structure = glycoct.loads(und_case2)
        sulfate_count = 0
        for node in structure:
            for i, subst in node.substituents():
                if subst.name == 'sulfate':
                    sulfate_count += 1
        self.assertEqual(sulfate_count, 5)

    def test_composition(self):
        gc = glypy.GlycanComposition.parse("{Hex:5; Hex2NAc:4; Neu5Ac:2; @sulfate: 1}")
        self.assertEqual(gc, glycoct.loads(glycoct.dumps(gc)))

    def test_ordered_serialization(self):
        for acc in ['G58143RL', 'G82388RB', 'G28839WC', 'G37369XO', 'G36221RT', 'G27293OK',
                    'G62831KM', 'G65832JS', 'G07337US', 'G01120GS']:
            text = raw_structures[acc].strip()
            structure = glycoct.loads(text)
            retext = glycoct.dumps(structure).strip()
            assert text == retext, "Failed to match %s" % acc


if __name__ == '__main__':
    unittest.main()
