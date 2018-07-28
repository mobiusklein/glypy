import unittest

from glypy import GlycanComposition
from glypy.io import glycoct, wurcs


G71237SD_wurcs = (
    'WURCS=2.0/6,18,17/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5]'
    '[Aad21122h-2a_2-6_5*NCC/3=O][a1221m-1a_1-5]/1-1-2-3-1-4-5-1-4-5-3-1-4-5-1-4-5-6/a4-b1_a6'
    '-r1_b4-c1_d2-e1_d4-h1_e4-f1_h4-i1_k2-l1_k6-o1_l4-m1_o4-p1_c?-d1_c?-k1_f?-g2_i?-j2_m?-n2_p?-q2')

G71237SD_glycoct = '''RES
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
23b:b-dglc-HEX-1:5
24s:n-acetyl
25b:b-dgal-HEX-1:5
26b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
27s:n-acetyl
28b:a-lgal-HEX-1:5|6:d
LIN
1:1d(2+1)2n
2:1o(4+1)3d
3:3d(2+1)4n
4:3o(4+1)5d
5:5o(-1+1)6d
6:6o(2+1)7d
7:7d(2+1)8n
8:7o(4+1)9d
9:9o(-1+2)10d
10:10d(5+1)11n
11:6o(6+1)12d
12:12d(2+1)13n
13:12o(4+1)14d
14:14o(-1+2)15d
15:15d(5+1)16n
16:5o(-1+1)17d
17:17o(2+1)18d
18:18d(2+1)19n
19:18o(4+1)20d
20:20o(-1+2)21d
21:21d(5+1)22n
22:17o(4+1)23d
23:23d(2+1)24n
24:23o(4+1)25d
25:25o(-1+2)26d
26:26d(5+1)27n
27:1o(6+1)28d'''


G35323LT_wurcs = ('WURCS=2.0/8,12,11/[AUd21122h][a11221h-1a_1-5_7*OPO/3O/3=O][a11221h-1a_1-5][a11222h-1a_1-5]'
                  '[a2122h-1b_1-5_2*NCC/3=O][a2112h-1b_1-5][a1221m-1a_1-5][a2122h-1a_1-5]/1-2-3-4-4-5-6-7-7-6-'
                  '8-8/a5-b1_b3-c1_c2-d1_d2-e1_d7-j1_e7-f1_f3-g1_f4-i1_g2-h1_j4-k1_k3-l1')

G35323LT_glycoct = '''RES
1b:x-dgro-dgal-NON-x:x|1:a|2:keto|3:d
2b:a-lgro-dman-HEP-1:5
3b:a-lgro-dman-HEP-1:5
4b:a-dgro-dman-HEP-1:5
5b:a-dgro-dman-HEP-1:5
6b:b-dglc-HEX-1:5
7s:n-acetyl
8b:b-dgal-HEX-1:5
9b:a-lgal-HEX-1:5|6:d
10b:a-lgal-HEX-1:5|6:d
11b:b-dgal-HEX-1:5
12b:a-dglc-HEX-1:5
13b:a-dglc-HEX-1:5
14s:phosphate
LIN
1:1o(5+1)2d
2:2o(3+1)3d
3:3o(2+1)4d
4:4o(2+1)5d
5:5o(7+1)6d
6:6d(2+1)7n
7:6o(3+1)8d
8:8o(2+1)9d
9:6o(4+1)10d
10:4o(7+1)11d
11:11o(4+1)12d
12:12o(3+1)13d
13:2o(7+1)14n'''


G41928NU_wurcs = ('WURCS=2.0/5,12,11/[h2112h_2*NCC/3=O][a2112h-1b_1-5][Aad21122h-2a_2-6_5*NCC/3=O]'
                  '[a2122h-1b_1-5_2*NCC/3=O][a1221m-1a_1-5]/1-2-3-4-5-2-4-2-5-5-4-2/a3-b1_a6-d1_b3'
                  '-c2_d3-e1_d4-f1_f3-g1_f6-k1_g3-h1_g4-j1_h2-i1_k4-l1')

G41928NU_glycoct = '''RES
1b:o-dgal-HEX-0:0|1:aldi
2s:n-acetyl
3b:b-dgal-HEX-1:5
4b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
5s:n-acetyl
6b:b-dglc-HEX-1:5
7s:n-acetyl
8b:a-lgal-HEX-1:5|6:d
9b:b-dgal-HEX-1:5
10b:b-dglc-HEX-1:5
11s:n-acetyl
12b:b-dgal-HEX-1:5
13b:a-lgal-HEX-1:5|6:d
14b:a-lgal-HEX-1:5|6:d
15b:b-dglc-HEX-1:5
16s:n-acetyl
17b:b-dgal-HEX-1:5
LIN
1:1d(2+1)2n
2:1o(3+1)3d
3:3o(3+2)4d
4:4d(5+1)5n
5:1o(6+1)6d
6:6d(2+1)7n
7:6o(3+1)8d
8:6o(4+1)9d
9:9o(3+1)10d
10:10d(2+1)11n
11:10o(3+1)12d
12:12o(2+1)13d
13:10o(4+1)14d
14:9o(6+1)15d
15:15d(2+1)16n
16:15o(4+1)17d'''


class TestWURCS(unittest.TestCase):
    def test_parse_wurcs(self):
        ref = glycoct.loads(G71237SD_glycoct)
        test = wurcs.loads(G71237SD_wurcs)
        self.assertEqual(ref, test)

        ref = glycoct.loads(G35323LT_glycoct)
        test = wurcs.loads(G35323LT_wurcs)
        self.assertEqual(ref, test)

        ref = glycoct.loads(G41928NU_glycoct)
        test = wurcs.loads(G41928NU_wurcs)
        self.assertEqual(ref, test)

    def test_write_wurcs(self):
        ref = glycoct.loads(G71237SD_glycoct)
        test = wurcs.loads(wurcs.dumps(ref))
        self.assertEqual(ref, test)

        ref = glycoct.loads(G41928NU_glycoct)
        test = wurcs.loads(wurcs.dumps(ref))
        self.assertEqual(ref, test)

    def test_serialize_partially_undefined_monosaccharides(self):
        monosaccharide_symbols = [
            "[uxxxxh]",
            "[a2112h-1x_1-5]",
            "[a2122h-1b_1-5_2*NCC/3=O]",
            "[a1221m-1a_1-5]",
            "[a2112h-1b_1-5]"
        ]
        G11275IL_wurcs = ("WURCS=2.0/5,6,5/[uxxxxh][a2112h-1x_1-5][a2122h-1b_1-5_2*NCC/3=O][a2112h-1b_1-5]"
                          "[a1221m-1a_1-5]/1-2-3-4-5-5/b3-c1_c3-d1_c4-f1_d2-e1_a?-b1")
        test = wurcs.loads(G11275IL_wurcs)
        string = wurcs.dumps(test)
        for symbol in monosaccharide_symbols:
            self.assertIn(symbol, string)

    def test_glycan_composition(self):
        source = '{Fuc:1; Gal:4; Man:3; Glc2NAc:6; Neu5Ac:4}'
        gc = GlycanComposition.parse(source)
        text = wurcs.dumps(gc)
        parser = wurcs.parser.WURCSParser(text)
        test = parser.parse()
        counts = parser.parse_counts()
        self.assertEqual(counts[-1], 0)
        self.assertEqual(gc, test)
        self.assertAlmostEqual(gc.mass(), test.mass())


if __name__ == '__main__':
    unittest.main()
