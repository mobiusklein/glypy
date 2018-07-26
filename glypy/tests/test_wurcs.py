import unittest

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


class TestWURCS(unittest.TestCase):
    def test_parse_wurcs(self):
        ref = glycoct.loads(G71237SD_glycoct)
        test = wurcs.loads(G71237SD_wurcs)
        self.assertEqual(ref, test)

    def test_write_wurcs(self):
        ref = glycoct.loads(G71237SD_glycoct)
        test = wurcs.loads(wurcs.dumps(ref))
        self.assertEqual(ref, test)


if __name__ == '__main__':
    unittest.main()
