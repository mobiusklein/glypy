import os
import unittest
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from glypy import monosaccharides, glycans
from glypy.io import glycoct
from glypy import plot


def outplot_dir():
    wd = os.getcwd()
    if "test_data" in os.listdir(wd):
        return "test_data" + os.sep
    elif "test_data" in os.listdir(wd + os.sep + '..'):
        return wd + os.sep + '..' + os.sep + "test_data"
    else:
        os.mkdir("test_data")
        return "test_data" + os.sep

plot_dir = outplot_dir()


glycan_structure = glycoct.loads('''
RES
1b:x-dglc-HEX-1:5
2s:n-acetyl
3b:b-dglc-HEX-1:5
4s:n-acetyl
5b:b-dman-HEX-1:5
6b:a-dman-HEX-1:5
7b:b-dglc-HEX-1:5
8s:n-acetyl
9b:a-lgal-HEX-1:5|6:d
10b:b-dgal-HEX-1:5
11b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
12s:n-glycolyl
13b:b-dglc-HEX-1:5
14s:n-acetyl
15b:b-dgal-HEX-1:5
16s:n-acetyl
17b:b-dglc-HEX-1:5
18s:n-acetyl
19b:a-dman-HEX-1:5
20b:b-dglc-HEX-1:5
21s:n-acetyl
22b:a-lgal-HEX-1:5|6:d
23b:b-dgal-HEX-1:5
24b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
25s:n-glycolyl
26b:b-dglc-HEX-1:5
27s:n-acetyl
28b:a-lgal-HEX-1:5|6:d
29b:b-dgal-HEX-1:5
30b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
31s:n-acetyl
32b:a-lgal-HEX-1:5|6:d
LIN
1:1d(2+1)2n
2:1o(4+1)3d
3:3d(2+1)4n
4:3o(4+1)5d
5:5o(3+1)6d
6:6o(2+1)7d
7:7d(2+1)8n
8:7o(3+1)9d
9:7o(4+1)10d
10:10o(3+2)11d
11:11d(5+1)12n
12:6o(4+1)13d
13:13d(2+1)14n
14:13o(4+1)15d
15:15d(2+1)16n
16:5o(4+1)17d
17:17d(2+1)18n
18:5o(6+1)19d
19:19o(2+1)20d
20:20d(2+1)21n
21:20o(3+1)22d
22:20o(4+1)23d
23:23o(3+2)24d
24:24d(5+1)25n
25:19o(6+1)26d
26:26d(2+1)27n
27:26o(3+1)28d
28:26o(4+1)29d
29:29o(3+2)30d
30:30d(5+1)31n
31:1o(6+1)32d
''')


class PlotTests(unittest.TestCase):
    def test_plot_monosaccharide(self):
        for name, residue in monosaccharides.items():
            plot.plot(residue)
            plt.clf()
            plt.close('all')

    def test_draw_tree(self):
        dt = plot.DrawTree(glycan_structure)
        self.assertEqual(dt.get(1).tree, glycan_structure.get(1))
        parent, child = dt.get_link_pair(1)
        parent = parent.tree
        child = child.tree
        self.assertEqual(parent, glycan_structure.get(1))

    def test_plot_iupac(self):
        plot.plot(glycan_structure, symbol_nomenclature='iupac')

    def test_plot_glycans(self):
        for name, structure in glycans.items():
            for layout in ['balanced', 'topological']:
                for symbol_nomenclature in ['cfg']:  # ['cfg', 'iupac']:
                    plot.plot(
                        structure, label=True, orientation='h',
                        symbol_nomenclature=symbol_nomenclature, layout=layout)
                    path = "%s/%s_%s_%s.png" % (plot_dir, name, layout, symbol_nomenclature)
                    plt.savefig(path)
                    plt.clf()
                    plt.close('all')

        plt.figure()
        plot.enumerate_tree(*plot.plot(glycan_structure))
        plt.savefig(plot_dir + os.sep + "enum_tree" + '.png')


if __name__ == '__main__':
    unittest.main()
