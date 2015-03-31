import sys
import pdb
import functools
import traceback
import unittest
import json
from collections import defaultdict

import pygly2
from pygly2.structure import constants, substituent, glycan
from pygly2.structure import link, named_structures, structure_composition
from pygly2.structure import crossring_fragments
from pygly2.io import glycoct, linear_code
from pygly2.io.nomenclature import identity, synonyms
from pygly2.utils import StringIO, identity as ident_op, multimap, pickle, ET, enum
from pygly2.composition import Composition, composition_transform, composition
from pygly2.algorithms import subtree_search
from pygly2.algorithms import similarity

from common import emit

Substituent = pygly2.Substituent
Glycan = glycan.Glycan

def debug_on(*exceptions):
    if not exceptions:
        exceptions = (AssertionError, )

    def decorator(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            try:
                return f(*args, **kwargs)
            except exceptions:
                info = sys.exc_info()
                traceback.print_exception(*info)
                pdb.post_mortem(info[2])
        return wrapper
    return decorator
monosaccharides = named_structures.monosaccharides
glycans = named_structures.glycans

# monosaccharide_masses = json.load(open("./test_data/monosaccharide_masses.json"))
monosaccharide_structures = json.load(
    open("./pygly2/structure/data/monosaccharides.json"))

wiki_masses = {
    "Iduronic Acid": 194.04,
    "Bacillosamine": 162.10,
    "Allose": 180.06,
}

#: Not common in any way other than
#: reused in many tests
common_glycan = '''
RES
1b:b-dglc-HEX-1:5
2b:b-dgal-HEX-1:5
3b:b-dglc-HEX-1:5
4s:n-acetyl
5b:a-lgal-HEX-1:5|6:d
6b:b-dgal-HEX-1:5
7b:b-dglc-HEX-1:5
8s:n-acetyl
9b:a-lgal-HEX-1:5|6:d
10b:b-dgal-HEX-1:5
LIN
1:1o(4+1)2d
2:2o(3+1)3d
3:3d(2+1)4n
4:3o(3+1)5d
5:3o(4+1)6d
6:6o(3+1)7d
7:7d(2+1)8n
8:7o(3+1)9d
9:7o(4+1)10d'''

branchy_glycan = '''
RES
1b:x-dglc-HEX-x:x
2s:n-acetyl
3b:b-dman-HEX-1:5
4b:a-dman-HEX-1:5
5b:b-dglc-HEX-1:5
6s:n-acetyl
7b:b-dgal-HEX-1:5
8b:b-dglc-HEX-1:5
9s:n-acetyl
10b:b-dgal-HEX-1:5
11b:a-dman-HEX-1:5
12b:b-dglc-HEX-1:5
13s:n-acetyl
14b:b-dgal-HEX-1:5
LIN
1:1d(2+1)2n
2:1o(4+1)3d
3:3o(3+1)11d
4:3o(6+1)4d
5:4o(2+1)8d
6:4o(6+1)5d
7:5d(2+1)6n
8:5o(4+1)7d
9:8d(2+1)9n
10:8o(4+1)10d
11:11o(2+1)12d
12:12d(2+1)13n
13:12o(4+1)14d'''

broad_n_glycan = '''
RES
1b:b-dglc-HEX-1:5
2s:n-acetyl
3b:b-dglc-HEX-1:5
4s:n-acetyl
5b:b-dman-HEX-1:5
6b:a-dman-HEX-1:5
7b:b-dglc-HEX-1:5
8s:n-acetyl
9b:b-dgal-HEX-1:5
10b:a-lgal-HEX-1:5|6:d
11b:b-dglc-HEX-1:5
12s:n-acetyl
13b:b-dgal-HEX-1:5
14b:a-dman-HEX-1:5
15b:b-dglc-HEX-1:5
16s:n-acetyl
17b:b-dgal-HEX-1:5
18b:b-dglc-HEX-1:5
19s:n-acetyl
20b:b-dgal-HEX-1:5
LIN
1:1d(2+1)2n
2:1o(4+1)3d
3:3d(2+1)4n
4:3o(4+1)5d
5:5o(3+1)14d
6:5o(6+1)6d
7:6o(2+1)11d
8:6o(6+1)7d
9:7d(2+1)8n
10:7o(3+1)10d
11:7o(4+1)9d
12:11d(2+1)12n
13:11o(4+1)13d
14:14o(2+1)18d
15:14o(4+1)15d
16:15d(2+1)16n
17:15o(4+1)17d
18:18d(2+1)19n
19:18o(4+1)20d'''

sulfated_glycan = '''
RES
1b:o-dgal-HEX-0:0|1:aldi
2b:b-dglc-HEX-1:5
3s:n-acetyl
4b:b-dgal-HEX-1:5
5b:b-dglc-HEX-1:5
6s:n-acetyl
7b:b-dgal-HEX-1:5
8b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
9s:n-acetyl
10s:sulfate
11s:sulfate
12s:sulfate
LIN
1:1o(3+1)2d
2:2d(2+1)3n
3:2o(4+1)4d
4:4o(3+1)5d
5:5d(2+1)6n
6:5o(4+1)7d
7:7o(6+2)8d
8:8d(5+1)9n
9:5o(6+1)10n
10:4o(6+1)11n
11:2o(6+1)12n
'''

complex_glycan = '''
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
'''


class GlycoCTParserTests(unittest.TestCase):
    _file_path = "./test_data/glycoct.txt"

    def test_parse_file(self):
        for g in glycoct.read(self._file_path):
            self.assertTrue(isinstance(g, glycan.Glycan))


class NamedStructureTests(unittest.TestCase):

    def test_accessors(self):
        self.assertEqual(named_structures.monosaccharides.Fucose,
                         named_structures.monosaccharides["Fucose"])
        self.assertNotEqual(named_structures.monosaccharides.Man,
                            named_structures.monosaccharides.Hex)
        self.assertEqual(named_structures.glycans["N-Linked Core"],
                         named_structures.glycans["N-Linked Core"])


class StructureCompositionTests(unittest.TestCase):

    def test_missing_composition(self):
        import warnings
        warnings.filterwarnings('ignore')
        self.assertEqual(
            {}, structure_composition.monosaccharide_composition['hexadodecose'])
        warnings.filterwarnings('always')


class MonosaccharideTests(unittest.TestCase):
    _file_path = "./test_data/glycoct.txt"
    glycan = iter(glycoct.read(_file_path)).next()

    def test_from_glycoct(self):
        s = self.glycan.root.to_glycoct()
        b = StringIO(s)
        g = iter(glycoct.read(b)).next()
        self.assertEqual(g.root.to_glycoct(), s)

    def test_named_structure_masses(self):
        for name, mass in wiki_masses.items():
            structure = named_structures.monosaccharides[name]
            self.assertAlmostEqual(mass, structure.mass(), 2)

    def test_named_structure_glycoct(self):
        for name, glycoct_str in monosaccharide_structures.items():
            structure = named_structures.monosaccharides[name]
            test, ref = (structure.to_glycoct(), glycoct_str)
            i = 0
            for j, k in zip(test, ref):
                if j != k:
                    test_loc = test.replace('\n', ' ')[i - 10:i + 10]
                    ref_loc = ref.replace('\n', ' ')[i - 10:i + 10]
                    raise AssertionError(
                        "{j} != {k} at {i} in {name}\n{test_loc}\n{ref_loc}".format(**locals()))
                i += 1

    def test_ring_limit_modification(self):
        structure = named_structures.monosaccharides['Hex']
        self.assertRaises(
            IndexError, lambda: structure.add_modification('d', 8))

    def test_occupancy_limit_modification(self):
        structure = named_structures.monosaccharides['Hex']
        structure.add_modification('d', 4)
        self.assertRaises(
            ValueError, lambda: structure.add_modification('d', 4))

    def test_ring_limit_substituent(self):
        structure = named_structures.monosaccharides['Hex']
        self.assertRaises(
            IndexError, lambda: structure.add_substituent('methyl', 8))

    def test_occupancy_limit_substituent(self):
        structure = named_structures.monosaccharides['Hex']
        structure.add_substituent('methyl', 4)
        self.assertRaises(
            ValueError, lambda: structure.add_substituent('methyl', 4))

    def test_add_remove_modifcations(self):
        for name, mass in wiki_masses.items():
            structure = named_structures.monosaccharides[name]
            ref = structure.clone()
            self.assertAlmostEqual(ref.mass(), structure.mass())
            open_sites, unknowns = structure.open_attachment_sites()
            n_sites = len(open_sites)
            comp_delta = n_sites * \
                structure_composition.modification_compositions[
                    constants.Modification.d]()
            for site in open_sites:
                structure.add_modification(constants.Modification.d, site)
            self.assertEqual(
                structure.total_composition(), ref.total_composition() + comp_delta)
            self.assertEqual([], structure.open_attachment_sites()[0])
            for site in open_sites:
                structure.drop_modification(site, constants.Modification.d)
            self.assertEqual(
                structure.total_composition(), ref.total_composition())

    def test_add_remove_substituents(self):
        for name, mass in wiki_masses.items():
            structure = named_structures.monosaccharides[name]
            ref = structure.clone()
            self.assertAlmostEqual(ref.mass(), structure.mass())
            open_sites, unknowns = structure.open_attachment_sites()
            n_sites = len(open_sites)
            mass_delta = substituent.Substituent(
                'methyl').mass() * n_sites - Composition("H2").mass * n_sites
            ping = True
            for site in open_sites:
                if ping:
                    structure.add_substituent(
                        substituent.Substituent('methyl'), position=site)
                    ping = False
                else:
                    structure.add_substituent('methyl', position=site)
                    ping = True
            self.assertAlmostEqual(structure.mass(), ref.mass() + mass_delta)
            for site in open_sites:
                structure.drop_substituent(
                    site, substituent.Substituent('methyl'))
            self.assertAlmostEqual(structure.mass(), ref.mass())

    def test_validate_drop_substituents(self):
        structure = named_structures.monosaccharides['Hex']
        self.assertRaises(IndexError, lambda: structure.drop_substituent(99))
        self.assertRaises(IndexError, lambda: structure.drop_substituent(4))
        self.assertRaises(IndexError, lambda: structure.add_substituent(
            "n-acetyl", 3, max_occupancy=4).add_substituent(
            "methyl", 3, max_occupancy=4).drop_substituent(3, "n-glycolyl"))

    def test_add_remove_monosaccharides(self):
        for name, mass in wiki_masses.items():
            structure = named_structures.monosaccharides[name]
            ref = structure.clone()
            self.assertAlmostEqual(ref.mass(), glycan.Glycan(structure).mass())
            open_sites, unknowns = structure.open_attachment_sites()
            n_sites = len(open_sites)
            mass_delta = named_structures.monosaccharides[
                "Hex"].mass() * n_sites - Composition("H2O").mass * n_sites
            for site in open_sites:
                structure.add_monosaccharide(
                    named_structures.monosaccharides["Hex"], position=site, child_position=3)
            self.assertAlmostEqual(
                glycan.Glycan(structure).mass(), ref.mass() + mass_delta)
            for site in open_sites:
                structure.drop_monosaccharide(site)
            self.assertAlmostEqual(glycan.Glycan(structure).mass(), ref.mass())

    def test_validate_drop_monosacharide(self):
        structure = named_structures.monosaccharides['Hex']
        self.assertRaises(
            IndexError, lambda: structure.drop_monosaccharide(99))
        self.assertRaises(ValueError, lambda: structure.drop_monosaccharide(4))
        self.assertRaises(ValueError, lambda: structure.add_monosaccharide(
            named_structures.monosaccharides["Hex"], 3, max_occupancy=4).add_monosaccharide(
            named_structures.monosaccharides["Hex"], 3, max_occupancy=4).drop_monosaccharide(3))

    def test_validate_enums(self):
        structure = named_structures.monosaccharides['Hex']

        def t():
            structure.anomer = 5
        self.assertRaises(KeyError, t)

        def t():
            structure.stem = "gibberish"
        self.assertRaises(KeyError, t)

        def t():
            structure.superclass = "monose"
        self.assertRaises(KeyError, t)

        def t():
            structure.configuration = 'not-real'
        self.assertRaises(KeyError, t)

    def test_validate_reducing_end(self):
        structure = named_structures.monosaccharides['Hex']
        structure.reducing_end = 1
        structure.reducing_end = 3

        def t():
            structure.reducing_end = -1
        self.assertRaises(ValueError, t)

        def t():
            structure.reducing_end = 99
        self.assertRaises(IndexError, t)


# class GlycanTests(unittest.TestCase):
#     _file_path = "./test_data/glycoct.txt"
# 
#     def test_from_glycoct(self):
#         for structure in glycoct.read(self._file_path):
#             self.assertEqual(
#                 structure, glycoct.loads(structure.to_glycoct()).next())

#     def test_fragments_preserve(self):
#         structure = glycoct.loads(branchy_glycan).next()
#         dup = structure.clone()
#         self.assertEqual(structure, dup)
#         list(dup.fragments('AXZCBY', 2))

#         self.assertEqual(structure, dup)

#     def test_branch_counts(self):
#         structure = glycoct.loads(branchy_glycan).next()
#         self.assertEqual(structure.count_branches(), 3)

#     def test_fragments_mass(self):
#         structure = glycoct.loads(common_glycan).next()
#         frags = list(structure.fragments('ZCBY', 1))
#         import json
#         frags_ref = json.load(open('test_data/fragments-example.json'))
#         container = multimap.MultiMap()
#         for frag in frags_ref:
#             container[frag['kind']] = frag['mass']

#         def almost_equal(a, b):
#             e = 0.0001
#             return (a - e) <= b <= (a + e)

#         for frag in frags:
#             structure.name_fragment(frag)
#             kind = frag[0]
#             mass = frag[-1]
#             candidates = container[kind] + container[kind[::-1]]
#             if len(candidates) == 0:
#                 raise AssertionError("No candidates found for {}".format(frag))
#             res = (any(almost_equal(mass, x) for x in candidates))
#             if not res:
#                 raise AssertionError(
#                     "{} found no matches in {}".format(frag, candidates))

#     def test_disjoint_subtrees(self):
#         structure = glycoct.loads(common_glycan).next()
#         f = open('test_data/test_disjoint_subtrees.pkl', 'rb')
#         refs = {}
#         while True:
#             try:
#                 for ref in pickle.load(f):
#                     refs[tuple(ref.link_ids)] = ref
#             except EOFError:
#                 break
#         for tree in structure.substructures(3):
#             ref = refs[tuple(tree.link_ids)]
#             link_id_eq = tree.link_ids == ref.link_ids
#             parent_tree_eq = tree.parent_tree == ref.parent_tree
#             if not parent_tree_eq or not link_id_eq:
#                 raise AssertionError("DisjointTrees did not match")
#                 self.assertEqual(list(tree.parent_tree), list(ref.parent_tree))

#     def test_reducing_end(self):
#         structure = glycoct.loads(common_glycan).next()
#         self.assertEqual(structure.reducing_end, None)
#         structure.reducing_end = 1
#         self.assertEqual(structure.reducing_end, 1)

#     def test_clone(self):
#         structure = glycoct.loads(common_glycan).next()
#         ref = structure.clone()
#         structure.reducing_end = 1
#         self.assertTrue(structure != ref)

#     def test_indexing(self):
#         structure = glycoct.loads(common_glycan).next()
#         ref = structure.clone()
#         for i, node in enumerate(structure.index):
#             self.assertEqual(node.id, ref[i].id)
#         structure.deindex()
#         for i, node in enumerate(structure.index):
#             self.assertNotEqual(node.id, ref[i].id)
#         structure.index = None
#         self.assertRaises(IndexError, lambda: structure[0])

#     def test_traversal(self):
#         structure = glycoct.loads(common_glycan).next()
#         structure[-
#                   1].add_monosaccharide(named_structures.monosaccharides['NeuGc'])
#         structure.reindex(method='dfs')
#         ref = structure.clone()
#         self.assertEqual(structure[-1], ref[-1])
#         structure.reindex(method='bfs')
#         self.assertNotEqual(structure[-1], ref[-1])

#     def test_traversal_by_name(self):
#         structure = glycoct.loads(common_glycan).next()
#         structure[-
#                   1].add_monosaccharide(named_structures.monosaccharides['NeuGc'])
#         structure.reindex(method='dfs')
#         ref = structure.clone()
#         self.assertEqual(structure, ref)
#         structure.reindex(method='depth_first_traversal')
#         self.assertEqual(structure, ref)
#         self.assertRaises(
#             AttributeError, lambda: structure.reindex(method='not_real_traversal'))

#     def test_leaves(self):
#         structure = glycoct.loads(common_glycan).next()
#         leaves = list(structure.leaves())
#         for node in leaves:
#             self.assertTrue(len(list(node.children())) == 0)
#         leaves = list(structure.leaves(bidirectional=True))
#         for node in leaves:
#             self.assertTrue(
#                 len(list(node.children())) == 0 or node == structure.root)

#     def test_custom_traversal_method(self):
#         def rev_sort_dfs(self, visited=None, from_node=None, *args, **kwargs):
#             node_stack = list([self.root])
#             visited = set()
#             while len(node_stack) > 0:
#                 node = node_stack.pop()
#                 if node.id in visited:
#                     continue
#                 visited.add(node.id)
#                 yield (node)
#                 node_stack.extend(reversed(list(terminal for pos, link in node.links.items()
#                                                 for terminal in link if terminal.id not in visited and
#                                                 len(link.child.substituent_links) < 1)))
#         structure = glycoct.loads(common_glycan).next()
#         structure[-3].add_monosaccharide(
#             named_structures.monosaccharides['Hex'], 4).add_substituent('methyl', 5)
#         structure.reindex()
#         ref = structure.clone()
#         self.assertEqual(structure[-1], ref[-1])
#         structure.reindex(method=rev_sort_dfs)
#         self.assertNotEqual(structure[-1], ref[-1])

#     def test_topological_equality(self):
#         base = glycoct.loads(branchy_glycan).next()
#         a = base.clone()
#         b = base.clone()
#         c = base.clone()
#         self.assertEqual(base, b)
#         self.assertEqual(a, b)
#         list(a.leaves())[0].add_monosaccharide(monosaccharides["NeuGc"])
#         list(b.leaves())[1].add_monosaccharide(monosaccharides["NeuGc"])
#         list(c.leaves())[2].add_monosaccharide(monosaccharides["NeuGc"])
#         d = base.clone()
#         d_children = list(d.leaves())
#         d_children[0].add_monosaccharide(monosaccharides["NeuGc"])
#         d_children[1].add_monosaccharide(monosaccharides["NeuGc"])
#         self.assertTrue(a.topological_equality(b))
#         self.assertFalse(a.topological_equality(c))
#         self.assertFalse(a.topological_equality(base))
#         self.assertFalse(a.topological_equality(d))
#         self.assertFalse(b.topological_equality(d))

#     def test_substructures_does_not_mutate(self):
#         structure = glycoct.loads(broad_n_glycan).next()
#         ref = structure.clone()
#         for substructure in structure.substructures(max_cleavages=2):
#             pass
#         self.assertEqual(structure, ref)

#     def test_crossring(self):
#         structure = glycoct.loads(branchy_glycan).next()
#         frag_data = {
#             '0,2A': 485.174,
#             '0,2X': 1317.470,
#             '0,3A': 455.163,
#             '0,3X': 1347.481,
#             '0,4A': 425.153,
#             '0,4X': 1377.491,
#             '1,3A': 425.153,
#             '1,3X': 1377.491,
#             '1,4A': 455.163,
#             '1,4X': 1347.481,
#             '1,5A': 864.322,
#             '1,5X': 938.322,
#             '2,4X': 1742.623,
#             '2,5A': 469.179,
#             '2,5X': 1333.465,
#             '3,5A': 439.168,
#             '3,5X': 1363.4769
#         }
#         for fragment in structure.fragments("AX"):
#             if fragment.link_ids == [3]:
#                 self.assertAlmostEqual(frag_data[fragment.kind], fragment.mass, 2)

#     def test_subtree_from(self):
#         structure = glycoct.loads(branchy_glycan).next()
#         child = structure.root.children().next()[1]
#         subtree = Glycan.subtree_from(structure, child)
#         temp = structure.clone()
#         temproot = temp.root.children().next()[1]
#         for link in temp.root.links.values():
#             link.break_link(refund=True)
#         temp.root = temproot
#         self.assertEqual(temp, subtree)
#         self.assertEqual(Glycan.subtree_from(structure, 1), temp)


class CrossRingTests(unittest.TestCase):


    def test_unroll_ring(self):
        linear = crossring_fragments.unroll_ring(monosaccharides.GlcNAc)
        for i in range(len(linear)):
            self.assertEqual(linear[i]['backbone'], Composition("CH2O"))
            if i == 1:
                self.assertTrue(len(linear[i]['substituent_links']) > 0)

    class FragmentRecord(object):
        def __init__(self, name, mass, permethylated_mass, cleave_1, cleave_2, kind):
            self.name = name
            self.mass = float(mass)
            self.permethylated_mass = float(permethylated_mass)
            self.cleave_1 = int(cleave_1)
            self.cleave_2 = int(cleave_2)
            self.kind = kind
            self.composition = {}

        def __repr__(self):
            return "FragmentRecord({name} {mass})".format(**self.__dict__)

    def load_all(self):
        tree = ET.parse('./test_data/Cross-ring-data.xml')
        sugars = {}
        by_name = {}
        for cls in tree.iterfind(".//class"):
            sugar_class = {}
            for case in cls.iterfind(".//sugar"):
                isomers = [tag.attrib['abbr'] for tag in case.findall(".//isomer")]
                frags = defaultdict(dict)
                for frag in case.iterfind(".//fragment"):
                    data = [frag.attrib[n] for n in ["name", "mass_mono", "mass_pm_mono", "cleav1", "cleav2", "type"]]
                    rec = self.FragmentRecord(*data)
                    rec.composition = dict(map(
                        lambda x: (x.attrib['element'], x.attrib['count']), frag.findall(".//composition")))
                    frags[rec.cleave_1, rec.cleave_2][rec.kind] = rec
                sugar_class[case.attrib['subclass']] = frags
                for name in isomers:
                    by_name[name] = frags
            sugars[cls.attrib['name']] = sugar_class
        return by_name

    def test_monosacchide_crossring(self):
        test_cases = [
            'Glc',
            "GlcNAc",
            "Xyl",
            'GalA',
            'Fuc',
            'IdoA',
            "KDN"
        ]
        store = self.load_all()
        for case in test_cases:
            test = store[case]
            target = monosaccharides[case]
            for k, v in test.items():
                target_d = {t.kind: t for t in crossring_fragments.crossring_fragments(target, k[0], k[1])}
                target_d_permethylated = {t.kind: t for t in crossring_fragments.crossring_fragments(
                    composition_transform.derivatize(target.clone(), "methyl"), k[0], k[1])}
                for kind in {"A", "X"}:
                    self.assertAlmostEqual(v[kind].mass, target_d[kind].mass(), 3)
                    self.assertAlmostEqual(v[kind].permethylated_mass, target_d_permethylated[kind].mass(), 3)


class SubstituentTests(unittest.TestCase):

    def test_substituent_substituent(self):
        parent = substituent.Substituent("n-acetyl")
        pmass = parent.mass()
        child = substituent.Substituent("methyl")
        cmass = child.mass()

        parent.add_substituent(child, 2)
        self.assertAlmostEqual(
            parent.mass(), (pmass + cmass - Composition(H=2).mass))
        self.assertEqual(
            parent.total_composition(), parent.composition + child.composition)
        self.assertEqual(parent.children().next()[1], child)
        self.assertEqual(child.parents().next()[1], parent)
        parent.drop_substiuent(2)
        self.assertAlmostEqual(parent.mass(), pmass)

    def test_equality(self):
        sub_1 = substituent.Substituent("n-acetyl")
        sub_2 = substituent.Substituent("n-acetyl")
        self.assertEqual(sub_1, sub_2)
        parent = named_structures.monosaccharides['Hex']
        parent.add_substituent(sub_1, 3)
        self.assertNotEqual(sub_1, sub_2)

    def test_clone(self):
        parent = substituent.Substituent("n-acetyl")
        child = substituent.Substituent("methyl")

        parent.add_substituent(child, 2)

        dup = parent.clone()
        self.assertEqual(parent, dup)
        self.assertNotEqual(child, child.clone())


class MultiMapTests(unittest.TestCase):

    def test_iterators(self):
        from collections import Counter
        mm = multimap.MultiMap(a=1, b=3)
        mm['a'] = 3
        self.assertTrue(set(mm.keys()) == {'a', 'b'})
        self.assertTrue(set(mm.items()) == {('a', 1), ('a', 3), ('b', 3)})
        self.assertTrue(Counter(mm.values()) == Counter({1: 1, 3: 2}))

    def test_ordered_multimap(self):
        mm = multimap.MultiMap(a=1, b=3)
        mm['a'] = 3
        omm = multimap.OrderedMultiMap(a=1, b=3)
        omm['a'] = 3
        self.assertEqual(mm, omm)
        omm['c'] = 1
        self.assertNotEqual(mm, omm)
        self.assertNotEqual(omm, mm)
        self.assertFalse('c' in mm)


class CompositionTests(unittest.TestCase):

    def test_derivativize_bare(self):
        permethylated_reduced_mass = 1716.9033
        glycan = glycoct.loads(common_glycan).next()
        glycan.reducing_end = 1
        composition_transform.derivatize(glycan, 'methyl')
        self.assertAlmostEqual(glycan.mass(), permethylated_reduced_mass, 3)

    def test_composition_equality(self):
        self.assertEqual(Composition("H2O"), Composition("H2O"))
        self.assertEqual(Composition("H2O"), Composition("(H2O)"))

    def test_composition_substraction(self):
        self.assertEqual(
            Composition("NH2O") - Composition("N"), Composition("H2O"))

    def test_isotope_parsing(self):
        self.assertFalse(Composition("O[18]") == Composition("O"))
        self.assertAlmostEqual(Composition("O[18]").mass, 17.999, 3)
        self.assertRaises(
            composition.ChemicalCompositionError, lambda: Composition("O[18.5]"))

    def test_inits(self):
        self.assertEqual(Composition(O=1, H=2), Composition(formula='H2O'))
        self.assertEqual(Composition(O=1, H=2), Composition("OH2"))
        self.assertRaises(
            composition.ChemicalCompositionError, lambda: Composition("O2.2"))

    def test_operators(self):
        case = Composition("H2O")
        self.assertEqual(case + {'O': 1}, Composition("H2O2"))
        self.assertEqual({'O': 1} + case, Composition("H2O2"))

        self.assertEqual(case - {'O': 1}, Composition("H2"))
        self.assertEqual({'O': 1, "H": 4} - case, Composition("H2"))

        self.assertEqual(case * 3, Composition("H6O3"))
        self.assertEqual(case * -1, -case)
        self.assertRaises(
            composition.ChemicalCompositionError, lambda: case * 5.2)

    def test_massing(self):
        case = Composition("H2O")
        mono_mass = case.calc_mass()
        avg_mass = case.calc_mass(average=True)
        protonated = case + {"H+": 1}
        prot_mass = protonated.calc_mass()
        self.assertNotEqual(mono_mass, avg_mass)
        self.assertAlmostEqual(mono_mass, 18.0105, 3)
        self.assertAlmostEqual(avg_mass, 18.01528, 3)
        self.assertNotEqual(mono_mass, prot_mass)
        self.assertAlmostEqual(prot_mass, 19.01784, 3)
        self.assertAlmostEqual(case.calc_mass(charge=1), 19.01784, 3)
        self.assertRaises(
            composition.ChemicalCompositionError, lambda: protonated.calc_mass(charge=1))


class ConstantTests(unittest.TestCase):

    def test_translate(self):
        self.assertTrue(
            constants.Modification.d == constants.Modification['d'])
        self.assertTrue(constants.Modification.d == constants.Modification[
                        constants.Modification.d.value])

    def test_compare(self):
        self.assertTrue(constants.Modification.d == 'd')
        self.assertTrue(
            constants.Modification.d == constants.Modification.d.value)
        self.assertNotEqual(
            constants.Modification.d, constants.Stem[constants.Modification.d.value])
        self.assertRaises(KeyError, lambda: constants.SuperClass[1])
        self.assertNotEqual(constants.Modification.d, constants.Modification.a)

    def test_multiname(self):
        Modification = constants.Modification
        self.assertTrue(Modification.d == 'd')
        self.assertFalse(Modification.d == 'Deoxidation')
        # Add new name to existing EnumValue
        Modification.d.add_name("Deoxidation")
        self.assertTrue(Modification.d == 'Deoxidation')

        # Add new name to class by re-using EnumValue
        self.assertFalse(Modification.d == "deoxy")
        Modification.deoxy = Modification.d
        self.assertTrue(Modification.d == "deoxy")
        self.assertEqual(Modification.deoxy, Modification.Deoxidation)

        # Add new name to class and EnumValue by interpolating value
        self.assertFalse(Modification.deoxy == 'deoxy2')
        Modification.deoxy2 = Modification.d.value
        self.assertTrue(Modification.deoxy == 'deoxy2')
        self.assertRaises(KeyError, lambda: Modification.d.add_name('a'))

    def test_wrap_value(self):
        constants.Modification.not_a_real_modification = -5
        self.assertTrue(isinstance(constants.Modification.not_a_real_modification, enum.EnumValue))
        self.assertEqual(constants.Modification.not_a_real_modification.resolve(
            constants.Modification), -5)

    def test_resolve_failure(self):
        mapping = {'d': 5}
        self.assertTrue(constants.Modification.d.resolve(mapping) == 5)
        self.assertRaises(KeyError, lambda: constants.Modification.aldi.resolve(mapping))

    def test_instantiate_error(self):
        self.assertRaises(Exception, enum.Enum)

    def test_crossclass_inequality(self):
        class E1(enum.Enum):
            A = 1

        class E2(enum.Enum):
            A = 1
        self.assertNotEqual(E1.A, E2.A)


class LinkTests(unittest.TestCase):

    def test_link_equality(self):
        parent = named_structures.monosaccharides['Hex']
        child = named_structures.monosaccharides['Hex']
        other = named_structures.monosaccharides['Hex']
        link_1 = link.Link(parent, child, parent_position=3,
                           child_position=3, parent_loss='H', child_loss='OH')
        link_2 = link.Link(child, other, parent_position=6,
                           child_position=3, parent_loss='H', child_loss='OH')
        self.assertEqual(link_1, link_1)
        self.assertNotEqual(link_1, link_2)
        self.assertFalse(link_1 is None)

    def test_loss_composition(self):
        parent = named_structures.monosaccharides['Hex']
        child = named_structures.monosaccharides['Hex']

        link_1 = link.Link(parent, child, parent_position=3,
                           child_position=3, parent_loss='H', child_loss='OH')

        self.assertEqual(link_1.parent_loss, Composition(formula="H"))
        self.assertEqual(link_1.child_loss, Composition(formula="OH"))

    def test_break_and_reconnect(self):
        parent = named_structures.monosaccharides['Hex']
        child = named_structures.monosaccharides['Hex']

        link_1 = link.Link(parent, child, parent_position=3,
                           child_position=3, parent_loss='H', child_loss='OH')
        link_1.break_link(refund=True, reorient_fn=ident_op)
        self.assertTrue(len(parent.links[3]) == 0)
        self.assertTrue(len(child.links[3]) == 0)

        link_1.reconnect(refund=True, reorient_fn=ident_op)
        self.assertTrue(len(parent.links[3]) == 1)
        self.assertTrue(len(child.links[3]) == 1)

    def test_traversal(self):
        parent = named_structures.monosaccharides['Hex']
        child = named_structures.monosaccharides['Hex']

        link_1 = link.Link(parent, child, parent_position=3,
                           child_position=3, parent_loss='H', child_loss='OH')
        self.assertTrue(link_1.is_parent(parent))
        self.assertTrue(link_1.is_child(child))
        self.assertEqual(link_1.to(parent), child)
        self.assertEqual(link_1.to(child), parent)
        self.assertRaises(KeyError, lambda: link_1.to(1))


class SubtreeSearchTests(unittest.TestCase):

    def test_subtree_inclusion(self):
        core = glycans['N-Linked Core']
        tree = glycoct.loads(broad_n_glycan).next()
        self.assertTrue(subtree_search.subtree_of(core, tree))
        self.assertTrue(subtree_search.subtree_of(tree, core) is None)

    def test_maximum_common_subtree(self):
        core = glycans['N-Linked Core']
        tree = glycoct.loads(branchy_glycan).next()
        res = subtree_search.maximum_common_subgraph(core, tree)
        self.assertEqual(res.score, 6.2)

    def test_is_n_glycan(self):
        core = glycans['N-Linked Core']
        tree = glycoct.loads(broad_n_glycan).next()
        result = (subtree_search.subtree_of(core, tree))
        self.assertTrue(result == 1)
        tree = glycoct.loads(complex_glycan).next()
        result = (subtree_search.subtree_of(core, tree, exact=False))
        self.assertTrue(result == 1)
        tree = glycoct.loads(branchy_glycan).next()
        result = (subtree_search.subtree_of(core, tree, exact=False))
        self.assertTrue(result is None)


class IdentifyTests(unittest.TestCase):

    def test_is_a_predicate(self):
        for name, monosaccharide in monosaccharides.items():
            self.assertTrue(identity.is_a(monosaccharide, name))

    @debug_on()
    def test_identify_as(self):
        for name, monosaccharide in monosaccharides.items():
            if monosaccharide == monosaccharides.Hex:
                continue
            pref_name = identity.identify(monosaccharide)
            if not (name == pref_name or name in synonyms.monosaccharides[pref_name]):
                raise AssertionError(
                    "{}".format((name, pref_name, synonyms.monosaccharides[pref_name])))

    def test_identify_substituents(self):
        self.assertTrue(
            identity.is_a(Substituent("n-acetyl"), Substituent('n-acetyl')))
        self.assertFalse(
            identity.is_a(Substituent('methyl'), Substituent('n-acetyl')))
        self.assertFalse(
            identity.is_a(monosaccharides.Man, Substituent('n-acetyl')))
        self.assertFalse(
            identity.is_a(Substituent('n-acetyl'), monosaccharides.Man))

    def test_get_preferred_name(self):
        self.assertTrue(identity.get_preferred_name('bdMan') == 'Man')

    def test_identify_failure(self):
        # Will fail because Hex is blacklisted
        self.assertRaises(
            identity.IdentifyException, lambda: identity.identify(monosaccharides.Hex))


class LinearCodeTests(unittest.TestCase):

    def test_translate(self):
        broad = glycoct.loads(broad_n_glycan).next()
        dup = linear_code.loads(linear_code.dumps(broad))
        self.assertEqual(broad, dup)

        # linear code doesn't know about modifications or
        # ring shape
        sulfated = glycoct.loads(sulfated_glycan).next()
        sulfated.root.drop_modification(1, 'aldi')
        sulfated.root.ring_start = 1
        sulfated.root.ring_end = 5
        dup = linear_code.loads(linear_code.dumps(sulfated))
        self.assertEqual(dup, sulfated)

        sulfated = glycoct.loads(sulfated_glycan).next()
        dup = linear_code.loads(linear_code.dumps(sulfated))
        self.assertNotEqual(sulfated, dup)


class SimilarityTests(unittest.TestCase):

    def test_deep_similarity(self):
        branchy = glycoct.loads(branchy_glycan).next()
        broad = glycoct.loads(broad_n_glycan).next()
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

if __name__ == '__main__':
    unittest.main()
