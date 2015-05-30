import sys
import pdb
import functools
import traceback
import unittest
import json
import itertools
from collections import defaultdict

import glypy
from glypy.structure import constants, substituent, glycan, monosaccharide
from glypy.structure import link, named_structures, structure_composition
from glypy.structure import crossring_fragments
from glypy.io import glycoct, linear_code
from glypy.io.nomenclature import identity, synonyms
from glypy.utils import StringIO, identity as ident_op, multimap, pickle, ET, enum
from glypy.composition import Composition, composition_transform
from glypy.algorithms import subtree_search
from glypy.algorithms import similarity

from common import emit, load

ReducedEnd = monosaccharide.ReducedEnd
Substituent = glypy.Substituent
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
    open("./glypy/structure/data/monosaccharides.json"))

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

    def test_parse_cyclical(self):
        structure = load("cyclical_glycan")
        self.assertAlmostEqual(structure.mass(), 810.2641170925)

    def test_parse_repeating(self):
        structure = load("repeating_glycan")
        self.assertAlmostEqual(structure.mass(), 678.18547284898)


class NamedStructureTests(unittest.TestCase):

    def test_accessors(self):
        self.assertEqual(named_structures.monosaccharides.Fucose,
                         named_structures.monosaccharides["Fucose"])
        self.assertNotEqual(named_structures.monosaccharides.Man,
                            named_structures.monosaccharides.Hex)
        self.assertEqual(named_structures.glycans["N-Linked Core"],
                         named_structures.glycans["N-Linked Core"])
        self.assertEqual(named_structures.motifs["N-Glycan hybrid 1"],
                         named_structures.motifs["N-Glycan hybrid 1"])


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
        composition = structure.total_composition()
        structure.reducing_end = ReducedEnd()
        self.assertEqual(structure.total_composition(), composition + Composition("H2"))
        structure.reducing_end = True
        self.assertEqual(structure.total_composition(), composition + Composition("H2"))
        self.assertEqual(structure.total_composition(), structure.clone().total_composition())
        self.assertEqual(structure.total_composition(), pickle.loads(pickle.dumps(structure)).total_composition())
        structure.reducing_end = None
        self.assertEqual(structure.total_composition(), composition)

    def test_low_level_traverse(self):
        branchy = load("branchy_glycan")
        t1 = monosaccharide.traverse(branchy.root)
        t2 = branchy.dfs()
        for a, b in itertools.izip(list(t1), list(t2)):
            self.assertEqual(a, b)

    def test_low_level_graph_clone(self):
        branchy = load("branchy_glycan")
        self.assertEqual(branchy.root, monosaccharide.graph_clone(branchy.root))


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
        parent.drop_substituent(2)
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

    def test_derivatize_pathway_setter(self):
        substituent.Substituent("not-real", None, Composition("H2O"), can_nh_derivatize=False, is_nh_derivatizable=True)
        self.assertRaises(KeyError, lambda: substituent.Substituent("not-real", None, Composition("H2O")))
        substituent.DerivatizePathway.register("not-real", False, True)
        substituent.Substituent("not-real", None, Composition("H2O"))


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
        link_1.break_link(refund=True)
        self.assertTrue(len(parent.links[3]) == 0)
        self.assertTrue(len(child.links[3]) == 0)

        link_1.reconnect(refund=True)
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
        sulfated.reducing_end = None
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
