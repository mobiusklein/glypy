import unittest
import hjson

from glypy.structure import named_structures, constants, monosaccharide, substituent, glycan
from glypy.composition import structure_composition, Composition, composition_transform
from glypy.io import glycoct

from .common import StringIO, load, pickle


Monosaccharide = monosaccharide.Monosaccharide


monosaccharide_structures = hjson.load(
    open("./glypy/structure/data/monosaccharides.hjson"))

wiki_masses = {
    "Iduronic Acid": 194.04,
    "Bacillosamine": 162.10,
    "Allose": 180.06,
    "Neu5Ac": 309.105,
    "Fru": 180.06,
}

ReducedEnd = monosaccharide.ReducedEnd


class MonosaccharideTests(unittest.TestCase):
    _file_path = "./test_data/glycoct.txt"
    glycan = next(iter(glycoct.read(_file_path)))

    def test_depth(self):
        self.assertEqual(monosaccharide.depth(load("complex_glycan").root), 7)

    def test_from_glycoct(self):
        s = self.glycan.root.serialize('glycoct')
        b = StringIO(s)
        g = next(iter(glycoct.read(b)))
        self.assertEqual(g.root.serialize('glycoct'), s)

    def test_named_structure_masses(self):
        for name, mass in wiki_masses.items():
            structure = named_structures.monosaccharides[name]
            self.assertAlmostEqual(mass, structure.mass(), 2)

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
        structure.reducing_end = True
        composition_transform.derivatize(structure, "methyl")
        self.assertEqual(structure.mass(), structure.clone().mass())

    def test_low_level_traverse(self):
        branchy = load("branchy_glycan")
        t1 = monosaccharide.traverse(branchy.root)
        t2 = branchy.iternodes()
        for a, b in zip(list(t1), list(t2)):
            self.assertEqual(a, b)

    def test_low_level_graph_clone(self):
        branchy = load("branchy_glycan")
        self.assertEqual(branchy.root, monosaccharide.graph_clone(branchy.root))

    def test_ring_shape(self):
        hexose = named_structures.monosaccharides.Hex
        self.assertEqual(hexose.ring_type, "pyranose")
        hexose.ring_start += 1
        self.assertEqual(hexose.ring_type, "furanose")
        hexose.ring_start = 0
        hexose.ring_end = 0
        self.assertEqual(hexose.ring_type, "open")
        hexose.ring_end = None
        self.assertEqual(hexose.ring_type, "x")

    def test_exact_ordering_equality(self):
        hexose = named_structures.monosaccharides.Hex
        self.assertEqual(hexose, hexose)
        res1, res2 = hexose.clone(), hexose.clone()

        self.assertEqual(res1, res2)
        res1.add_substituent("n_acetyl", 2)
        res2.add_substituent("n_acetyl", 2)
        self.assertEqual(res1, res2)
        res2.drop_substituent(2)
        res2.add_substituent("n_acetyl", 3)
        self.assertNotEqual(res1, res2)

        res3 = res1.clone()
        res1.add_monosaccharide(hexose.clone(), 6)
        res3.add_monosaccharide(hexose.clone(), 6)
        self.assertEqual(res1, res3)

        res3.drop_monosaccharide(6)
        res3.add_monosaccharide(hexose.clone(), 5)
        self.assertNotEqual(res1, res3)

    def test_children(self):
        n_core = named_structures.glycans["N-Linked Core"]
        n_core.set_reducing_end(True)
        for pos, child in n_core.root.children():
            self.assertNotEqual(type(child), substituent.Substituent)
        for pos, child in n_core.root.reducing_end.children():
            self.assertNotEqual(type(child), Monosaccharide)

    def test_total_composition(self):
        hexose = named_structures.monosaccharides.Hex
        self.assertEqual(hexose.total_composition(), {"C": 6, "O": 6, "H": 12})
        hexose.reducing_end = True
        self.assertEqual(hexose.total_composition(), {"C": 6, "O": 6, "H": 14})

    def test_serialize_does_not_mutate(self):
        ref = self.glycan.clone()
        ids = [node.id for node in ref]
        for i, node in enumerate(ref):
            self.assertEqual(ids[i], node.id)
            self.assertEqual(node, self.glycan[i])
            str(node)
            self.assertEqual(ids[i], node.id)
            self.assertEqual(node, self.glycan[i])

    def test_pickle(self):
        for mono in named_structures.monosaccharides.values():
            self.assertEqual(mono, pickle.loads(pickle.dumps(mono)))


class ReducedEndTests(unittest.TestCase):
    def test_equality(self):
        self.assertEqual(ReducedEnd(), ReducedEnd())
        self.assertNotEqual(ReducedEnd(), ReducedEnd("H[2]H"))


if __name__ == '__main__':
    unittest.main()
