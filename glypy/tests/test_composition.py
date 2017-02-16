import pickle
import unittest

from glypy.composition import composition, composition_transform
from glypy.structure import monosaccharide

from .common import load

ReducedEnd = monosaccharide.ReducedEnd


def make_composition_suite(composition_type):
    class CompositionTests(unittest.TestCase):
        def test_derivatize_bare(self):
            permethylated_reduced_mass = 1716.9033
            glycan = load("common_glycan")
            glycan.reducing_end = ReducedEnd()
            composition_transform.derivatize(glycan, 'methyl')
            self.assertAlmostEqual(glycan.mass(), permethylated_reduced_mass, 3)

        def test_strip_derivatize(self):
            glycan = load("common_glycan")
            glycan.reducing_end = ReducedEnd()
            mass = glycan.mass()
            composition_transform.derivatize(glycan, 'methyl')
            self.assertNotEqual(mass, glycan.mass())
            composition_transform.strip_derivatization(glycan)
            self.assertAlmostEqual(glycan.mass(), mass, 3)

        def test_composition_equality(self):
            self.assertEqual(composition_type("H2O"), composition_type("H2O"))
            self.assertEqual(composition_type("H2O") * 2, composition_type("(H2O)2"))

        def test_composition_substraction(self):
            self.assertEqual(
                composition_type("NH2O") - composition_type("N"), composition_type("H2O"))
            c = composition_type("NH2O")
            c -= composition_type("N")
            self.assertEqual(c, composition_type("H2O"))
            c += {"N": 1}
            self.assertEqual(c, composition_type("NH2O"))

        def test_isotope_parsing(self):
            self.assertFalse(composition_type("O[18]") == composition_type("O"))
            self.assertAlmostEqual(composition_type("O[18]").mass, 17.999, 3)
            self.assertRaises(
                composition.ChemicalCompositionError, lambda: composition_type("O[18.5]"))

        def test_inits(self):
            self.assertEqual(composition_type(O=1, H=2), composition_type(formula='H2O'))
            self.assertEqual(composition_type(O=1, H=2), composition_type("OH2"))
            self.assertEqual(composition_type("(N)(C[12]H3)2(H)"), composition_type("NC[12]2H6H"))
            self.assertRaises(
                composition.ChemicalCompositionError, lambda: composition_type("O2.2"))

        def test_operators(self):
            case = composition_type("H2O")
            self.assertEqual(case + {'O': 1}, composition_type("H2O2"))
            self.assertEqual({'O': 1} + case, composition_type("H2O2"))

            self.assertEqual(case - {'O': 1}, composition_type("H2"))
            self.assertEqual({'O': 1, "H": 4} - case, composition_type("H2"))

            self.assertEqual(case * 3, composition_type("H6O3"))
            self.assertEqual(case * -1, -case)
            self.assertRaises(
                composition.ChemicalCompositionError, lambda: case * 5.2)

        def test_massing(self):
            case = composition_type("H2O")
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

        def test_pickling(self):
            case = composition_type("H2O")
            self.assertEqual(case, pickle.loads(pickle.dumps(case)))

        def test_most_probable_isotopic_composition(self):
            case = composition_type("H2O")
            comp, probability = composition.most_probable_isotopic_composition(case)
            self.assertEqual(comp, {'H[1]': 2, 'O[16]': 1})
            self.assertAlmostEqual(probability, 0.997, 3)

    return CompositionTests


from glypy.composition.composition import PComposition
PCompositionTests = make_composition_suite(PComposition)
try:
    from glypy.composition.composition import CComposition
    CCompositionTests = make_composition_suite(CComposition)
except ImportError:
    pass

if __name__ == '__main__':
    unittest.main()
