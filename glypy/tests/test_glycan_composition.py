import unittest

import glypy
from glypy.composition import composition_transform
from glypy.composition import glycan_composition
from glypy import monosaccharides, Substituent, glycans

from common import load

water_mass = glypy.Composition("H2O").mass

GlycanComposition = glycan_composition.GlycanComposition


class IUPACLiteTests(unittest.TestCase):
    def test_parse_identity(self):
        items = ["Glc2NAc", "Neu5NAc", "Fuc", "6-dHex"]
        for item in items:
            term = glycan_composition.from_iupac_lite(item)
            deparse = glycan_composition.to_iupac_lite(term)
            self.assertEqual(item, deparse)


class MonosaccharideResidueTests(unittest.TestCase):
    def test_from_monosaccharide(self):
        for base in ["GlcNAc", "NeuAc", "Fuc"]:
            base = monosaccharides[base]
            residue = glycan_composition.MonosaccharideResidue.from_monosaccharide(base)
            self.assertAlmostEqual(base.mass() - water_mass, residue.mass(), 3)

    def test_clone_coherence(self):
        for base in ["GlcNAc", "NeuAc", "Fuc"]:
            base = monosaccharides[base]
            residue = glycan_composition.MonosaccharideResidue.from_monosaccharide(base)
            dup = residue.clone()
            self.assertAlmostEqual(dup.mass(), residue.mass(), 3)

    def test_derivatize(self):
        reference = {"GlcNAc": 245.1263, "Fuc": 174.0892, "NeuAc": 361.1737}
        for name in ["GlcNAc", "NeuAc", "Fuc"]:
            base = monosaccharides[name]
            residue = glycan_composition.MonosaccharideResidue.from_monosaccharide(base)
            composition_transform.derivatize(base, "methyl")
            composition_transform.derivatize(residue, "methyl")
            self.assertAlmostEqual(
                base.mass() - (water_mass + 2 * (Substituent("methyl").mass() - glypy.Composition("H").mass * 2)),
                residue.mass(), 3)
            self.assertAlmostEqual(reference[name], residue.mass(), 3)

            base = monosaccharides[name]
            composition_transform.derivatize(base, "methyl")
            residue = glycan_composition.MonosaccharideResidue.from_monosaccharide(base)

            self.assertAlmostEqual(
                base.mass() - (water_mass + 2 * (Substituent("methyl").mass() - glypy.Composition("H").mass * 2)),
                residue.mass(), 3)

            self.assertAlmostEqual(reference[name], residue.mass(), 3)


class GlycanCompositionTests(unittest.TestCase):
    def test_from_glycan(self):
        glyc = glycans["N-Linked Core"]
        comp = GlycanComposition.from_glycan(glyc)
        self.assertAlmostEqual(glyc.mass(), comp.mass(), 3)

    def test_derivatize(self):
        glyc = glycans["N-Linked Core"]
        composition_transform.derivatize(glyc, "methyl")
        comp = GlycanComposition.from_glycan(glyc)
        self.assertAlmostEqual(glyc.mass(), comp.mass(), 3)

        comp = GlycanComposition.from_glycan(glycans["N-Linked Core"])
        composition_transform.derivatize(comp, "methyl")
        self.assertAlmostEqual(glyc.mass(), comp.mass(), 3)

    def test_serialize(self):
        ref = '{Man:3; Glc2NAc:2}'
        glyc = glycans["N-Linked Core"]
        comp = GlycanComposition.from_glycan(glyc)
        self.assertEqual(comp.serialize(), ref)

        ref2 = '{Hex:3; Hex2NAc:2}'
        comp.drop_stems()
        self.assertEqual(comp.serialize(), ref2)

    def test_parse(self):
        ref = '{Man:3; Glc2NAc:2}'
        glyc = glycans["N-Linked Core"]
        comp = GlycanComposition.from_glycan(glyc)

        comp2 = GlycanComposition.parse(ref)
        self.assertEqual(comp, comp2)

    def test_update(self):
        ref = '{Man:3; Glc2NAc:2}'
        comp = GlycanComposition.parse(ref)
        self.assertEqual(comp["Man"], 3)
        comp.update({"Man": 1})
        self.assertEqual(comp["Man"], 1)
        comp.update(Man=2)
        self.assertEqual(comp["Man"], 2)

    def test_arithmetic(self):
        ref = '{Man:3; Glc2NAc:2}'
        comp = GlycanComposition.parse(ref)

        comp["Man"] += 5
        self.assertEqual(comp["Man"], 8)
        comp2 = GlycanComposition(Neu5NAc=2, Fuc=1, Man=1)
        self.assertEqual(comp2 + comp, GlycanComposition(Neu5NAc=2, Fuc=1, Man=9, Glc2NAc=2))
        comp -= {'Man': 5}
        self.assertEqual(comp, GlycanComposition(Man=4, Glc2NAc=2) - GlycanComposition(Man=1, Glc2NAc=0))
        comp *= 2
        self.assertEqual(comp, (GlycanComposition(Man=4, Glc2NAc=2) - GlycanComposition(Man=1, Glc2NAc=0)) * 2)

    def test_total_composition(self):
        ref = '{Man:3; Glc2NAc:2}'
        comp = GlycanComposition.parse(ref)
        glyc = glycans["N-Linked Core"]
        self.assertEqual(comp.total_composition(), glyc.total_composition())
