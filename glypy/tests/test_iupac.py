import unittest


import glypy
from glypy.composition import composition_transform
from glypy.tests import common
from glypy.io import iupac, glycoct

monosaccharides = common.monosaccharides


class IUPACTests(unittest.TestCase):
    def test_monosaccharide_to_iupac(self):
        iupac_residue = iupac.to_iupac(monosaccharides.Fucose)
        reference = "a-L-Fucp"
        self.assertEqual(iupac_residue, reference)

    def test_hexa(self):
        text = "a-D-HexpA"
        result = iupac.from_iupac(text)
        self.assertAlmostEqual(result.mass(), 194.042, 2)

    def test_autoname_substituent(self):
        h = monosaccharides.Hex
        h.add_substituent("phosphate")
        iupac.substituents_map_to.pop("phosphate")
        t = iupac.to_iupac(h)
        r = iupac.from_iupac(t)
        self.assertEqual(h, r)

    def test_special_cases(self):
        for text, mass in [('?-D-Neup', 267.095), ('?-D-Neup5Ac', 309.105), ('a-D-Neup5NGc', 325.100)]:
            self.assertAlmostEqual(iupac.from_iupac(
                iupac.to_iupac(iupac.from_iupac(text))).mass(), mass, 2)

    def test_ring_type(self):
        for text, start, stop in [("a-D-Hexp", 1, 5), ("a-D-Hexf", 1, 4), ("a-D-Hexo", 0, 0)]:
            mono = iupac.from_iupac(text)
            self.assertEqual(mono.ring_start, start)
            self.assertEqual(mono.ring_end, stop)

    def test_alternate_superclass(self):
        text = "a-D-2-deoxy-araHex"
        mono = iupac.from_iupac(text)
        self.assertEqual(mono.stem[0], 'ara')

    def test_monosaccharide_from_iupac(self):
        text = "a-L-Fucp"
        reference = monosaccharides.Fucose
        result = iupac.from_iupac(text)
        self.assertEqual(result, reference)

    def test_glycan_to_iupac(self):
        reference = 'a-L-Fucp-(1-6)-[a-D-Neup5Ac-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1-3)]b-D-Glcp2NAc-(1-6)-[a-D-Neup5Gc-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1-3)]b-D-Glcp2NAc-(1-2)]a-D-Manp-(1-6)-[b-D-Glcp2NAc-(1-4)][b-D-Galp2NAc-(1-4)-b-D-Glcp2NAc-(1-4)-[a-D-Neup5Gc-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1-3)]b-D-Glcp2NAc-(1-2)]a-D-Manp-(1-3)]b-D-Manp-(1-4)-b-D-Glcp2NAc-(1-4)]?-D-Glcp2NAc'
        structure = common.load("complex_glycan")
        self.assertEqual(iupac.to_iupac(structure), reference)

    def test_glycan_from_iupac(self):
        text = 'a-L-Fucp-(1-6)-[a-D-Neup5Ac-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1-3)]b-D-Glcp2NAc-(1-6)-[a-D-Neup5Gc-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1-3)]b-D-Glcp2NAc-(1-2)]a-D-Manp-(1-6)-[b-D-Glcp2NAc-(1-4)][b-D-Galp2NAc-(1-4)-b-D-Glcp2NAc-(1-4)-[a-D-Neup5Gc-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1-3)]b-D-Glcp2NAc-(1-2)]a-D-Manp-(1-3)]b-D-Manp-(1-4)-b-D-Glcp2NAc-(1-4)]?-D-Glcp2NAc'
        reference = common.load("complex_glycan")
        self.assertEqual(glycoct.canonicalize(iupac.from_iupac(text)), glycoct.canonicalize(reference))

    def test_glycan_from_iupac_file_like(self):
        text = '''>test
        a-L-Fucp-(1-6)-[a-D-Neup5Ac-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1-3)]b-D-Glcp2NAc-(1-6)-[a-D-Neup5Gc-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1
        -3)]b-D-Glcp2NAc-(1-2)]a-D-Manp-(1-6)-[b-D-Glcp2NAc-(1-4)][b-D-Galp2NAc-(1-4)-b-D-Glcp2NAc-(1-4)-[a-D-Neup5Gc-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp
        -(1-3)]b-D-Glcp2NAc-(1-2)]a-D-Manp-(1-3)]b-D-Manp-(1-4)-b-D-Glcp2NAc-(1-4)]?-D-Glcp2NAc'''
        reader = iupac.IUPACParser.loads(text, 'fasta')
        structure = next(reader)

        text = 'a-L-Fucp-(1-6)-[a-D-Neup5Ac-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1-3)]b-D-Glcp2NAc-(1-6)-[a-D-Neup5Gc-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1-3)]b-D-Glcp2NAc-(1-2)]a-D-Manp-(1-6)-[b-D-Glcp2NAc-(1-4)][b-D-Galp2NAc-(1-4)-b-D-Glcp2NAc-(1-4)-[a-D-Neup5Gc-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1-3)]b-D-Glcp2NAc-(1-2)]a-D-Manp-(1-3)]b-D-Manp-(1-4)-b-D-Glcp2NAc-(1-4)]?-D-Glcp2NAc'
        equiv = next(iupac.IUPACParser.loads(text, 'line'))
        self.assertEqual(equiv, structure)


class DerivatizationAwareIUPACTests(unittest.TestCase):
    def test_monosaccharide_parse(self):
        text = '?-?-Hexp2NAc^Me'
        parser = iupac.DerivatizationAwareMonosaccharideDeserializer()
        obj, _ = parser(text)
        ref = composition_transform.derivatize(glypy.monosaccharides.HexNAc, 'methyl')
        self.assertEqual(obj, ref)

    def test_monosaccharide_serialize(self):
        obj = composition_transform.derivatize(glypy.monosaccharides.HexNAc, 'methyl')
        ref = '?-?-Hexp2NAc^Me'
        serializer = iupac.DerivatizationAwareMonosaccharideSerializer()
        text = serializer(obj)
        self.assertEqual(text, ref)

    def test_derivatized_glycan_parse(self):
        ref = composition_transform.derivatize(
            glypy.motifs["N-Glycan complex 1"], "methyl")
        serializer = iupac.GlycanSerializer(iupac.DerivatizationAwareMonosaccharideSerializer())
        text = serializer(ref)
        deserializer = iupac.GlycanDeserializer(iupac.DerivatizationAwareMonosaccharideDeserializer())
        obj = deserializer(text)
        self.assertEqual(obj, ref)

    def test_compatible_with_non_derivatized(self):
        serializer = iupac.GlycanSerializer(iupac.DerivatizationAwareMonosaccharideSerializer())
        deserializer = iupac.GlycanDeserializer(iupac.DerivatizationAwareMonosaccharideDeserializer())
        self.assertEqual(
            deserializer(
                serializer(
                    glypy.motifs["N-Glycan complex 1"])),
            glypy.motifs["N-Glycan complex 1"])


if __name__ == '__main__':
    unittest.main()
