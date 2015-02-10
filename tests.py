import unittest
import json
import logging

#logging.basicConfig(level="DEBUG")

import pygly2

from pygly2.structure import monosaccharide, constants, substituent, glycan, link, named_structures, structure_composition
from pygly2.io import glycoct
from pygly2.utils import enum, StringIO
from pygly2.composition import Composition, composition_transform


class GlycoCTParserTests(unittest.TestCase):
    _file_path = "./test_data/glycoct.txt"

    def test_parse_file(self):
        for g in glycoct.read(self._file_path):
            self.assertTrue(isinstance(g, glycan.Glycan))


# monosaccharide_masses = json.load(open("./test_data/monosaccharide_masses.json"))
monosaccharide_structures = json.load(open("./pygly2/structure/data/monosaccharides.json"))

wiki_masses = {
    "Iduronic Acid": 194.04,
    "Bacillosamine": 162.10,
    "Allose": 180.06,
}

#: Not common in any way other than
#: reused in many tests
common_glycan = '''
RES
1b:x-dglc-HEX-x:x
2b:b-dglc-HEX-1:5
3b:b-dglc-HEX-1:5
4b:b-dglc-HEX-1:5
5b:b-dglc-HEX-1:5
6b:b-dglc-HEX-1:5
7b:b-dglc-HEX-1:5
LIN
1:1o(6+1)2d
2:2o(3+1)3d
3:2o(6+1)4d
4:4o(6+1)5d
5:5o(3+1)6d'''


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
                    test_loc = test.replace('\n',' ')[i-10:i+10]
                    ref_loc = ref.replace('\n',' ')[i-10:i+10]
                    raise AssertionError("{j} != {k} at {i} in {name}\n{test_loc}\n{ref_loc}".format(**locals()))
                i += 1

    def test_ring_limit_modification(self):
        structure = named_structures.monosaccharides['Hex']
        self.assertRaises(IndexError, lambda : structure.add_modification('d', 8))

    def test_occupancy_limit_modification(self):
        structure = named_structures.monosaccharides['Hex']
        structure.add_modification('d', 4)
        self.assertRaises(ValueError, lambda : structure.add_modification('d', 4))

    def test_ring_limit_substituent(self):
        structure = named_structures.monosaccharides['Hex']
        self.assertRaises(IndexError, lambda : structure.add_substituent('methyl', 8))

    def test_occupancy_limit_substituent(self):
        structure = named_structures.monosaccharides['Hex']
        structure.add_substituent('methyl', 4)
        self.assertRaises(ValueError, lambda : structure.add_substituent('methyl', 4))

    def test_add_remove_modifcations(self):
        for name, mass in wiki_masses.items():
            structure = named_structures.monosaccharides[name]
            ref = structure.clone()
            self.assertAlmostEqual(ref.mass(), structure.mass())
            open_sites, unknowns = structure.open_attachment_sites()
            n_sites = len(open_sites)
            comp_delta = n_sites * structure_composition.modification_compositions[constants.Modification.aldi]
            for site in open_sites:
                structure.add_modification(constants.Modification.aldi, site)
            self.assertEqual(structure.total_composition(), ref.total_composition() + comp_delta)
            self.assertEqual([], structure.open_attachment_sites()[0])
            for site in open_sites:
                structure.drop_modification(site, constants.Modification.aldi)
            self.assertEqual(structure.total_composition(), ref.total_composition())


class GlycanTests(unittest.TestCase):
    _file_path = "./test_data/glycoct.txt"

    def test_from_glycoct(self):
        for glycan in glycoct.read(self._file_path):
            self.assertAlmostEqual(glycan.mass(), iter(glycoct.loads(glycan.to_glycoct())).next().mass())

    def test_fragments_preserve(self):
        for glycan in glycoct.read(self._file_path):
            dup = glycan.clone()
            self.assertEqual(glycan, dup)
            fragments = list(dup.fragments(kind='BY'))

            self.assertEqual(glycan, dup)

    def test_clone(self):
        glycan = glycoct.loads(common_glycan).next()
        ref = glycan.clone()
        glycan.reducing_end = 1        
        self.assertTrue(glycan != ref)

class CompositionTests(unittest.TestCase):

    def test_derivativize_bare(self):
        permethylated_reduced_mass = 1286.6718
        glycan = glycoct.loads(common_glycan).next()
        ref = glycan.clone()
        glycan.reducing_end = 1
        composition_transform.derivatize(glycan, 'methyl')
        self.assertAlmostEqual(glycan.mass(), permethylated_reduced_mass, 3)

    def test_composition_equality(self):
        self.assertEqual(Composition("H2O"), Composition("H2O"))

    def test_composition_substraction(self):
        self.assertEqual(Composition("NH2O") - Composition("N"), Composition("H2O"))

    def test_isotope_parsing(self):
        self.assertFalse(Composition("O[18]") == Composition("O"))
        self.assertAlmostEqual(Composition("O[18]").mass, 17.999, 3)




if __name__ == '__main__':
    unittest.main()
