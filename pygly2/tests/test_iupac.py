import unittest

from pygly2.tests import common
from pygly2.io import iupac

monosaccharides = common.monosaccharides


class IUPACTests(unittest.TestCase):
    def test_monosaccharide_to_iupac(self):
        iupac_residue = iupac.to_iupac(monosaccharides.Fucose)
        reference = "a-L-Fucp"
        self.assertEqual(iupac_residue, reference)

    def test_monosaccharide_from_iupac(self):
        text = "a-L-Fucp"
        reference = monosaccharides.Fucose
        result = iupac.from_iupac(text)
        self.assertEqual(result, reference)

    def test_glycan_to_iupac(self):
        reference = 'a-L-Fucp-(1-6)-[a-D-Neu5Acp-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1-3)]b-D-Glcp2NAc-(1-6)-[a-D-Neu5Gcp-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1\
-3)]b-D-Glcp2NAc-(1-2)]a-D-Manp-(1-6)-[b-D-Glcp2NAc-(1-4)][b-D-Galp2NAc-(1-4)-b-D-Glcp2NAc-(1-4)-[a-D-Neu5Gcp-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp\
-(1-3)]b-D-Glcp2NAc-(1-2)]a-D-Manp-(1-3)]b-D-Manp-(1-4)-b-D-Glcp2NAc-(1-4)]?-D-Glcp2NAc'
        structure = common.load("complex_glycan")
        self.assertEqual(iupac.to_iupac(structure), reference)

    def test_glycan_from_iupac(self):
        text = 'a-L-Fucp-(1-6)-[a-D-Neu5Acp-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1-3)]b-D-Glcp2NAc-(1-6)-[a-D-Neu5Gcp-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1\
-3)]b-D-Glcp2NAc-(1-2)]a-D-Manp-(1-6)-[b-D-Glcp2NAc-(1-4)][b-D-Galp2NAc-(1-4)-b-D-Glcp2NAc-(1-4)-[a-D-Neu5Gcp-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp\
-(1-3)]b-D-Glcp2NAc-(1-2)]a-D-Manp-(1-3)]b-D-Manp-(1-4)-b-D-Glcp2NAc-(1-4)]?-D-Glcp2NAc'
        reference = common.load("complex_glycan")
        self.assertEqual(iupac.from_iupac(text), reference)
