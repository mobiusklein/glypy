import unittest

from pygly2.tests import common
from pygly2.io import iupac

monosaccharides = common.monosaccharides

class IUPACTests(unittest.TestCase):
    def test_to_iupac(self):
        iupac_residue = iupac.to_iupac(monosaccharides.Fucose)
        reference = "a-L-Fucp"
        self.assertEqual(iupac_residue, reference)
