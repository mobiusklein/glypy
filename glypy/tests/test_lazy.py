import unittest

from glypy.utils import lazy
from glypy.structure.named_structures import MonosaccharideIndex


class ProxyTest(unittest.TestCase):
    def test_init(self):
        index = lazy.ProxyObject(MonosaccharideIndex)
        self.assertEquals(index._initializer, MonosaccharideIndex)
        self.assertEquals(index._source, None)
        fucose = index.Fucose
        self.assertTrue(isinstance(index._source, MonosaccharideIndex))
        index.dHex = fucose
        self.assertEquals(index.dHex, index.Fucose)
        index['ldeoxyGal'] = fucose
        self.assertEquals(index['ldeoxyGal'], index.Fucose)
