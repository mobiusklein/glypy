import unittest

from glypy import monosaccharides, Substituent
from glypy.io.nomenclature import synonyms, identity


class IdentifyTests(unittest.TestCase):

    def test_is_a_predicate(self):
        for name, mono in monosaccharides.items():
            result = identity.is_a(mono, name)
            self.assertTrue(result)

    def test_identify_as(self):
        for name, mono in monosaccharides.items():
            if name in {"Hex", "Pen", "Oct", "Hep", "Non"}:
                continue
            pref_name = identity.identify(mono)
            if monosaccharides[pref_name] == mono:
                if not (name == pref_name or name in synonyms.monosaccharides[pref_name]):
                    raise AssertionError(
                        "{}".format((name, pref_name, synonyms.monosaccharides[pref_name])))

    def test_identify_index(self):
        identifier = identity.MonosaccharideIdentifier()
        for name, mono in monosaccharides.items():
            pref_name = identifier.identify(mono)
            if monosaccharides[pref_name] == mono:
                if not (name == pref_name or name in synonyms.monosaccharides[pref_name]):
                    raise AssertionError(
                        "{}".format((name, pref_name, synonyms.monosaccharides[pref_name])))

    def test_identify_substituents(self):
        self.assertTrue(
            identity.is_a(Substituent("n-acetyl"), Substituent("n-acetyl")))
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

    def test_precision(self):
        self.assertFalse(identity.is_a(monosaccharides.Kdn, monosaccharides.NeuAc))
        self.assertFalse(identity.is_a(monosaccharides.NeuAc, monosaccharides.Kdn))

    def test_grouping_axes(self):
        tree = identity.residue_list_to_tree((monosaccharides).values())


if __name__ == '__main__':
    unittest.main()
