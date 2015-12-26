import unittest

import glypy

from glypy.io import glycomedb
from glypy import tree, root


try:
    is_online = True
    response = glycomedb.get(576)
    response.raise_for_status()
except:
    is_online = False


def skip_not_online(fn):
    if not is_online:
        return unittest.skip("Not able to reach host")(fn)
    else:
        return fn


class GlycomeDBClientTests(unittest.TestCase):

    @classmethod
    def setupClass(cls):
        glycomedb.set_cache(":memory:")

    def test_get(self):
        case = glycomedb.get(576)
        self.assertAlmostEqual(case.mass(), 2516.914472, 3)

    def test_get_record(self):
        case = glycomedb.get_record(576)
        self.assertAlmostEqual(case.mass(), 2516.914472, 3)
        self.assertTrue(case.is_n_glycan)
        self.assertEqual(case.id, 576)
        self.assertEqual(case.dbxref[0].id, "719")
        self.assertEqual(case.taxa[0].tax_id, '9913')
        self.assertEqual(root(case), root(glycomedb.get(576)))

    def test_search_by_species(self):
        r = glycomedb.search_by_species(9606)
        g = r.next()
        self.assertEqual(g.id, 28)
        self.assertEqual(g.num_species, '3')
        case = g.get()
        self.assertAlmostEqual(case.mass(), 709.26405752285, 3)
        i = 0
        for _ in r:
            i += 1
            if i > 50:
                break
        self.assertGreater(i, 1)

    @unittest.skip(True)
    def test_minimum_common_substructure(self):
        case = glycomedb.get(576)
        g = glycomedb.search_minimum_common_substructure(case).next()
        self.assertEqual(g.id, 576)
        self.assertEqual(g.score, '20,000')

if __name__ == '__main__':
    unittest.main()
