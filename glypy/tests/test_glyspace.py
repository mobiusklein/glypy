import unittest
from glypy.io import glyspace
from glypy import tree, root


try:
    is_online = True
    response = glyspace.requests.get(glyspace.sparql_endpoint)
    response.raise_for_status()
except:
    is_online = False


def skip_not_online(fn):
    if not is_online:
        return unittest.skip("Not able to reach host")(fn)
    else:
        return fn


@skip_not_online
class GlySpaceClientTest(unittest.TestCase):
    def test_get(self):
        record = glyspace.get("G00034ND")
        self.assertEqual(record.has_primary_id, u'G00034ND')
        self.assertEqual(round(tree(record).mass()), 4104.0)
        self.assertEqual(record.exact_match.structure_, record.structure_)
        # matches = glyspace.get(glyspace.NSSKOS.exactMatch)
        # self.assertTrue(len(matches.exact_match) > 0)

    def test_from_taxon(self):
        self.assertTrue(len(glyspace.client.from_taxon(9606, limit=10)) == 10)
