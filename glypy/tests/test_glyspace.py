import unittest
from glypy.io import glyspace
from glypy import tree, root

from pytest import mark

try:
    is_online = True
    response = glyspace.requests.get(glyspace.sparql_endpoint)
    response.raise_for_status()
except:
    is_online = False


def skip_not_online(fn):
    if not is_online:
        return mark.xfail(mark.skip("Not able to reach host")(fn))
    else:
        return fn


example_n3 = '''@prefix dcterms: <http://purl.org/dc/terms/> .
@prefix glycan: <http://purl.jp/bio/12/glyco/glycan#> .
@prefix glycoinfo: <http://rdf.glycoinfo.org/glycan/> .
@prefix glycomedb: <http://rdf.glycome-db.org/glycan/> .
@prefix glytoucan: <http://www.glytoucan.org/glyco/owl/glytoucan#> .
@prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .
@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .
@prefix skos: <http://www.w3.org/2004/02/skos/core#> .
@prefix xml: <http://www.w3.org/XML/1998/namespace> .
@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .

glycoinfo:G80903UK a glycan:saccharide ;
    glycan:has_glycosequence <http://rdf.glycoinfo.org/glycan/G80903UK/glycoct>,
        <http://rdf.glycoinfo.org/glycan/G80903UK/iupac_condensed>,
        <http://rdf.glycoinfo.org/glycan/G80903UK/iupac_extended>,
        <http://rdf.glycoinfo.org/glycan/G80903UK/wurcs/2.0> ;
    glycan:has_motif glycoinfo:G00031MO,
        glycoinfo:G00032MO,
        glycoinfo:G00033MO,
        glycoinfo:G00034MO,
        glycoinfo:G00041MO,
        glycoinfo:G00042MO ;
    glycan:has_resource_entry <http://rdf.glycoinfo.org/glycan/resource-entry/G80903UK/glyconnect/3641>,
        <http://rdf.glycoinfo.org/glycan/resource-entry/G80903UK/unicarbkb/10>,
        <http://rdf.glycoinfo.org/glycan/resource-entry/null/glytoucan/G80903UK> ;
    glytoucan:has_derivatized_mass <http://rdf.glycoinfo.org/derivatization_type_node/877.31754966> ;
    glytoucan:has_primary_id "G80903UK" .
'''


class LocalClient(glyspace.GlyTouCanRDFClient):
    def __init__(self):
        self.accession_ns = glyspace.NSGlycoinfo
        glyspace.ConjunctiveGraph.__init__(self)
        self.cache = dict()


class RDFClientTest(unittest.TestCase):
    def test_graph_reconstruction(self):
        client = LocalClient()
        client.parse(data=example_n3, format='n3')
        ref = client.get("G80903UK")
        graph = glyspace.ConjunctiveGraph().parse(
            data=ref.to_graph().serialize(format='n3'),
            format='n3')
        self.assertTrue(graph.isomorphic(client))


@skip_not_online
class GlyTouCanRDFClientTest(unittest.TestCase):
    def test_get(self):
        record = glyspace.get("G00034ND")
        self.assertEqual(record.has_primary_id, u'G00034ND')
        self.assertEqual(round(tree(record).mass()), 4104.0)

    def test_from_taxon(self):
        self.assertTrue(len(glyspace.client.from_taxon(9606, limit=10)) == 10)

    def test_with_motif(self):
        self.assertTrue(len(glyspace.structures_with_motif("G00026MO", limit=10)) == 10)

    def test_structure(self):
        self.assertAlmostEqual(glyspace.structure("G00026MO").mass(), 910.327, 2)


if __name__ == '__main__':
    unittest.main()
