from rdflib import Graph

g = Graph()
g.bind('glycan', "http://purl.jp/bio/12/glyco/glycan#")
g.parse("https://raw.githubusercontent.com/ReneRanzinger/GlycoRDF/master/ontology/glycan.owl", format="turtle")
g.serialize(open("glyspace.ttl", "wb"), 'turtle')
