import pkg_resources
import logging
from rdflib import Namespace, Graph, Literal
from rdflib.namespace import RDF, RDFS, FOAF

logging.getLogger("rdflib.term").addHandler(logging.NullHandler())
logging.getLogger("rdflib.term").propagate = False


class PopulatedNamespace(Namespace):
    def __new__(cls, value, members):
        inst = super(cls, cls).__new__(cls, value)
        inst._members = members
        return inst

    def __dir__(self):
        members = dir(super(PopulatedNamespace, self))
        members = sorted(set(list(members) + self._members))
        return members


def parse_owl(stream):
    from lxml import etree

    NSGlycan = None

    def iri(element):
        return NSGlycan[element.attrib['IRI'].strip("#")]

    def abbrev_iri(element):
        return NSGlycan[element.attrib['abbreviatedIRI'].strip("#")]

    tree = etree.parse(stream)
    root = tree.getroot()

    BASE_URI = "http://purl.jp/bio/12/glyco/glycan"

    entities = []
    for elt in root.iter("{http://www.w3.org/2002/07/owl#}Declaration"):
        decl = elt[0]
        if "IRI" in decl.attrib:
            entities.append(decl.attrib['IRI'].strip("#"))
    NSGlycan = PopulatedNamespace(BASE_URI + "#", entities)

    GlycoRDF = Graph()
    for elt in root.iter("{http://www.w3.org/2002/07/owl#}ClassAssertion"):
        cls = iri(elt.find("{http://www.w3.org/2002/07/owl#}Class"))
        inst = iri(elt.find("{http://www.w3.org/2002/07/owl#}NamedIndividual"))
        GlycoRDF.add((inst, RDF.type, cls))

    for elt in root.iter("{http://www.w3.org/2002/07/owl#}ObjectPropertyAssertion"):
        prop = iri(elt.find("{http://www.w3.org/2002/07/owl#}ObjectProperty"))
        subj, val = map(iri, elt.findall(
            "{http://www.w3.org/2002/07/owl#}NamedIndividual"))
        GlycoRDF.add((subj, prop, val))

    for elt in root.iter("{http://www.w3.org/2002/07/owl#}DataPropertyAssertion"):
        prop = iri(elt.find("{http://www.w3.org/2002/07/owl#}DataProperty"))
        subj = iri(elt.find("{http://www.w3.org/2002/07/owl#}NamedIndividual"))
        val = elt.find("{http://www.w3.org/2002/07/owl#}Literal").text
        GlycoRDF.add((subj, prop, Literal(val)))

    for elt in root.iter("{http://www.w3.org/2002/07/owl#}SubClassOf"):
        try:
            cls1, cls2 = map(iri, elt.findall(
                "{http://www.w3.org/2002/07/owl#}Class"))
        except (KeyError, ValueError):
            continue
        GlycoRDF.add((cls1, RDFS.subClassOf, cls2))

    for elt in root.iter("{http://www.w3.org/2002/07/owl#}AnnotationAssertion"):
        annotation = elt.find(
            "{http://www.w3.org/2002/07/owl#}AnnotationProperty").attrib['abbreviatedIRI']
        ns, tp = annotation.split(":")
        if ns == "rdf":
            rel = RDF[tp]
        elif ns == 'rdfs':
            rel = RDFS[tp]
        elif ns == 'foaf':
            rel = FOAF[tp]
        else:
            raise ValueError(annotation)
        try:
            subj = NSGlycan[elt.find("{http://www.w3.org/2002/07/owl#}IRI").text.strip("#")]
        except AttributeError:
            continue
        val = elt.find("{http://www.w3.org/2002/07/owl#}Literal").text
        GlycoRDF.add((subj, rel, Literal(val)))

    return NSGlycan, GlycoRDF


NSGlycan, glycordf = parse_owl(
    pkg_resources.resource_stream(__name__, "data/glycan.owl"))


__all__ = [
    "NSGlycan",
    "glycordf",
]
