import pkg_resources
import logging
from rdflib import Namespace

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

    tree = etree.parse(stream)
    root = tree.getroot()

    entities = []
    for elt in root.iter("{http://www.w3.org/2002/07/owl#}Declaration"):
        decl = elt[0]
        if "IRI" in decl.attrib:
            entities.append(decl.attrib['IRI'].strip("#"))
    return PopulatedNamespace("http://purl.jp/bio/12/glyco/glycan#", entities)


NSGlycan = parse_owl(pkg_resources.resource_stream(__name__, "data/glycan.owl"))


__all__ = [
    "NSGlycan",
]
