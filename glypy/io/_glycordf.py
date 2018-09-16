import pkg_resources
import logging
from rdflib import Graph, Namespace

logging.getLogger("rdflib.term").addHandler(logging.NullHandler())
logging.getLogger("rdflib.term").propagate = False


NSGlycan = ("http://purl.jp/bio/12/glyco/glycan#")
glycordf = Graph().parse(pkg_resources.resource_stream(__name__, "data/glycordf.ttl"), format='turtle')
_entities = set(glycordf.subjects()) | set(glycordf.objects())
_entities = [
    entity.replace(NSGlycan, "") for entity in _entities if (NSGlycan in entity)
]


class PopulatedNamespace(Namespace):
    def __new__(cls, value, members):
        inst = super(cls, cls).__new__(cls, value)
        inst._members = members
        return inst

    def __dir__(self):
        members = dir(super(PopulatedNamespace, self))
        members = sorted(set(list(members) + self._members))
        return members


NSGlycan = PopulatedNamespace(NSGlycan, _entities)


__all__ = [
    "NSGlycan",
    "glycordf",
]
