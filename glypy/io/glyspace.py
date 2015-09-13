import logging
import re

import requests
from rdflib import ConjunctiveGraph, Namespace, URIRef, Literal
from rdflib.namespace import split_uri

from collections import defaultdict

from glypy.io import glycoct

# http://glytoucan.org/glyspace/documentation/apidoc.html
# http://code.glytoucan.org/system/glyspace

logger = logging.getLogger(__name__)

ssl_verification = False


sparql_endpoint = "http://ts.glytoucan.org/sparql"


def execute_sparql(query):
    payload = {
        "query": query,
        "format": "json",
        "debug": "on",
        "timeout": 10
    }
    r = requests.post(sparql_endpoint, params=payload)
    r.raise_for_status()
    return r.json()


_get_all_glycans_sparql = '''
PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
PREFIX glytoucan:  <http://www.glytoucan.org/glyco/owl/glytoucan#>


SELECT DISTINCT ?Saccharide ?PrimaryId ?Sequence
WHERE {
?Saccharide glytoucan:has_primary_id ?PrimaryId .
?Saccharide glycan:has_glycosequence ?GlycoSequence .
?GlycoSequence glycan:has_sequence ?Sequence .
?GlycoSequence glycan:in_carbohydrate_format glycan:carbohydrate_format_glycoct
}
ORDER BY ?PrimaryId
'''


def download_all_sequences():
    data = execute_sparql(_get_all_glycans_sparql)
    return data


NSGlyTouCan = Namespace("http://www.glytoucan.org/glyco/owl/glytoucan#")
NSGlycan = Namespace("http://purl.jp/bio/12/glyco/glycan#")
NSGlycoinfo = Namespace("http://rdf.glycoinfo.org/glycan/")
NSGlycomeDB = Namespace("http://rdf.glycome-db.org/glycan/")
NSSKOS = Namespace("http://www.w3.org/2004/02/skos/core#")


def _camel_to_snake(name):
    return re.sub(r'([a-z0-9])([A-Z])', r'\1_\2', re.sub(r'(.)([A-Z][a-z]+)', r"\1_\2", name)).lower()


class ReferenceEntity(object):
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

        self.__keys = kwargs.keys()

    def __iter__(self):
        for k in self.__keys:
            yield k, getattr(self, k)

    def __getitem__(self, k):
        return getattr(self, k)

    def __setitem__(self, k, v):
        setattr(self, k, v)

    def __repr__(self):
        tokens = [self.__class__.__name__]
        for k, v in self:
            if isinstance(v, list):
                v_str = ',\n    '.join(map(str, v[:20]))
                if len(v) > 20:
                    v_str += '\n     ...\n'
            else:
                v_str = v
            tokens.append("%s: %s" % (k, v_str))
        return "\n".join(tokens)

    def __root__(self):
        return self.structure_.__root__()

    def __tree__(self):
        return self.structure_


class BoundURIRef(URIRef):
    '''
    A subclass of :class:`rdflib.term.URIRef` which bakes in a way to
    fetch the referenced subgraph (in the semantic web sense) as a :class:`.ReferenceEntity` by keeping
    a reference (in the memory address sense) to the :class:`GlySpaceRDFClient` instance
    which fetched the object it is attached to.

    It has some convenience features for interactive use such as implict
    resource fetching when checking for attribute completions.

    Attributes
    ----------
    _bind_source: GlySpaceRDFClient
    _result_ref: ReferenceEntity
        A reference to the ReferenceEntity fetched from this URI.
        Acts as a cache

    '''
    __slots__ = ("_bind_source", "_result_ref")

    def __new__(cls, value, base=None, source=None):
        rt = URIRef.__new__(cls, value, base)
        rt._bind_source = source
        rt._result_ref = None
        return rt

    def get(self, simplify=True, refresh=False):
        """Get the referenced entity either from :attr:`_bind_source`
        or the cached reference in :attr:`_result_ref`

        Parameters
        ----------
        simplify : bool, optional
            As in :meth:`.GlySpaceRDFClient.get`
        refresh : bool, optional
            If `True`, always request the URI's semantic reference,
            ignoring :attr:`_result_ref`

        Returns
        -------
        ReferenceEntity
        """
        if self._result_ref is None or refresh:
            result = self._bind_source.get(self, simplify=simplify)
            self._result_ref = result
            return result
        else:
            return self._result_ref

    __call__ = get

    def __getattr__(self, name):
        '''
        A convenience method to forward missed attribute lookup
        to the referenced entity. Calls :meth:`__call__`, which
        may initiate a network request.
        '''
        # IPython's displayhook tries to find the attribute _ipython_display_
        # which naturally is not found and initiates a network request. This
        # causes undesirable lag when using the interactive terminal, so short-
        # circuit.
        if name.startswith("_ipython"):
            raise AttributeError()
        reference = self()
        try:
            return getattr(reference, name)
        except AttributeError:
            raise AttributeError("Neither %s nor its reference have an attribute %s" % (self, name))

    def __dir__(self):
        if self._result_ref is None:
            self()
        return dir(self._result_ref)

    def __eq__(self, other):
        result = str(self) == str(other)
        return result


class GlySpaceRDFClient(ConjunctiveGraph):
    predicate_processor_map = {}
    _predicates_seen = set()

    @classmethod
    def register_predicate_processor(cls, predicate):
        def wrapper(f):
            cls.predicate_processor_map[predicate] = f
            return f
        return wrapper

    def __init__(self):
        super(GlySpaceRDFClient, self).__init__(store="SPARQLStore")
        self.open(sparql_endpoint)
        self.bind("glytoucan", NSGlyTouCan)
        self.bind("glycomedb", NSGlycomeDB)
        self.bind("glycan", NSGlycan)
        self.bind("glycoinfo", NSGlycoinfo)

    def accession_to_uriref(self, accession):
        return NSGlycoinfo[accession]

    def get(self, uriref, simplify=True):
        if not isinstance(uriref, URIRef):
            uriref = self.accession_to_uriref(uriref)
        results = defaultdict(list)
        for subject, predicate, obj in set(self.triples((uriref, None, None))):
            predicate_name = _camel_to_snake(split_uri(predicate)[1])
            self._predicates_seen.add(predicate)
            if isinstance(obj, Literal):
                obj = obj.toPython()
            elif isinstance(obj, URIRef):
                obj = BoundURIRef(obj, source=self)
            if predicate in self.predicate_processor_map:
                obj = self.predicate_processor_map[predicate](results, obj)
            if obj is not None:
                results[predicate_name].append(obj)

        # If there were no results, the query might be a predicate, so try to find all the
        # pairs that satisfy it.
        if len(results) == 0:
            predicate_name = _camel_to_snake(split_uri(uriref)[1])
            for subject, predicate, obj in set(self.triples((None, uriref, None))):
                if isinstance(obj, Literal):
                    obj = obj.toPython()
                elif isinstance(obj, URIRef):
                    obj = BoundURIRef(obj, source=self)

                if isinstance(subject, Literal):
                    subject = subject.toPython()
                elif isinstance(subject, URIRef):
                    subject = BoundURIRef(subject, source=self)

                results[predicate_name].append((subject, obj))

        if simplify:
            results = {k: v if len(v) > 1 else v[0] for k, v in results.items()}
        results = ReferenceEntity(**results)
        return results


_client = GlySpaceRDFClient()

get = _client.get


@GlySpaceRDFClient.register_predicate_processor(NSGlycan.has_glycosequence)
def has_glycosequence_processor(data, uri):
    reference = uri()
    if reference.in_carbohydrate_format == NSGlycan.carbohydrate_format_glycoct:
        # trailing underscore in case a URI would claim "structure"
        data["structure_"] = [glycoct.loads(reference.has_sequence)]
    return uri


class ImageResource(dict):
    def __init__(self):
        dict.__init__(self)

    def register_uri(self, uriref):
        reference = uriref()
        symbol_format = split_uri(reference.has_symbol_format)[1].replace("symbol_format_", "")
        url = uriref.toPython()
        self[symbol_format] = url

    def get(self, symbol_format):
        url = self[symbol_format]
        return requests.get(url, verify=ssl_verification)

    def __repr__(self):
        return "ImageResource(%s)" % ', '.join(self.keys())


@GlySpaceRDFClient.register_predicate_processor(NSGlycan.has_image)
def has_image_processor(data, uri):
    if "has_image" in data:
        image_resource_manager = data["has_image"][0]
    else:
        image_resource_manager = data["has_image"] = [ImageResource()]
        image_resource_manager = image_resource_manager[0]
    image_resource_manager.register_uri(uri)
    return None
