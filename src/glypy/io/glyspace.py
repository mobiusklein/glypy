'''The :mod:`glyspace` module provides an API for interacting with `GlyTouCan <https://glytoucan.org/>`_
and `UnicarbKB <http://www.unicarbkb.org/>`_. The interaction may be done by executing prepared
:term:`SPARQL` queries through the provided methods, executing user-provided :term:`SPARQL` queries,
:term:`RDF` Graph-operation supported by `RDFLib <https://rdflib.readthedocs.io/en/stable/>`_,
or using :term:`RDF`-object mapping via :meth:`GlySpaceRDFClient.get` and related methods.

'''

import logging
import re
import warnings

from numbers import Number
from collections import defaultdict, OrderedDict
from six import text_type

try:
    from urllib import quote
except ImportError:
    from urllib.parse import quote

import requests
with warnings.catch_warnings():
    from rdflib import ConjunctiveGraph, Namespace, URIRef, Literal, Graph
    from rdflib.namespace import split_uri
    from glypy.io import glycoct, iupac, wurcs, _glycordf

# http://glytoucan.org/glyspace/documentation/apidoc.html
# http://code.glytoucan.org/system/glyspace

logger = logging.getLogger(__name__)

ssl_verification = False


sparql_endpoint = "http://ts.glytoucan.org/sparql"

_get_all_monosaccharides_sparql = r'''
PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
PREFIX glytoucan:  <http://www.glytoucan.org/glyco/owl/glytoucan#>

SELECT DISTINCT ?Saccharide ?Sequence
WHERE {
    ?Saccharide glycan:has_glycosequence ?GlycoSequence .
    ?GlycoSequence glycan:has_sequence ?Sequence .
    ?GlycoSequence glycan:in_carbohydrate_format glycan:carbohydrate_format_glycoct .
    ?Saccharide a glycan:monosaccharide
}
'''


_get_all_glycans_sparql = r'''
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

NSGlycan = _glycordf.NSGlycan
NSGlyTouCan = Namespace("http://www.glytoucan.org/glyco/owl/glytoucan#")
NSGlycoinfo = Namespace("http://rdf.glycoinfo.org/glycan/")
NSGlycomeDB = Namespace("http://rdf.glycome-db.org/glycan/")
NSSKOS = Namespace("http://www.w3.org/2004/02/skos/core#")
NSUniprotCore = Namespace("http://purl.uniprot.org/core/")
NSUniprotEntity = Namespace("http://purl.uniprot.org/uniprot/")
NSTaxonomy = Namespace("http://purl.uniprot.org/taxonomy/")
NSDCTerms = Namespace("http://purl.org/dc/terms/")
NSGlycoconjugate = Namespace("http://purl.jp/bio/12/glyco/conjugate#")
NSUnicarbKB = Namespace("http://rdf.unicarbkb.org/")
NSUnicarbKBStructure = Namespace("http://rdf.unicarbkb.org/structure/")
NSUnicarbKBProtein = Namespace("http://unicarbkb.org/unicarbkbprotein/")
NSFaldo = Namespace("http://www.biohackathon.org/resource/faldo/")
NSPato = Namespace("http://purl.obolibrary.org/obo/uo.owl")
NSUnicorn = Namespace("http://purl.jp/bio/12/glyco/unicorn/")
NSOwl = Namespace("http://www.w3.org/2002/07/owl#")
NSSIO = Namespace("http://semanticscience.org/resource/")
NSFOAF = Namespace("http://xmlns.com/foaf/0.1/")
NSRDF = Namespace("http://www.w3.org/1999/02/22-rdf-syntax-ns#")

GlycoRDF = _glycordf.glycordf

_uniprot_taxon_uri = u"http://www.uniprot.org/taxonomy/{}.rdf"


def _camel_to_snake(name):
    '''
    Many properties and predicate names are in camelCase, which is
    not consistent with Python naming conventions.
    '''
    return re.sub(r'([a-z0-9])([A-Z])', r'\1_\2', re.sub(r'(.)([A-Z][a-z]+)', r"\1_\2", name)).lower()


def _uri_to_identifier(uri):
    name = _camel_to_snake(split_uri(uri)[1])
    return name


class LRUDict(object):  # pragma: no cover
    def __init__(self, *args, **kwargs):
        maxsize = kwargs.pop("maxsize", 24)
        self.store = OrderedDict()
        self.maxsize = maxsize
        self.purge()

    def __len__(self):
        return len(self.store)

    def popitem(self, last=True):
        return self.store.popitem(last=last)

    def pop(self, key, default=None):
        return self.store.pop(key, default)

    def purge(self):
        overflow = max(0, len(self) - self.maxsize)
        for _ in range(overflow):
            self.popitem(last=False)

    def __repr__(self):
        return "LRUDict(%r)" % (dict(self.store),)

    def __contains__(self, key):
        return key in self.store

    def __iter__(self):
        return iter(self.store)

    def keys(self):
        return self.store.keys()

    def values(self):
        return self.store.values()

    def items(self):
        return self.store.items()

    def __getitem__(self, key):
        value = self.store[key]
        self._mark_used(key)
        return value

    def __setitem__(self, key, value):
        self.store[key] = value
        self.purge()

    def _mark_used(self, key):
        value = self.store.pop(key, None)
        self.store[key] = value


class PredicateDescriptor(text_type):

    """A specialization of the unicode text type for representing a string which
    is a convenient suffix of a complete URI denoting a predicate.

    This type is used to attach a more verbose URI to a short string used for a
    :class:`ReferenceEntity`'s attributes.

    This class should be immutable, and caches all unique URIs.
    """

    __slots__ = ("source", )
    _cache = {}

    @classmethod
    def bind(cls, uri):
        return cls(_uri_to_identifier(uri), uri)

    def __new__(cls, value, source):
        if source in cls._cache:
            return cls._cache[source]
        named = super(PredicateDescriptor, cls).__new__(cls, value)
        named.source = source
        cls._cache[source] = named
        return named


class ReferenceEntity(object):
    '''
    A ReferenceEntity is a generic type to for storing results from
    a semantic query. Its attribute names are usually derived from
    predicates, and their values may be scalar or lists. It is usually
    constructed from :meth:`get` from a :class:`BoundURIRef` instance.

    A preprocessor may add new attributes to a ReferenceEntity during
    construction, such as `structure_` with a |Glycan| instance value
    when the referenced entity satisfies `glycan:has_glycosequence`.
    '''

    query_type = 'subject'

    def __init__(self, uriref, **kwargs):
        self._uriref = uriref
        for k, v in kwargs.items():
            setattr(self, k, v)

        self._keys = {k for k in kwargs.keys() if isinstance(k, PredicateDescriptor)}
        self._object_of = None

    def __iter__(self):
        for k in self._keys:
            yield k, getattr(self, k)

    def __getitem__(self, k):  # pragma: no cover
        return getattr(self, k)

    def __setitem__(self, k, v):  # pragma: no cover
        setattr(self, k, v)
        if isinstance(k, PredicateDescriptor):
            self._keys.add(k)

    def __repr__(self):  # pragma: no cover
        tokens = [self.__class__.__name__]
        for k, v in self:
            if isinstance(v, list):
                v_str = ',\n    '.join(map(str, v[:20]))
                if len(v) > 20:
                    v_str += '\n    ...\n'
            else:
                v_str = v
            tokens.append("%s: %s" % (k, v_str))
        return "\n".join(tokens)

    def __root__(self):
        return self.structure_.__root__()

    def __tree__(self):
        return self.structure_

    @property
    def _graph(self):
        return self._uriref._graph

    @property
    def object_of(self):
        if self._object_of is None:
            self._object_of = self._graph.get(self._uriref, query_type='object')
        return self._object_of

    @property
    def uriref(self):
        return self._uriref

    def to_graph(self, graph=None, deep=False, visited=None):
        if graph is None:
            graph = Graph()
            for prefix, ns in self._graph.namespaces():
                graph.bind(prefix, ns)

        if visited is None:
            visited = set()
        key = (self.uriref, self.query_type)
        if key in visited:
            return graph
        else:
            visited.add(key)

        for predicate, value in self:
            if isinstance(value, (text_type, Number)):
                is_uri = isinstance(value, URIRef)
                if not is_uri:
                    value_n = Literal(value)
                else:
                    value_n = URIRef(value)
                graph.add((self.uriref, URIRef(predicate.source), value_n))
                if is_uri and deep:
                    try:
                        value.to_graph(graph, deep=deep, visited=visited)
                    except AttributeError:
                        self._graph.get(value).to_graph(graph, deep=deep, visited=visited)
            elif isinstance(value, list):
                value_collection = value
                for value in value_collection:
                    if isinstance(value, (text_type, Number)):
                        is_uri = isinstance(value, URIRef)
                        if not is_uri:
                            value_n = Literal(value)
                        else:
                            value_n = URIRef(value)
                        graph.add((self.uriref, URIRef(predicate.source), value_n))
                        if is_uri and deep:
                            try:
                                value.to_graph(graph, deep=deep, visited=visited)
                            except AttributeError:
                                self._graph.get(value).to_graph(graph, deep=deep, visited=visited)
                    else:
                        raise TypeError(
                            "Could not translate %r of type %r for predicate %s" % (
                                value, type(value), predicate.source))
            else:
                raise TypeError(
                    "Could not translate %r of type %r for predicate %s" % (
                        value, type(value), predicate.source))

        if deep:
            self.object_of.to_graph(graph, deep=deep, visited=visited)
        return graph


class PredicateCollection(ReferenceEntity):
    query_type = 'predicate'

    def to_graph(self, graph=None, deep=False, visited=None):
        if graph is None:
            graph = Graph()
            for prefix, ns in self._graph.namespaces():
                graph.bind(prefix, ns)
        if visited is None:
            visited = set()
        key = (self.uriref, self.query_type)
        if key in visited:
            return graph
        else:
            visited.add(key)
        for subj, obj in self:
            if isinstance(subj, (unicode, str, int, float)):
                is_uri_subj = isinstance(subj, URIRef)
                if not is_uri_subj:
                    subj = Literal(subj)
            if isinstance(obj, (unicode, str, int, float)):
                is_uri_obj = isinstance(obj, URIRef)
                if not is_uri_obj:
                    obj = Literal(obj)
            graph.add((subj, self.uriref, obj))
            if is_uri_subj and deep:
                try:
                    subj.to_graph(graph, deep=deep, visited=visited)
                except AttributeError:
                    self._graph.get(subj).to_graph(graph, deep=deep, visited=visited)
            if is_uri_obj:
                try:
                    obj.to_graph(graph, deep=deep, visited=visited)
                except AttributeError:
                    self._graph.get(obj).to_graph(graph, deep=deep, visited=visited)

        return graph


class ObjectOfEntity(ReferenceEntity):
    query_type = 'object'

    def to_graph(self, graph=None, deep=False, visited=None):
        if graph is None:
            graph = Graph()
            for prefix, ns in self._graph.namespaces():
                graph.bind(prefix, ns)

        if visited is None:
            visited = set()
        key = (self.uriref, self.query_type)
        if key in visited:
            return graph
        else:
            visited.add(key)

        for predicate, subj in self:
            if isinstance(subj, (unicode, str, int, float)):
                is_uri = isinstance(subj, URIRef)
                if not is_uri:
                    subj = Literal(subj)
                graph.add((subj, URIRef(predicate.source), self.uriref))
                if is_uri and deep:
                    try:
                        subj.to_graph(graph, deep=deep, visited=visited)
                    except AttributeError:
                        self._graph.get(subj).to_graph(graph, deep=deep, visited=visited)
            elif isinstance(subj, list):
                subj_collection = subj
                for subj in subj_collection:
                    if isinstance(subj, (unicode, str, int, float)):
                        is_uri = isinstance(subj, URIRef)
                        if not is_uri:
                            subj = Literal(subj)
                        graph.add((subj, URIRef(predicate.source), self.uriref))
                        if is_uri and deep:
                            try:
                                subj.to_graph(graph, deep=deep, visited=visited)
                            except AttributeError:
                                self._graph.get(subj).to_graph(graph, deep=deep, visited=visited)
                    else:
                        raise TypeError(
                            "Could not translate %r of type %r for predicate %s" % (
                                subj, type(subj), predicate.source))
            else:
                raise TypeError(
                    "Could not translate %r of type %r for predicate %s" % (
                        subj, type(subj), predicate.source))
        return graph


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
    _graph: :class:`RDFClientBase`
    _result_ref: :class:`ReferenceEntity`
        A reference to the ReferenceEntity fetched from this URI.
        Acts as a cache

    '''
    __slots__ = ("_graph", "_result_ref")

    def __new__(cls, value, base=None, source=None):
        rt = URIRef.__new__(cls, value, base)
        rt._graph = source
        rt._result_ref = None
        return rt

    def __hash__(self):
        return URIRef.__hash__(self)

    def get(self, simplify=True, refresh=False, query_type='auto'):
        """Get the referenced entity either from :attr:`_graph`
        or the cached reference in :attr:`_result_ref`

        Parameters
        ----------
        simplify : bool, optional
            As in :meth:`RDFClientBase.get`
        refresh : bool, optional
            If `True`, always request the URI's semantic reference,
            ignoring :attr:`_result_ref`

        Returns
        -------
        ReferenceEntity
        """
        if self._result_ref is None or refresh:
            result = self._graph.get(self, simplify=simplify, query_type=query_type)
            self._result_ref = result
            return result
        else:
            return self._result_ref

    __call__ = get

    def _repr_pretty_(self, p, cycle):
        return p.pretty(str(self))

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

    def __dir__(self):  # pragma: no cover
        if self._result_ref is None:
            self()
        return dir(self._result_ref)

    def __eq__(self, other):
        '''
        Overrides the equality method of :class:`rdflib.term.URIRef` which
        does exact :func:`type` comparison before comparing contents to short-circuit
        on non-URIs.
        '''
        result = str(self) == str(other)
        return result

    def n3(self, namespace_manager=None):
        try:
            text = super(BoundURIRef, self).n3(namespace_manager)
            return text
        except Exception:
            if namespace_manager:
                return namespace_manager.normalizeUri(self)
            else:
                return "<%s>" % self


class ChainFunctionDict(defaultdict):
    '''
    A wrapper around :class:`defaultdict(list)` which
    keys on :class:`rdflib.term.URIRef` strings. Added values
    should be callables which will be invoked in the order given on
    an entity state `dict` and each `object` which is linked by the key `predicate`.
    '''
    def __init__(self, *args, **kwargs):
        defaultdict.__init__(self, list, *args, **kwargs)

    def __call__(self, predicate, state, obj):
        ret_val_was_none = False
        non_none_ret_val = None
        for f in self[predicate]:
            ret_val = f(state, obj)
            if ret_val is not None:
                non_none_ret_val = ret_val
            else:
                ret_val_was_none = True
        if ret_val_was_none:
            return None
        else:
            return non_none_ret_val


class RDFClientBase(ConjunctiveGraph):
    _predicates_seen = set()
    predicate_processor_map = ChainFunctionDict()

    @classmethod
    def register_predicate_processor(cls, predicate):
        """Decorator to register a callable processor for a `URIRef` `predicate` with this type's
        :attr:`predicate_processor_map`. The actual decorated callable is returned unchanged.

        Parameters
        ----------
        predicate : rdflib.term.URIRef or str
            The type of URI to add the decorated callable
            to the processor chain of

        Returns
        -------
        callable
        """
        def wrapper(f):
            cls.predicate_processor_map[predicate].append(f)
            return f
        return wrapper

    def __init__(self, sparql_endpoint, accession_ns, cache_size=100):
        super(RDFClientBase, self).__init__(store="SPARQLStore")
        self.open(sparql_endpoint)
        self.accession_ns = accession_ns
        self.cache = LRUDict(maxsize=cache_size)

    def accession_to_uriref(self, accession):
        """Utility method to translate free strings into full URIs
        derived from this instance's :attr:`accession_ns`

        Parameters
        ----------
        accession : str
            A regular string comprised of just the accession number of
            an entity.

        Returns
        -------
        rdflib.term.URIRef
        """
        return self.accession_ns[str(accession)]

    def get(self, uriref, simplify=True, query_type='auto'):
        """Download all related information for `uriref` from the remote
        data source.

        Collects all the triples from the remote data source where `uriref` is
        the subject. If `uriref` is not the subject of any triples, it is re-queried
        as a predicate, storing the subject-object pairs.

        Any objects (and subjects) which are themselves :class:`rdflib.term.URIRef` instances
        will be converted into :class:`BoundURIRef` which will silently fetch the relevant
        entity from the remote source.

        If the predicate matches a processor rules, instead of it's object value being
        stored, the object will be transformed by each rule in the processor chain.

        Parameters
        ----------
        uriref: str or rdflib.term.URIRef
            A subject, predicate, or a database accession number to transform through
            :meth:`accession_to_uriref`
        simplify: bool, optional
            If true, any predicate with a single value will be a scalar,
            and any other will be a list.

        Returns
        -------
        ReferenceEntity
            An object representing the subject whose attributes are named after
            predicates with their objects as values.
        """
        if not isinstance(uriref, URIRef):
            uriref = self.accession_to_uriref(uriref)
        if (uriref, query_type) in self.cache:
            return self.cache[uriref, query_type]
        results = defaultdict(list)
        if query_type == 'auto' or query_type == 'subject':
            for subject, predicate, obj in set(self.triples((uriref, None, None))):
                predicate_name = PredicateDescriptor.bind(predicate)
                self._predicates_seen.add(predicate)
                if isinstance(obj, Literal):
                    obj = obj.toPython()
                elif isinstance(obj, URIRef):
                    obj = BoundURIRef(obj, source=self)
                if predicate in self.predicate_processor_map:
                    obj = self.predicate_processor_map(predicate, results, obj)
                if obj is not None:
                    results[predicate_name].append(obj)
            if results:
                query_type = 'subject'

        # If there were no results, the query might be a predicate, so try to find all the
        # pairs that satisfy it.
        if (len(results) == 0 and query_type == 'auto') or query_type == 'predicate':
            try:
                predicate_name = PredicateDescriptor.bind(uriref)
            except Exception:
                predicate_name = None
            if predicate_name is not None:
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
            if results:
                query_type = 'predicate'

        if query_type == 'object':
            for subject, predicate, obj in set(self.triples((None, None, uriref))):
                predicate_name = PredicateDescriptor.bind(predicate)
                self._predicates_seen.add(predicate)
                if isinstance(subject, Literal):
                    subject = subject.toPython()
                elif isinstance(subject, URIRef):
                    subject = BoundURIRef(subject, source=self)
                if predicate in self.predicate_processor_map:
                    subject = self.predicate_processor_map(predicate, results, subject)
                if subject is not None:
                    results[predicate_name].append(subject)

        if simplify:
            results = {k: v if len(v) > 1 else v[0] for k, v in results.items()}
        if query_type in ('subject', 'auto'):
            entity_type = ReferenceEntity
        elif query_type == 'predicate':
            entity_type = PredicateCollection
        elif query_type == 'object':
            entity_type = ObjectOfEntity
        else:
            raise ValueError("Could not determine query type %s" % (query_type,))
        results = entity_type(BoundURIRef(uriref, source=self), **results)
        self.cache[uriref, query_type] = results
        return results


class UniprotRDFClient(RDFClientBase):
    predicate_processor_map = ChainFunctionDict()
    _predicates_seen = set()
    _sparql_endpoint_uri = "http://sparql.uniprot.org/sparql/"

    def __init__(self):
        super(UniprotRDFClient, self).__init__(self._sparql_endpoint_uri, NSUniprotEntity)
        self.bind("up", NSUniprotCore)
        self.bind("taxonomy", NSTaxonomy)


class GlyTouCanRDFClient(RDFClientBase):
    '''
    An RDF Client for glySpace. The default namespace is `glycoinfo`, and
    the following namespaces are bound:

    .. code-block:: python

        glytoucan = NSGlyTouCan
        glycomedb = NSGlycomeDB
        glycan = NSGlycan
        glycoinfo = NSGlycoinfo
        skos = NSSKOS


    Attributes
    ----------
    predicate_processor_map: :class:`ChainFunctionDict`
        A dictionary keeping track of the chain of processors registered for
        each predicate.
    _sparql_endpoint_uri: str
        The web address to use as the remote backend for SPARQL queries.
        Passed on to :class:`rdflib.ConjunctiveGraph` and the `SPARQLStore`
        storage plugin.

    '''

    predicate_processor_map = ChainFunctionDict()
    _predicates_seen = set()
    _sparql_endpoint_uri = sparql_endpoint

    def __init__(self):
        super(GlyTouCanRDFClient, self).__init__(self._sparql_endpoint_uri, NSGlycoinfo)
        self.bind("glytoucan", NSGlyTouCan)
        self.bind("glycomedb", NSGlycomeDB)
        self.bind("glycan", NSGlycan)
        self.bind("glycoinfo", NSGlycoinfo)
        self.bind("skos", NSSKOS)
        self.bind("dcterms", NSDCTerms)

    def from_taxon(self, taxon, limit=None):
        r"""Fetch all accession numbers for all structures
        from the given taxonomic identifier, up to `limit` records.

        Equivalent to the following SPARQL

        .. code-block:: sparql

            SELECT DISTINCT ?saccharide WHERE {
                ?saccharide a glycan:saccharide .
                ?saccharide skos:exactMatch ?gdb .
                ?gdb glycan:has_reference ?ref .
                ?ref glycan:is_from_source ?source .
                ?source glycan:has_taxon ?taxon
                FILTER REGEX(str(?taxon), "http://www.uniprot.org/taxonomy/<taxonomy-accession>.rdf")
            }

        The REGEX filter is used because at current taxonomic information in Glycome-DB
        is encoded as a string instead of a URI.

        Parameters
        ----------
        taxon : str or int
            A string or number which corresponds to the taxonomy database
            id for the taxon of interest.
        limit : int, optional
            The maximum number of results to retrieve.

        Returns
        -------
        list of BoundURIRef
        """
        sparql = r'''
        SELECT DISTINCT ?accession ?taxon_id
        WHERE{
            ?saccharide  glytoucan:has_primary_id ?accession .
            ?saccharide glycan:is_from_source ?source.
            ?source a glycan:Source .
            VALUES ?taxon_id {"%s"}
            ?source dcterms:identifier ?taxon_id
        }
        '''
        query_string = sparql % (str(taxon), )
        if limit is not None:
            query_string += " limit %d" % limit
        results = self.query(query_string)
        k = results.vars[0]
        return [BoundURIRef(row[k], source=self) for row in results.bindings]

    def structures_with_motif(self, motif, limit=None):
        r"""Fetch all accession numbers and structures for all structures
        which contain the given motif accession, up to `limit` records.

        Equivalent to the following SPARQL

        .. code-block:: sparql

            SELECT DISTINCT ?saccharide WHERE {
                ?saccharide a glycan:saccharide .
                ?saccharide glycan:has_motif <motif-accession>
            }

        Parameters
        ----------
        motif : str
            The accession number of the motif of interest.
        limit : int, optional
            The maximum number of results to retrieve.

        Returns
        -------
        list of :class:`ReferenceEntity`
        """
        sparql = r'''
        SELECT DISTINCT ?saccharide WHERE {
                ?saccharide glycan:has_motif %s
        }
        '''
        if isinstance(motif, URIRef):
            motif_str = self.qname(motif)
        else:
            motif_str = self.qname(NSGlycoinfo[motif])
        query_string = sparql % motif_str
        if limit is not None:
            query_string += " limit %d" % limit
        results = self.query(query_string)
        k = results.vars[0]
        return [ReferenceEntity(row[k]) for row in results.bindings]

    def structure(self, *accessions):
        accumulator = []
        sparql = r'''
        SELECT DISTINCT ?saccharide ?glycoct WHERE {
            ?saccharide glycan:has_glycosequence ?sequence .
            FILTER CONTAINS(str(?sequence), "glycoct") .
            ?sequence glycan:has_sequence ?glycoct .
            FILTER ("%s" = str(?saccharide))
        }
        '''
        for accession in accessions:
            if isinstance(accession, URIRef):
                accession_str = str(accession)
            else:
                accession_str = str(NSGlycoinfo[accession])
            query_string = sparql % accession_str
            results = self.query(query_string)
            g = results.vars[1]
            glycoct_string = results.bindings[0][g]
            structure = glycoct.loads(glycoct_string)
            accumulator.append(structure)
        if len(accumulator) == 1:
            return accumulator[0]
        else:
            return accumulator


class UnicarbKBRDFClient(RDFClientBase):
    predicate_processor_map = ChainFunctionDict()
    _predicates_seen = set()
    _sparql_endpoint_uri = 'http://203.101.226.16:40935/unicarbkbv2.0.1/query'

    def __init__(self):
        super(UnicarbKBRDFClient, self).__init__(self._sparql_endpoint_uri, NSUnicarbKBStructure)
        self.bind("glytoucan", NSGlyTouCan)
        self.bind("glycomedb", NSGlycomeDB)
        self.bind("glycan", NSGlycan)
        self.bind("glycoinfo", NSGlycoinfo)
        self.bind("skos", NSSKOS)
        self.bind("uniprot", NSUniprotCore)
        self.bind("gco", NSGlycoconjugate)
        self.bind('ukp', NSUnicarbKBProtein)
        self.bind("uk", NSUnicarbKB)
        self.bind("sio", NSSIO)
        self.bind("unicorn", NSUnicorn)
        self.bind("owl", NSOwl)
        self.bind("dcterms", NSDCTerms)
        self.bind("faldo", NSFaldo)


GlySpaceRDFClient = GlyTouCanRDFClient


glytoucan_client = client = GlyTouCanRDFClient()
unikarbkb_client = UnicarbKBRDFClient()


get = client.get
query = client.query
structure = client.structure
from_taxon = client.from_taxon
structures_with_motif = client.structures_with_motif
triples = client.triples


@GlyTouCanRDFClient.register_predicate_processor(NSGlycan.has_glycosequence)
@UnicarbKBRDFClient.register_predicate_processor(NSGlycan.has_glycosequence)
def has_glycosequence_processor(state, uri):
    """Detect and extract GlycoCT sequence data and parse
    into a |Glycan| object.

    Parameters
    ----------
    state : ReferenceEntity or dict
        The key-value store to add annotation to.
    uri : rdflib.term.URIRef
        The `URIRef` to load structure data from.

    Returns
    -------
    BoundURIRef
    """
    reference = uri(query_type='subject')
    # print(uri, reference, state)
    try:
        in_format = (reference.in_carbohydrate_format)
    except AttributeError:
        return uri
    if in_format == (NSGlycan.carbohydrate_format_glycoct):
        # trailing underscore in case a URI would claim "structure"
        state["structure_"] = [glycoct.loads(reference.has_sequence)]
    if "structure_" not in state and in_format == (NSGlycan.carbohydrate_format_wurcs):
        try:
            state["structure_"] = [wurcs.loads(reference.has_sequence)]
        except (wurcs.WURCSError, KeyError):
            pass
    if "structure_" not in state and in_format == (NSGlycan.carbohydrate_format_iupac_extended):
        try:
            state["structure_"] = [iupac.loads(reference.has_sequence)]
        except (iupac.IUPACError, KeyError):
            pass
    return uri


class ImageResource(dict):
    """ImageResource wraps image resource URIs provided by the data source,
    unifying the image fetching API.

    ImageResource inherits from |dict|, and its keys
    are pairs of (symbol_format, image_format), e.g. ("cgf", "png").

    The wrapping :meth:`get` method will fetch the relevant URI from
    storage and fetches the image data, which can be written to file.

    Attributes
    ----------
    _uriref_list: list
        List of wrapped resource URIs.
    """
    def __init__(self):
        dict.__init__(self)
        self._uriref_list = []

    def _materialize(self):
        for uriref in self._uriref_list:
            reference = uriref()
            image_format = reference.format.split("/")[1]
            symbol_format = split_uri(reference.has_symbol_format)[1].replace("symbol_format_", "")
            url = uriref.toPython()
            self[symbol_format, image_format] = url

    def register_uri(self, uriref):
        self._uriref_list.append(uriref)
        reference = uriref()
        image_format = reference.format.split("/")[1]
        symbol_format = split_uri(reference.has_symbol_format)[1].replace("symbol_format_", "")
        url = uriref.toPython()
        self[symbol_format, image_format] = url

    def get(self, symbol_format, image_format='png'):
        """Fetch image data from a remote source.

        Parameters
        ----------
        symbol_format : str
            The symbol format of the image to fetch.
        image_format : str, optional
            The image format to request. Defaults to 'png'. The supported
            formats will vary depending upon the remote source.

        Returns
        -------
        name : TYPE
            Description
        """
        url = self[symbol_format, image_format]
        return requests.get(url, verify=ssl_verification)

    def __repr__(self):
        return "ImageResource(%s)" % ', '.join(map(str, self.keys()))


# @GlyTouCanRDFClient.register_predicate_processor(NSGlycan.has_image)
def has_image_processor(state, uri):
    '''
    Intercept raw URIRefs to image generation services and store them in
    the :class:`ImageResource` instance associated with `state`. Since the
    predicate name is being used by the ImageResource, return `None` so
    the URIRef doesn't try to overwrite it.
    '''
    if "has_image" in state:
        image_resource_manager = state["has_image"][0]
    else:
        image_resource_manager = state["has_image"] = [ImageResource()]
        image_resource_manager = image_resource_manager[0]
    image_resource_manager.register_uri(uri)
    return None
