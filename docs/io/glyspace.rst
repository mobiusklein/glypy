glySpace Clients
================

.. currentmodule:: glypy.io.glyspace

The :mod:`glyspace` module provides an API for interacting with `GlyTouCan <https://glytoucan.org/>`_
and `UnicarbKB <http://www.unicarbkb.org/>`_. The interaction may be done by executing prepared
:term:`SPARQL` queries through the provided methods, executing user-provided :term:`SPARQL` queries,
:term:`RDF` Graph-operation supported by `RDFLib <https://rdflib.readthedocs.io/en/stable/>`_,
or using :term:`RDF`-object mapping via :meth:`~.RDFClientBase.get` and related methods.

The module contains a pre-created instance of :class:`GlyTouCanRDFClient` named `glytoucan_client` whose
:meth:`GlyTouCanRDFClient.get`, :meth:`GlyTouCanRDFClient.structure`, :meth:`GlyTouCanRDFClient.from_taxon`,
and :meth:`GlyTouCanRDFClient.structures_with_motif` methods are available as top-level functions of the
module. There is also a :class:`UnicarbKBRDFClient` object called `unikarbkb_client` ready to use.


RDF Namespaces that are pre-created
-----------------------------------
.. code-block:: python

    NSGlyTouCan = Namespace("http://www.glytoucan.org/glyco/owl/glytoucan#")
    NSGlycan = Namespace("http://purl.jp/bio/12/glyco/glycan#")
    NSGlycoinfo = Namespace("http://rdf.glycoinfo.org/glycan/")
    NSGlycomeDB = Namespace("http://rdf.glycome-db.org/glycan/")
    NSSKOS = Namespace("http://www.w3.org/2004/02/skos/core#")
    NSUniprotCore = Namespace("http://purl.uniprot.org/core/")
    NSUniprotEntity = Namespace("http://purl.uniprot.org/uniprot/")
    NSTaxonomy = Namespace("http://purl.uniprot.org/taxonomy/")


RDF Graph Clients
~~~~~~~~~~~~~~~~~

.. autoclass:: glypy.io.glyspace.RDFClientBase

.. autoclass:: glypy.io.glyspace.GlyTouCanRDFClient
    :members: get, query, __getitem__, triples, from_taxon, structures_with_motif

.. autoclass:: glypy.io.glyspace.UnicarbKBRDFClient


RDF-Object Mapping Components
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: glypy.io.glyspace.BoundURIRef
    :members: get, __call__, __getattr__

.. autoclass:: glypy.io.glyspace.ReferenceEntity
    :members:


