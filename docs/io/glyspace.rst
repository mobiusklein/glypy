===============
glySpace Client
===============

.. currentmodule:: glypy.io.glyspace

This module implements a client for communicating with remote data stores part of the **glySpace** project.
Currently, this communicates with the *SPARQL* endpoint hosted by https://glytoucan.org/.


The module contains a pre-created instance of :class:`GlySpaceRDFClient` named `client` whose
:meth:`get`, :meth:`structure`, :meth:`from_taxon`, and :meth:`structures_with_motif` methods
are available as top-level functions of the module.


-----------------------------------
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


.. autoclass:: glypy.io.glyspace.GlySpaceRDFClient
    :members: get, query, __getitem__, triples, from_taxon, structures_with_motif


.. autoclass:: glypy.io.glyspace.BoundURIRef
    :members: get, __call__, __getattr__

.. autoclass:: glypy.io.glyspace.ReferenceEntity
    :members:


