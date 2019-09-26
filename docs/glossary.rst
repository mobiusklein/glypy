Glossary of Terms
-----------------

.. glossary::

    Glycan
        A molecule composed of saccharide units and substituents

    Monosaccharide
        A single saccharide unit which is composed of a carbon backbone with hydroxyl
        side chains or substituent groups.

    Topological Matching
        When discussing a *topological* matching traversal, the algorithm does not check
        whether the precise connection locations of edges match, just that there exists
        an edge which connects a matching parent node to a matching child node between two
        trees/graphs. This is particularly useful when unknown linkages are present. For
        example, ``a-D-Neup5Ac-(2-3)-b-D-Galp-(1-3)-b-D-Glcp2NAc`` would match
        ``a-D-Neup5Ac-(2-6)-b-D-Galp-(1-3)-b-D-Glcp2NAc`` despite the difference in the
        linkage between ``a-D-Neup5Ac`` and ``b-D-Galp``. The alternative is a
        :term:`Exact Matching` traversal.

    Exact Matching
        When discussing an *exact* matching traversal, the algorithm requires that the
        positions of all connections of matching nodes between two trees/graphs must also
        match to be considered a full match. For example, ``a-D-Neup5Ac-(2-3)-b-D-Galp-(1-3)-b-D-Glcp2NAc``
        would not match ``a-D-Neup5Ac-(2-6)-b-D-Galp-(1-3)-b-D-Glcp2NAc`` because of the
        difference in the linkage between ``a-D-Neup5Ac`` and ``b-D-Galp``.

    RDF
        Resource Description Framework (RDF) is a web technology standard for specifying
        a semantic description of data. It is used by the glySpace project and affiliated
        databases to represent glycans, glycoconjugates, and related entities.

    SPARQL
        SPARQL is a query language for :term:`rdf`.
