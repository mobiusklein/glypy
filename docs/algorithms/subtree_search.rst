Substructure Search Methods
===========================

.. automodule:: glypy.algorithms.subtree_search

.. contents::
    :local:

Inclusion Comparison Algorithms
-------------------------------
Inclusion-based search asks whether one graph is completely contained in another sub-graph,
using either :term:`Exact Matching` or :term:`Topological Matching` according to the chosen
approach. These comparisons can be made fuzzy, using :func:`~.similarity.commutative_similarity`
to evaluate matches between individual nodes in a graph.

Direct Comparators
^^^^^^^^^^^^^^^^^^
These direct algorithms compare two structures starting from the provided root nodes, using
their particular matching strategy. They are useful for immediate comparisons of similarity,
returning the maximal :func:`~.similarity.commutative_similarity` of their traversal, but they
do not consider non-root-to-root comparisons.

    .. autofunction:: topological_inclusion

    .. autofunction:: exact_ordering_inclusion

Substructure Inclusion-based Algorithms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Building off the `Direct Comparators`_, these algorithms answer higher level questions about
whether one structure is included in another.

    .. autofunction:: subtree_of

    .. autofunction:: find_matching_subtree_roots

Partial Subgraph Walks
**********************
    .. autofunction:: walk_with


Common Substructure-based Algorithms
--------------------------------------------

    .. autofunction:: maximum_common_subgraph

    .. autofunction:: n_saccharide_similarity

    .. autofunction:: distinct_fragments

Implementation Details of :func:`maximum_common_subgraph`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    .. autoclass:: MaximumCommonSubgraphSolver

    .. autoclass:: MaximumCommonSubtreeResults


Treelets
--------

    .. autofunction:: treelets

    .. autofunction:: treelet_enrichment

Implementation Details
^^^^^^^^^^^^^^^^^^^^^^
Treelet methods are built atop a more complex object system

    .. autoclass:: Treelet

    .. autoclass:: TreeletIterator
