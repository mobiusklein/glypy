'''A set of algorithms for finding similar structures and substructures.
'''

from .inclusion import (
    topological_inclusion, exact_ordering_inclusion, subtree_of,
    find_matching_subtree_roots, walk_with, TopologicalInclusionMatcher)

from .common_subgraph import (
    maximum_common_subgraph, treelets, n_saccharide_similarity, distinct_fragments,
    treelet_enrichment)

from .common_subgraph import (
    Treelet, TreeletIterator, TreeletEnrichmentTest,
    MaximumCommonSubgraphSolver, MaximumCommonSubtreeResults)
