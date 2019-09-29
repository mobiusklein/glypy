Heuristic Similarity
====================

A collection of routines for doing fuzzy matching of monosaccharides, as well as a set of predicates
for classifying common properties of monosaccharides.

.. automodule:: glypy.algorithms.similarity
    :no-members:

.. contents::
    :local:

Core Heuristic
--------------

    .. autofunction:: monosaccharide_similarity

Commutative Options
^^^^^^^^^^^^^^^^^^^

    .. autofunction:: commutative_similarity

    .. autofunction:: commutative_similarity_score

    .. autofunction:: commutative_similarity_score_with_tolerance


Predicates
----------

    .. autofunction:: has_substituent
    .. autofunction:: has_modification
    .. autofunction:: has_monosaccharide
    .. autofunction:: is_reduced
    .. autofunction:: is_amine
    .. autofunction:: is_aminated
    .. autofunction:: is_generic_monosaccharide
    .. autofunction:: is_derivatized

Convenience Specializations
^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Some common predicates have been pre-bound using :class:`functools.partial`. Their names
    should hopefully be self-explanatory.

    .. function:: has_fucose
    .. function:: has_n_acetyl
    .. function:: is_acidic
    .. function:: is_sulfated


