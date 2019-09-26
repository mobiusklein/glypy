Glycan Structures
-----------------

Represent polysaccharide molecules and their associated functions

.. currentmodule:: glypy.structure.glycan

.. automodule:: glypy.structure.glycan
    :no-members:

    .. autoclass:: Glycan
        :no-members:

    Indexing
    ========

        .. automethod:: Glycan.get

        .. automethod:: Glycan.order

        .. automethod:: Glycan.reroot

        .. automethod:: Glycan.reindex

        .. automethod:: Glycan.deindex

        .. automethod:: Glycan.label_branches

    Traversal
    =========
        Glycan structures may be linear or branching, and can be traversed many ways. By
        default, a :meth:`Glycan.depth_first_traversal` is used, which will fully traverse
        one branch before visiting another, but other methods are available. Some methods
        simply control the behavior of the iterator but do not control the order of iteration,
        and take a ``method`` argument where either the name of the traversal method or a callable
        is specified.

        :class:`Glycan` objects implement the :class:`~.Iterable` interface, and their
        :meth:`__iter__` method :meth:`Glycan.depth_first_traversal`.

            .. automethod:: Glycan.depth_first_traversal

            .. automethod:: Glycan.breadth_first_traversal

            .. automethod:: Glycan.indexed_traversal

            .. automethod:: Glycan.iternodes

            .. automethod:: Glycan.iterlinks

            .. automodule:: Glycan._get_traversal_method

    Specialized Traversals
    ~~~~~~~~~~~~~~~~~~~~~~

            .. automethod:: Glycan.leaves


    Canonicalization
    ~~~~~~~~~~~~~~~~
        The same glycan structure can be constructed/written multiple ways, but they should all have the
        same representation. That representation is derived by applying a *canonicalization* algorithm to
        the structure, which will sort the branches of each node according to the order they should be
        traversed in.

        If a structure has been constructed manually, the user should call :meth:`Glycan.canonicalize`
        before assuming that identical structures will have the same traversal paths.

        .. automethod:: Glycan.canonicalize

    Ambiguous Structures
    ====================
        When a structure has unknown or ambiguous connections between is nodes, :class:`~.AmbiguousLink` instances
        may be used to express the possible options, or their locations may be expressed with an *unknown position*
        constant, represented with ``-1``. Two methods are included to detect these scenarios, and one is used to iterate
        over possible configuration states described by :class:`~.AmbiguousLink`.

        Support for ambiguous connections is only partial. For instance, :mod:`glypy` can read ``UND`` sections from
        :title:`GlycoCT`, but does not attempt to render them.

        .. automethod:: Glycan.ambiguous_links

        .. automethod:: Glycan.has_undefined_linkages

        .. automethod:: Glycan.iterconfiguration

    Serialization
    =============
        There are many ways to write glycan structures as text. By default, :mod:`glypy` will render
        :class:`~.Glycan` instances using :title:`GlycoCT`, but the :meth:`Glycan.serialize` method can
        be used to specify different serialization formats. For more information on those options, see
        :mod:`glypy.io`.

        When converting a :class:`Glycan` to a string, :meth:`Glycan.serialize` will be used with its
        default argument.

        .. automethod:: Glycan.serialize

        .. automethod:: Glycan.register_serializer

        .. automethod:: Glycan.available_serializers

    Sub-Structures
    ==============


