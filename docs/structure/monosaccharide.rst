Monosaccharide
===================

Represents individual saccharide residues and their associated functions. These
are the basic unit of structural representation, possesing graph node-like properties.

.. currentmodule:: glypy.structure.monosaccharide


.. automodule:: glypy.structure.monosaccharide
    :no-members:

Monosaccharide Objects
----------------------

.. autoclass:: Monosaccharide
    :no-members:

.. contents:: Monosaccharide Methods
    :local:

Connection Enumeration
^^^^^^^^^^^^^^^^^^^^^^
        .. automethod:: Monosaccharide.parents
        .. automethod:: Monosaccharide.children
        .. automethod:: Monosaccharide.substituents

Adding and Removing Connections and Modifications
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        .. automethod:: Monosaccharide.add_monosaccharide
        .. automethod:: Monosaccharide.add_substituent
        .. automethod:: Monosaccharide.add_modification

        .. automethod:: Monosaccharide.drop_monosaccharide
        .. automethod:: Monosaccharide.drop_substituent
        .. automethod:: Monosaccharide.drop_modification

Position Occupancy
^^^^^^^^^^^^^^^^^^
        .. automethod:: Monosaccharide.is_occupied
        .. automethod:: Monosaccharide.open_attachment_sites
        .. automethod:: Monosaccharide.total_attachement_sites
        .. automethod:: Monosaccharide.occupied_attachment_sites

Equality Comparison
^^^^^^^^^^^^^^^^^^^
        Monosaccharide objects support equality comparison operators, ``==`` and ``!=``. They also support hashing,
        using the :func:`hash` value of :attr:`Monosaccharide.id`.

        .. automethod:: Monosaccharide.exact_ordering_equality
        .. automethod:: Monosaccharide.topological_equality
        .. automethod:: Monosaccharide.__eq__
        .. automethod:: Monosaccharide.__hash__

Serialization
^^^^^^^^^^^^^
        .. automethod:: Monosaccharide.serialize

        .. automethod:: Monosaccharide.register_serializer

        .. automethod:: Monosaccharide.available_serializers

Mass Spectrometry Utilities
^^^^^^^^^^^^^^^^^^^^^^^^^^^
        .. automethod:: Monosaccharide.total_composition
        .. automethod:: Monosaccharide.mass

Miscellaneous
^^^^^^^^^^^^^

        .. automethod:: Monosaccharide.clone


Explicit Uncyclized Reducing Ends and Labels
--------------------------------------------
    .. autoclass:: ReducedEnd
        :members:
