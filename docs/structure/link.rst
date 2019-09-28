Linkage and Bonds
=================

Represents the connections between residues, such as the bond between |Substituent| and |Monosaccharide| or the glycocidic bond between two |Monosaccharide| residues.

.. currentmodule:: glypy.structure.link

.. automodule:: glypy.structure.link
    :no-members:

    Explicit Linkages
    -----------------

    .. autoclass:: Link
        :no-members:

    .. automethod:: Link.__init__

    Connection State Management
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    .. automethod:: Link.apply

    .. automethod:: Link.break_link

    Allow Failure
    *************

    .. automethod:: Link.try_apply

    .. automethod:: Link.try_break_link

    Helpers
    *******

    .. automethod:: Link.refund


    Traversal
    ^^^^^^^^^

    .. automethod:: Link.to

    .. automethod:: Link.__iter__

    Connection Testing
    ^^^^^^^^^^^^^^^^^^
    These methods help determine if some object is related to a :class:`Link` instance, or
    can classify the type of linkage represented.

    .. automethod:: Link.is_parent
    .. automethod:: Link.is_child

    .. automethod:: Link.is_attached
    .. automethod:: Link.is_substituent_link
    .. automethod:: Link.is_bridge_link

    Ambiguity Testing
    *****************
    These methods are useful for switching behaviors when linkage is ambiguous. The base :class:`Link` class
    provides dummy implementations of these methods, but the :class:`AmbiguousLink` subclass implements them.

    .. automethod:: Link.is_ambiguous
    .. automethod:: Link.has_ambiguous_linkage
    .. automethod:: Link.has_ambiguous_termini

    Utilities
    ^^^^^^^^^
    .. automethod:: Link.clone

    :class:`Link` objects support equality comparison, but this tests the equality of their termini. To test
    whether two links' other properties are equal, see, :meth:`Link.trait_equality`

    .. automethod:: Link.__eq__
    .. automethod:: Link.trait_equality



    Ambiguous Linkages with Multiple Positions or Terminals
    -------------------------------------------------------

    Sometimes a link may be ambiguous, with several options for connection positions,
    or terminal residues. When representing these cases, :mod:`glypy` uses an :class:`AmbiguousLink`
    instead of a regular :class:`Link`. They inherit all the behaviors of :class:`Link`, but
    add several new fields for storing the possible options for each attribute, and provide
    means for switching amongst different configurations.

    .. autoclass:: AmbiguousLink
        :no-members:

    Configuration Combinations
    ^^^^^^^^^^^^^^^^^^^^^^^^^^

    .. automethod:: AmbiguousLink.reconfigure

    .. automethod:: AmbiguousLink.iterconfiguration

    .. automethod:: AmbiguousLink.find_open_position



