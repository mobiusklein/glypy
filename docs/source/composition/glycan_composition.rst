===================
Glycan Composition
===================

.. currentmodule:: glypy.composition.glycan_composition


There are some applications where we aren't interested in actual glycan structures, but just the
types and quantities of monosaccharides that compose them. The :class:`MonosaccharideResidue` and
:class:`GlycanComposition` serve this purpose.

A |Monosaccharide| is meant to be able to precisely describe where all of the bonds from
the carbon backbone are. A :class:`MonosaccharideResidue` abstracts away the notion of
position, and automatically deduct a water molecule from their :attr:`composition` to
account for a single incoming and a single outgoing glycosidic bond. Because they do not try to
completely describe the physical configuration of the molecule, :class:`MonosaccharideResidue`
removes information about ring type, anomericty, configuration, and optionally stem type. The level
of detail discarded is customizable in the :meth:`MonosaccharideResidue.from_monosaccahride` class method.


A :class:`GlycanComposition` is just a bag of :class:`MonosaccharideResidue`, similar to |Composition|.
Its keys may be either :class:`MonosaccharideResidue` instances or strings which can be parsed by
:func:`from_iupac_lite`, and its values are integers. They may also be written to and from a string using
:meth:`GlycanComposition.serialize` and :meth:`GlycanComposition.parse`.




.. autoclass:: GlycanComposition
    :members:

.. autoclass:: MonosaccharideResidue
    :members:
    :exclude-members: name

    .. py:method:: name()

        An alias of :func:`to_iupac_lite` called on `self`

.. autofunction:: to_iupac_lite
.. autofunction:: from_iupac_lite
.. autofunction:: from_glycan
.. autofunction:: parse
