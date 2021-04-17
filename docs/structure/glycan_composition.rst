Glycan Compositions and Residues
================================

.. currentmodule:: glypy.structure.glycan_composition


:class:`GlycanComposition`, :class:`MonosaccharideResidue`, and :class:`SubstituentResidue` are
useful for working with bag-of-residues where topology and connections are not relevant, but
the aggregate composition is known. These types work with a subset of the IUPAC three letter code
for specifying compositions.


A |Monosaccharide| is meant to be able to precisely describe where all of the bonds from
the carbon backbone are. A :class:`MonosaccharideResidue` abstracts away the notion of
position, and automatically deduct a water molecule from their :attr:`composition` to
account for a single incoming and a single outgoing glycosidic bond. Because they do not try to
completely describe the physical configuration of the molecule, :class:`MonosaccharideResidue`
removes information about ring type, anomericty, configuration, and optionally stem type. The level
of detail discarded is customizable in the :meth:`MonosaccharideResidue.from_monosaccahride` class method.


A :class:`GlycanComposition` is just a bag of :class:`MonosaccharideResidue` and :class:`SubstituentResidue`,
similar to |Composition|. Its keys may be either :class:`MonosaccharideResidue` instances,
:class:`SubstituentResidue` instances or strings which can be parsed by :func:`from_iupac_lite`, and its values
are integers. They may also be written to and from a string using :meth:`~.GlycanComposition.serialize` and
:meth:`~.GlycanComposition.parse`.

>>> g = GlycanComposition(Hex=3, HexNAc=2)
>>> g["Hex"]
3
>>> r = MonosaccharideResidue.from_iupac_lite("Hex")
>>> r
MonosaccharideResidue(Hex)
>>> g[r]
3
>>> import glypy
>>> abs(g.mass() - glypy.motifs["N-Glycan core basic 1"].mass()) < 1e-5
True
>>> g2 = GlycanComposition(Hex=5)
>>> g["@n-acetyl"] = -2 # Remove two n-acetyl groups from the composition
>>> abs(g.mass() - g2.mass()) < 1e-5
True

IUPAClite
---------

:title-reference:`IUPAClite` is a dialect of IUPAC for describing monosaccharides
while omitting some precise structural information from the grammar. It also includes
a compact line notation for glycan compositions.

A monosaccharide is denoted using IUPAC notation, omitting ring shape, anomeric state, chirality,
and modification positions (optionally). For example, ``a-D-Manp`` would be written ``Man``, or
``b-D-Glcp2NAc`` would be written ``GlcNAc``.

You can also use generic base types like ``Hex`` or ``Pen`` for example, to denote a six or five
carbon monosaccharide. The notation is composable, so you can specify an arbitrarily modified
monosaccharide, like ``HexNAc(S)`` to specify a sulfated HexNAc, using the parenthesized convention
that separates substituent groups, or ``dHexN`` for a deoxy-Hexosamine.

You can also define "floating" substituent groups by prefixing their full lowercase
names with an ``@``-sign, like ``@sulfate`` for sulfate or ``@acetyl`` for an acetyl group. Lastly,
it is also possible to denote an arbitrary named group using the notation ``#<name>#<chemical-formula>``,
though this should be used only when no other option is available.

See :func:`from_iupaclite` and :func:`to_iupaclite` implementations of monomer reading and writing.

A glycan composition is written as one or more ``<monosaccharide>:<count>`` occurrences separated by
a "; " (semi-colon + space), enclosed in "{ }". See :meth:`GlycanComposition.parse` and
:meth:`GlycanComposition.serialize` for implementation. A few examples are shown below:

.. code-block:: python
    :caption: IUPAClite glycan compositions

    {Hex:5; HexNAc:4; Neu5Ac:1}
    {Hex:5; HexNAc:4; Neu5Ac:2}
    {Fuc:1; Hex:5; HexNAc:4; Neu5Ac:2}
    {Fuc:2; Hex:6; HexNAc:5; Neu5Ac:1}
    {Fuc:1; Hex:6; HexNAc:5; Neu5Ac:2}


Residues
--------

.. autoclass:: MonosaccharideResidue
    :members:


Frozen Residues
~~~~~~~~~~~~~~~

:class:`MonosaccharideResidue` operations may require :class:`str` conversions which can be expensive.
Instead, use :class:`FrozenMonosaccharideResidue`, which once created is immutable, and substantially faster.

.. autoclass:: FrozenMonosaccharideResidue
    :members:

Substituent Residues
~~~~~~~~~~~~~~~~~~~~
.. autoclass:: SubstituentResidue
    :members:

Glycan Composition
------------------

.. autoclass:: GlycanComposition
    :members:

    .. automethod:: GlycanComposition.__init__

Frozen Composition
~~~~~~~~~~~~~~~~~~

:class:`GlycanComposition` objects automatically convert :class:`str` arguments to :class:`MonosaccharideResidue`
instances, which as previously mentioned, can be slow. If key objects will not be modified, the
:class:`FrozenGlycanComposition` is considerably faster for all operations. If both the keys themselves *and* the
values will not be modified after creation, the :class:`HashableGlycanComposition` is also useful and *hashable*.

.. autoclass:: FrozenGlycanComposition
    :members:

.. autoclass:: HashableGlycanComposition
    :members:

IUPAClite
---------

.. autofunction:: to_iupac_lite
.. autofunction:: from_iupac_lite
