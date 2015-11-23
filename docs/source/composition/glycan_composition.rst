===================
Glycan Composition
===================

.. currentmodule:: glypy.composition.glycan_composition


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


A :class:`GlycanComposition` is just a bag of :class:`MonosaccharideResidue` and :class:`SubstituentResidue`, similar to |Composition|.
Its keys may be either :class:`MonosaccharideResidue` instances, :class:`SubstituentResidue` instances or strings which can be parsed by
:func:`from_iupac_lite`, and its values are integers. They may also be written to and from a string using
:meth:`GlycanComposition.serialize` and :meth:`GlycanComposition.parse`.

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




.. autoclass:: GlycanComposition
    :members:

.. autoclass:: MonosaccharideResidue
    :members:
    :exclude-members: name

.. autoclass:: SubstituentResidue
    :members:

.. autofunction:: to_iupac_lite
.. autofunction:: from_iupac_lite
.. autofunction:: from_glycan
.. autofunction:: parse
