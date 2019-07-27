glypy Glycan Analysis and Glycoinformatics Library for Python
-------------------------------------------------------------

|https://img.shields.io/travis/mobiusklein/glypy.svg| |Documentation
Status|


Glycobiology is the study of the biological functions, properties, and
structures of carbohydrate biomolecules, also called *glycans*. These
large, tree-like molecules are complex, having a wide variety of
building blocks as well as modifications and substitutions on those
building blocks.

:mod:`glypy` is a Python library providing code for reading, writing, and
manipulating glycan structures, glycan compositions, monosaccharides, and
their substituents. It also includes interfaces to popular glycan structure
databases, `GlyTouCan <https://glytoucan.org/>`_ and `UnicarbKB <http://www.unicarbkb.org/>`_
using :term:`SPARQL` queries and an RDF-object mapper.

Example Use Cases
~~~~~~~~~~~~~~~~~

1. Traverse structures using either canonical or residue-level rule
   ordering.
2. Operate on monosaccharide and substituents as nodes and bonds as
   edges.
3. Add, remove, and modify these structures to alter glycan properties.
4. Identify substructures and motifs, classifying glycans.
5. Evaluate structural similarities with one of several ordering and
   comparator methods.
6. Plot tree structures with MatPlotLib, rendering using a configurable
   symbol nomenclature, such as SNFG, CFG, or IUPAC. Layout using vector
   graphics for lossless scaling.
7. Calculate the mass of a native or derivatized glycan.
8. Generate glycosidic and cross ring cleavage fragments for a
   collection of glycan structures for performing MS/MS database search.
9. Perform substructure similarity searches with exact ordering or
   topological comparison and exact or fuzzy per-residue matching to
   classify a structure as an N-linked glycan.
10. Annotate MS spectra with glycan structures, labeling which peaks
    matched a database entry.
11. Download all N-Glycans from `GlyTouCan <https://glytoucan.org/>`__
12. Find all glycans in a list which contain a particular subtree, or
    find common subtrees in a database of glycans, performing treelet
    enrichment analysis.
13. Synthesize all possible glycans using a set of enzymes starting from
    a set of seed structures.

.. |https://img.shields.io/travis/mobiusklein/glypy.svg| image:: https://img.shields.io/travis/mobiusklein/glypy.svg
   :target: https://travis-ci.org/mobiusklein/glypy
.. |Documentation Status| image:: https://readthedocs.org/projects/glypy/badge/?version=master
   :target: http://glypy.readthedocs.org/en/master/?badge=master


Citing
~~~~~~

If you use :mod:`glypy` in a publication please cite:

    Klein, J., & Zaia, J. (2019). glypy - An open source glycoinformatics library.
    Journal of Proteome Research.
    https://doi.org/10.1021/acs.jproteome.9b00367
