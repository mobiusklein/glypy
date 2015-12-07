glypy Glycan Analysis and Glycoinformatics Library for Python
==============================================================


|https://img.shields.io/travis/mobiusklein/glypy.svg| |Documentation Status| 


Glycobiology is the study of the biological functions, properties, and
structures of carbohydrate biomolecules, also called *glycans*. These
large, tree-like molecules are complex, having a wide variety of
building blocks as well as modifications and substitutions on those
building blocks.

Much in the same way other bioinformatics libraries provide ways to
represent DNA, RNA, or Protein sequences, this library attempts to
provide a representation of glycans. Much of the variation found in the
building blocks of these structures, monosaccharides, are caused by
substitutions of functional groups on a common core structure.

Features
--------

Read in and write out common glycan structure formats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. GlycoCT{condensed} (i/o)
2. GlycoCT{XML} (i)
3. GlycoMinds Linear Code (i/o)
4. IUPAC Three Letter Code (i/o)

Manipulate glycan data structures like trees
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

2. Traverse the structure with common algorithms like *breadth-first*
   and *depth-first*, or use node-level information to choose a
   customized path.
3. Operate on monosaccharide and substituents as nodes and on bonds as
   edges.
4. Add, remove, and modify these structures to alter glycan properties.
5. Identify substructures and motifs, classifying glycans.
6. Score structural similarities with one of several ordering and
   comparator methods.
7. Plot tree structures with Matplotlib, rendering against any viable
   backend using configurable symbol nomenclature, such as Consortium
   for Functional Glycomics (CFG) or IUPAC text. Specialized SVG
   labeling for better web-interactivity.

Example Cases
~~~~~~~~~~~~~

3. Calculate the mass of a native or derivatized glycan.
4. Generate glycosidic and cross ring cleavage fragments for a
   collection of glycan structures for performing MS/MS database search.
5. Perform substructure similarity searches with exact ordering or
   topological comparison and exact or fuzzy per-residue matching to
   classify a structure as an N-linked glycan.
6. Annotate MS spectra with glycan structures, labeling which peaks
   matched a database entry.

.. |https://img.shields.io/travis/mobiusklein/glypy.svg| image:: https://img.shields.io/travis/mobiusklein/glypy.svg
   :target: https://travis-ci.org/mobiusklein/glypy
.. |Documentation Status| image:: https://readthedocs.org/projects/glypy/badge/?version=master
   :target: http://glypy.readthedocs.org/en/master/?badge=master
