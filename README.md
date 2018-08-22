[![https://img.shields.io/travis/mobiusklein/glypy.svg](https://img.shields.io/travis/mobiusklein/glypy.svg)](https://travis-ci.org/mobiusklein/glypy)
[![Documentation Status](https://readthedocs.org/projects/glypy/badge/?version=latest&style=flatsquare)](https://glypy.readthedocs.io/en/latest/?badge=latest)
# glypy
Glycan Analysis and Glycoinformatics Library for Python

Glycobiology is the study of the biological functions, properties, and structures of carbohydrate biomolecules,
also called *glycans*. These large, tree-like molecules are complex, having a wide variety of building blocks
as well as modifications and substitutions on those building blocks.

Much in the same way other bioinformatics libraries provide ways to represent DNA, RNA, or Protein sequences,
this library attempts to provide a representation of glycans. Much of the variation found in the
building blocks of these structures, monosaccharides, are caused by substitutions of functional groups on a
common core structure.

## Features

### Read In and Write Out Common Glycan Structure Formats and Sources
1. GlycoCT{condensed} (I/O)
2. GlycoCT{XML} (I)
3. GlycoMinds Linear Code (I/O)
4. IUPAC Three Letter Code (I/O)
5. Retreive data from `glySpace` using the web services provided by [GlyTouCan](https://glytoucan.org/), or run `SPARQL` queries directly on their Triplestore.

### Manipulate Glycan Structures
1. Traverse structures using either canonical or residue-level rule ordering.
2. Operate on monosaccharide and substituents as nodes and bonds as edges.
3. Add, remove, and modify these structures to alter glycan properties.
4. Identify substructures and motifs, classifying glycans.
5. Evaluate structural similarities with one of several ordering and comparator methods.
6. Plot tree structures with MatPlotLib, rendering using a configurable symbol nomenclature, such as SNFG, CFG, or IUPAC. Layout using vector graphics for perfect scaling.

### Example Use Cases
1. Calculate the mass of a native or derivatized glycan.
2. Generate glycosidic and cross ring cleavage fragments for a collection of glycan structures for performing MS/MS database search.
3. Perform substructure similarity searches with exact ordering or topological comparison and exact or fuzzy per-residue matching to classify a structure as an N-linked glycan.
4. Annotate MS spectra with glycan structures, labeling which peaks matched a database entry.
5. Download all N-Glycans from [GlyTouCan](https://glytoucan.org/)
6. Find all glycans in a list which contain a particular subtree, or find common subtrees in a database of glycans, performing treelet enrichment analysis.
7. Synthesize all possible glycans using a set of enzymes starting from a set of seed structures.
