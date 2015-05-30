# glypy
Glycan Analysis and Glycominformatics Library for Python

## Features

### Read in and write out common glycan structure formats
1. GlycoCT{condensed} (i/o)
1. GlycoCT{XML} (i)
1. GlycoMinds Linear Code (i/o)
1. IUPAC Three Letter Code (i)

### Manipulate glycan data structures like trees
2. Traverse the structure with common algorithms like *breadth-first* and *depth-first*.
2. Operate on monosaccharide and substituents as nodes and on bonds as edges.
2. Add, remove, and modify these structures to alter glycan properties.
3. Identify substructures and motifs, classifying glycans.
4. Score structural similarities with one of several ordering and comparator methods.
6. Plot tree structures with Matplotlib, rendering against any viable backend using configurable symbol nomenclature, such as Consortium for Functional Glycomics (CFG). Specialized SVG labeling for better web-interactivity.

### Example uses
3. Calculate the mass of a native or derivatized glycan.
3. Generate glycosidic and cross ring cleavage fragments for a collection of glycan structures for performing MS/MS database search.
3. Perform substructure similarity searches with exact ordering or topological comparison and exact or fuzzy per-residue matching to classify a structure as an N-linked glycan.
3. Annotate MS spectra with glycan structures, labeling which peaks matched a database entry.
