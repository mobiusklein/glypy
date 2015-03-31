# pygly2
Glycan Analysis and Glycominformatics Library for Python

## Features

### Read in and write out common glycan structure formats
1. GlycoCT{condensed} (i/o)
1. GlycoCT{XML} (i)
1. GlycoMinds Linear Code (i/o)
1. IUPAC Three Letter Code (i)

### Manipulate glycan data structures like trees
2. Traverse the structure with common algorithms like *bread-first* and *depth-first*.
2. Operate on monosaccharide and substituents as nodes and on bonds as edges.
2. Add, remove, and modify these structures to alter glycan properties
2. Plot tree structures with Matplotlib

### Example uses
3. Calculate the mass of a native or derivatized glycan
3. Generate glycosidic cleavage fragments, including internal fragments, for matching with tandem MS including cross ring cleavage.
3. Perform substructure similarity searches with exact ordering or topological comparison and exact or fuzzy per-residue matching
3. Annotate MS spectra with glycan structures
