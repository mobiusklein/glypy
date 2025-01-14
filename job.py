import glypy

print(glypy)
from glypy.structure import glycan_composition

x = glycan_composition.GlycanComposition(Hex=6, HexNAc=2)
print(x)