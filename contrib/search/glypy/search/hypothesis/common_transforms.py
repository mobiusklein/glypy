from glypy.composition import composition_transform
from glypy.structure.monosaccharide import ReducedEnd


# Mass Transform Functions

def monoisotopic_mass(row):
    return row.mass()


def reduced_mass(row):
    row.structure.set_reducing_end(ReducedEnd())
    return monoisotopic_mass(row)


def derivatized_mass(row, derivative):
    return composition_transform.derivatize(row.structure, derivative).mass()


def permethelylated_mass(row):
    return derivatized_mass(row, "methyl")


def deuteroreduced_permethylated_mass(row):
    row.structure.set_reducing_end(ReducedEnd("H[2]H"))
    return permethelylated_mass(row)
