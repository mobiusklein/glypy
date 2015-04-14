import argparse
import sys
import csv

from pygly2.algorithms import database
from pygly2.composition import composition_transform
from pygly2.structure.monosaccharide import ReducedEnd


headings = [["Molecular Weight", "C", "Composition"], ["Adduct/Replacement", "Adduct Amount"]]


# Mass Transform Functions

def monoisotopic_mass(row):
    return row.mass()


def reduced_mass(row):
    row.structure.set_reducing_end(ReducedEnd())
    return monoisotopic_mass(row)


def permethelylated_mass(row):
    return derivatized_mass(row, "methyl")


def deuteroreduced_permethylated_mass(row):
    row.structure.set_reducing_end(ReducedEnd("H[2]H"))
    return permethelylated_mass(row)


def derivatized_mass(row, derivative):
    return composition_transform.derivatize(row.structure.clone(), derivative).mass() 


# Residue Composition Functions

def get_all_residue_types(db):
    residues = set()
    for record in db:
        residues |= set(record.monosaccharides)
    return residues


def pack_composition_string(record, residues):
    composition = record.monosaccharides
    return "[{}]".format(
        ';'.join(str(composition.get(residue, 0)) for residue in residues))


def prepare_row(row, residues, adduct_mass, num_adducts, mass_fn=lambda x: x.mass()):
    composition = row.monosaccharides
    return [mass_fn(row) + (adduct_mass * num_adducts), 0, pack_composition_string(row, residues)] +\
    [composition.get(residue, 0) for residue in residues] +\
    [adduct_mass if adduct_mass > 0.0 else "/0", num_adducts]


def hypothesis(db, outstream=sys.stdout, adduct_mass=0., num_adducts=0, mass_fn=monoisotopic_mass):
    residues = get_all_residue_types(db)
    columns = headings[0] + list(residues) + headings[1]
    writer = csv.writer(outstream)
    writer.writerow(columns)
    for row in db:
        writer.writerow(map(str, prepare_row(row, residues, adduct_mass, num_adducts, mass_fn=mass_fn)))
    try:
        outstream.close()
    except:
        pass


def main():
    db_path = sys.argv[1]
    hypothesis(database.RecordDatabase(db_path))


if __name__ == '__main__':
    main()
