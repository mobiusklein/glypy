import argparse
import os

from pygly2.algorithms import database

from .matching import find_matches, DEFAULT_MS2_MATCH_TOLERANCE, DEFAULT_MS1_MATCH_TOLERANCE
from .spectra import bupid_topdown_deconvoluter, spectra
from .report import render


def main(structure_database, observed_data,
         ms1_match_tolerance=DEFAULT_MS1_MATCH_TOLERANCE,
         ms2_match_tolerance=DEFAULT_MS2_MATCH_TOLERANCE,
         ion_types="ABCXYZ"):
    if isinstance(structure_database, str):
        structure_database = database.RecordDatabase(structure_database)
    if isinstance(observed_data, str):
        if os.path.splitext(observed_data)[1] == ".yaml":
            observed_db = bupid_topdown_deconvoluter.BUPIDYamlParser(observed_data).to_db()
        elif os.path.splitext(observed_data)[1] == ".db":
            observed_db = spectra.MSMSSqlDB(observed_data)
        else:
            raise Exception("Cannot load data: {}".format(observed_data))
    elif isinstance(observed_data, spectra.MSMSSqlDB):
        observed_db = observed_data
    else:
        raise Exception("Cannot load data: {}".format(observed_data))
    matches = []
    for structure in structure_database:
        results = find_matches(structure, observed_db, ms1_match_tolerance,
                               ms2_match_tolerance, ion_types)
        if results.intact_structures_searched > 0:
            matches.append(results)
    return matches


app = argparse.ArgumentParser("pygly-ms2")
app.add_argument("-s", "--structure-database")
app.add_argument("-d", "--observed-data")
app.add_argument("-t1", "--ms1-tolerance", type=float, default=DEFAULT_MS1_MATCH_TOLERANCE)
app.add_argument("-t2", "--ms2-tolerance", type=float, default=DEFAULT_MS2_MATCH_TOLERANCE)
app.add_argument("-i", "--ion-types", action="append", default=[], help='Control which ion types (ABCXYZ) are considered. Defaults to all of them.')
app.add_argument("-o", "--output", default=None)


def taskmain():
    args = app.parse_args()
    # Flatten list
    args.ion_types = [c for i in args.ion_types for c in i]
    if len(args.ion_types) == 0:
        args.ion_types = ["A", "B", "C", "X", "Y", "Z"]
    print args
    matches = main(args.structure_database, args.observed_data,
                   args.ms1_tolerance, args.ms2_tolerance, args.ion_types)
    if args.output is None:
        args.output = os.path.splitext(args.structure_database)[0] + ".results.html"
    outfile = open(args.output, "w")
    outfile.write(render(matches, args.__dict__))
    outfile.close()
