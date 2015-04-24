import argparse
import os

from pygly2.algorithms import database
from pygly2.utils import pickle

from .matching import find_matches, DEFAULT_MS2_MATCH_TOLERANCE, DEFAULT_MS1_MATCH_TOLERANCE, MassShift, NoShift
from .spectra import bupid_topdown_deconvoluter, spectra
from .report import render


def main(structure_database, observed_data,
         ms1_match_tolerance=DEFAULT_MS1_MATCH_TOLERANCE,
         ms2_match_tolerance=DEFAULT_MS2_MATCH_TOLERANCE,
         shifts=None,
         ion_types="ABCXYZ"):
    if shifts is None:
        shifts = [NoShift]
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
        results = find_matches(structure, observed_db,
                               shifts=shifts,
                               ms1_match_tolerance=ms1_match_tolerance,
                               ms2_match_tolerance=ms2_match_tolerance,
                               ion_types=ion_types)
        if results.intact_structures_searched > 0:
            matches.append(results)
    return matches


app = argparse.ArgumentParser("pygly-ms2")
app.add_argument("-s", "--structure-database")
app.add_argument("-d", "--observed-data")
app.add_argument("-t1", "--ms1-tolerance", type=float, default=DEFAULT_MS1_MATCH_TOLERANCE)
app.add_argument("-t2", "--ms2-tolerance", type=float, default=DEFAULT_MS2_MATCH_TOLERANCE)
app.add_argument("-i", "--ion-types", action="append", default=[], help='Control which ion types (ABCXYZ) are considered. Defaults to all of them.')
app.add_argument("-m", "--mass-shift", action='append', nargs=2, default=[])
app.add_argument("-o", "--output", default=None)


def taskmain():
    args = app.parse_args()
    # Flatten list
    args.ion_types = [c for i in args.ion_types for c in i.split(',')]
    if len(args.ion_types) == 0:
        args.ion_types = ["A", "B", "C", "X", "Y", "Z"]
    if len(args.mass_shift) == 0:
        args.mass_shift = [NoShift]
    else:
        shifts = []
        for name, mass in args.mass_shift:
            shifts.append(
                MassShift(name.replace("'", ""), float(mass.replace("'", "")))
                )
        args.mass_shift = shifts + [NoShift]
    matches = main(args.structure_database, args.observed_data,
                   shifts=args.mass_shift,
                   ms1_match_tolerance=args.ms1_tolerance,
                   ms2_match_tolerance=args.ms2_tolerance,
                   ion_types=args.ion_types)
    if args.output is None:
        args.output = os.path.splitext(args.structure_database)[0] + ".results.html"
    outfile = open(args.output, "w")
    outfile.write(render(matches, args.__dict__))
    outfile.close()
    store_file = open(os.path.splitext(args.output)[0] + ".pkl", 'wb')
    pickle.dump({"matches": matches, "settings": args.__dict__}, store_file)
    store_file.close()


def rerender(data_path, output_path=None):
    if output_path is None:
        output_path = os.path.splitext(data_path)[0] + ".html"
    results = pickle.load(open(data_path))
    with open(output_path, 'w') as outfile:
        outfile.write(render(**results))


def rerendermain():
    import sys
    rerender(*sys.argv[1:])
