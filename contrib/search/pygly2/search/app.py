import re
import argparse
import os
import sys
import logging
try:
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                        datefmt="%H:%M:%S")
    logger = logging.getLogger()
except:
    pass
from pygly2.algorithms import database
from pygly2.utils import pickle

from .matching import (find_matches, DEFAULT_MS2_MATCH_TOLERANCE,
                       DEFAULT_MS1_MATCH_TOLERANCE, MassShift,
                       NoShift, collect_matched_scans, ResultsRecord,
                       ResultsDatabase)
from .spectra import bupid_topdown_deconvoluter, spectra
from .report import render


def main(structure_database, observed_data,
         ms1_match_tolerance=DEFAULT_MS1_MATCH_TOLERANCE,
         ms2_match_tolerance=DEFAULT_MS2_MATCH_TOLERANCE,
         shifts=None,
         ion_types="ABCXYZ", settings={}):
    if shifts is None:
        shifts = [NoShift]
    if isinstance(structure_database, str):
        structure_database = database.RecordDatabase(structure_database)
    if isinstance(observed_data, str):
        if os.path.splitext(observed_data)[1] == ".yaml":
            logger.info("Indexing Observed Ions")
            observed_db = bupid_topdown_deconvoluter.BUPIDYamlParser(observed_data).to_db()
        elif os.path.splitext(observed_data)[1] == ".db":
            observed_db = spectra.MSMSSqlDB(observed_data)
        else:
            raise Exception("Cannot load data: {}".format(observed_data))
    elif isinstance(observed_data, spectra.MSMSSqlDB):
        observed_db = observed_data
    else:
        raise Exception("Cannot load data: {}".format(observed_data))
    spectral_matches_list = []
    matches_db = ResultsDatabase("test.db", flag='w')
    spectral_match_db = spectra.MSMSSqlDB("test.db")
    matches = []
    logger.info("Begin Matching")
    for structure in structure_database:
        results, spectral_matches = find_matches(structure, observed_db,
                                                 shifts=shifts,
                                                 ms1_match_tolerance=ms1_match_tolerance,
                                                 ms2_match_tolerance=ms2_match_tolerance,
                                                 ion_types=ion_types)
        matches_db.load_data([ResultsRecord.from_base(results)], set_id=False)
        spectral_matches_list.append(spectral_matches)

        if results.intact_structures_searched > 0:
            matches.append(results)
    logger.info("Matching Complete")
    experimental_statistics = {}

    scans_matched, scans_not_matched = map(list, collect_matched_scans(matches, observed_db))
    experimental_statistics["count_scans_matched"] = len(scans_matched)
    experimental_statistics["count_scans_not_matched"] = len(scans_not_matched)

    matches_db.set_metadata("experimental_statistics", experimental_statistics)
    matches_db.set_metadata("scans_matched", scans_matched)
    matches_db.set_metadata("scans_not_matched", scans_not_matched)
    matches_db.set_metadata("spectral_matches", spectral_matches_list)
    matches_db.set_metadata("settings", settings)

    return matches, experimental_statistics, scans_matched, scans_not_matched, spectral_matches_list


app = argparse.ArgumentParser("pygly-ms2")
app.add_argument("-s", "--structure-database", help='Path to the structure databse to search against')
app.add_argument("-d", "--observed-data", help='Path to the observed ion data to search against')
app.add_argument("-t1", "--ms1-tolerance", type=float, default=DEFAULT_MS1_MATCH_TOLERANCE, help='PPM match tolerance for MS1 Matching')
app.add_argument("-t2", "--ms2-tolerance", type=float, default=DEFAULT_MS2_MATCH_TOLERANCE, help='PPM match tolerance for MS2 Matching')
app.add_argument("-i", "--ion-types", action="append", default=[],
                 help='Control which ion types (A,B,C,X,Y,Z and multiples of them for internal fragments)\
                 are considered as a comma separated list. Defaults to A,B,C,X,Y,Z.')
app.add_argument("-m", "--mass-shift", action='append', nargs=2, default=[])
app.add_argument("-o", "--output", default=None)


def taskmain():
    args = app.parse_args()
    logger.debug(args)
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
                MassShift(re.sub(r"'|\\", '', name), float(re.sub(r"'|\\", '', mass)))
                )
        args.mass_shift = shifts + [NoShift]
    results = main(args.structure_database, args.observed_data,
                   shifts=args.mass_shift,
                   ms1_match_tolerance=args.ms1_tolerance,
                   ms2_match_tolerance=args.ms2_tolerance,
                   ion_types=args.ion_types, settings=args.__dict__)
    matches, experimental_statistics, scans_matched, scans_not_matched, spectral_match_db = results
    if args.output is None:
        args.output = os.path.splitext(args.structure_database)[0] + ".results"
    store_file = open(args.output + ".pkl", 'wb')
    packed_results = {
        "matches": matches,
        "experimental_statistics": experimental_statistics,
        "settings": args.__dict__,
        "scans_matched": scans_matched,
        "scans_not_matched": scans_not_matched,
        "spectral_match_db": spectral_match_db
    }
    pickle.dump(packed_results, store_file)
    store_file.close()
    outfile = open(args.output + ".html", "w")
    outfile.write(render(**packed_results))
    outfile.close()

    store_file.close()
    logger.debug("Done")


def rerender(data_path, output_path=None):
    if output_path is None:
        output_path = os.path.splitext(data_path)[0] + ".html"
    results = pickle.load(open(data_path))
    with open(output_path, 'w') as outfile:
        outfile.write(render(**results))


def rerendermain():
    rerender(*sys.argv[1:])
