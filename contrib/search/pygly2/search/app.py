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
from .spectra import bupid_topdown_deconvoluter, spectrum_model
from .report import render


def main(structure_database, observed_data,
         ms1_match_tolerance=DEFAULT_MS1_MATCH_TOLERANCE,
         ms2_match_tolerance=DEFAULT_MS2_MATCH_TOLERANCE,
         shifts=None,
         ion_types="ABCXYZ", settings=None):
    if settings is None:
        settings = {}
    if shifts is None:
        shifts = [NoShift]
    if isinstance(structure_database, str):
        structure_database = database.RecordDatabase(structure_database)
    if isinstance(observed_data, str):
        if os.path.splitext(observed_data)[1] == ".yaml":
            logger.info("Indexing Observed Ions")
            observed_db = bupid_topdown_deconvoluter.BUPIDYamlParser(observed_data).to_db()
        elif os.path.splitext(observed_data)[1] == ".db":
            observed_db = spectrum_model.MSMSSqlDB(observed_data)
        else:
            raise Exception("Cannot load data: {}".format(observed_data))
    elif isinstance(observed_data, spectrum_model.MSMSSqlDB):
        observed_db = observed_data
    else:
        raise Exception("Cannot load data: {}".format(observed_data))
    store_file = settings.get("output") + '.db'
    try:
        os.remove(store_file)
    except:
        pass
    matches_db = ResultsDatabase(store_file, flag='w')
    spectral_match_db = spectrum_model.MSMSSqlDB(store_file)
    for line in observed_db.connection.iterdump():
        spectral_match_db.executescript(line)
    matches = []
    logger.info("Begin Matching")
    for structure in structure_database:
        structure = ResultsRecord.from_base(structure)
        results, spectral_matches = find_matches(structure, observed_db,
                                                 shifts=shifts,
                                                 ms1_match_tolerance=ms1_match_tolerance,
                                                 ms2_match_tolerance=ms2_match_tolerance,
                                                 ion_types=ion_types)
        matches_db.load_data([results], set_id=False)
        matches_db.bind(results)
        spectral_match_db.load_data(spectral_matches)
    logger.info("Matching Complete")
    matches_db.commit()
    matches_db.apply_indices()
    experimental_statistics = {}

    scans_matched, scans_not_matched = map(list, collect_matched_scans(matches, observed_db))
    experimental_statistics["count_scans_matched"] = len(scans_matched)
    experimental_statistics["count_scans_not_matched"] = len(scans_not_matched)

    matches_db.set_metadata("experimental_statistics", experimental_statistics)
    matches_db.set_metadata("scans_matched", scans_matched)
    matches_db.set_metadata("scans_not_matched", scans_not_matched)
    matches_db.set_metadata("settings", settings)

    return matches, experimental_statistics, scans_matched, scans_not_matched


app = argparse.ArgumentParser("pygly-ms2")
app.add_argument("-s", "--structure-database", help='Path to the structure databse to search against', required=True)
app.add_argument("-d", "--observed-data", help='Path to the observed ion data to search against', required=True)
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
    if args.output is None:
        args.output = os.path.splitext(args.structure_database)[0] + ".results"
    results = main(args.structure_database, args.observed_data,
                   shifts=args.mass_shift,
                   ms1_match_tolerance=args.ms1_tolerance,
                   ms2_match_tolerance=args.ms2_tolerance,
                   ion_types=args.ion_types, settings=args.__dict__)
    matches, experimental_statistics, scans_matched, scans_not_matched = results
    packed_results = {
        "matches": matches.from_sql(
            matches.execute(
                "SELECT * FROM {table_name} WHERE scan_count > 0 ORDER BY mass ASC")),
        "experimental_statistics": experimental_statistics,
        "settings": args.__dict__,
        "scans_matched": scans_matched,
        "scans_not_matched": scans_not_matched,
    }
    outfile = open(args.output + ".html", "w")
    outfile.write(render(**packed_results))
    outfile.close()

    logger.debug("Done")


def rerender(data_path, output_path=None):
    if output_path is None:
        output_path = os.path.splitext(data_path)[0] + ".html"
    db = ResultsDatabase(data_path)
    metadata = db.get_metadata()
    matches = db.from_sql(db.execute("SELECT * FROM {table_name} WHERE scan_count > 0 ORDER BY mass ASC"))
    with open(output_path, 'w') as outfile:
        outfile.write(render(matches=matches, **metadata))


def rerendermain():
    rerender(*sys.argv[1:])
