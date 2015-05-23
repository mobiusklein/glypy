import argparse
import re
import os

from . import grouping_peaks

MassShift = grouping_peaks.MassShift
NoShift = grouping_peaks.NoShift
parse_to_database = grouping_peaks.parse_to_database
MSMSSqlDB = grouping_peaks.MSMSSqlDB
RecordDatabase = grouping_peaks.RecordDatabase


DEFAULT_MS1_MATCH_TOLERANCE = 1e-5

app = argparse.ArgumentParser("glycan-ms1")
app.add_argument("-s", "--structure-database", help='Path to the structure databse to search against', required=True)
app.add_argument("-d", "--observed-data", help='Path to the observed ion data to search against', required=True)
app.add_argument("-t1", "--ms1-tolerance", type=float, default=DEFAULT_MS1_MATCH_TOLERANCE, help='PPM match tolerance for MS1 Matching')
app.add_argument("-m", "--mass-shift", action='append', nargs=2, default=[])
app.add_argument("-o", "--output", default=None)


def taskmain(args):
    shifts = []
    for name, mass in args.mass_shift:
        shifts.append(
            MassShift(re.sub(r"'|\\", '', name), float(re.sub(r"'|\\", '', mass)))
            )
    args.mass_shift = shifts + [NoShift]

    main(args.observed_data,
         args.structure_database,
         args.mass_shift,
         args.ms1_tolerance,
         args.output)


def main(ms1_isos, hypothesis, adducts=None, ms1_match_tolerance=DEFAULT_MS1_MATCH_TOLERANCE, output_path=None):
    if adducts is None:
        adducts = [NoShift]
    if not isinstance(ms1_isos, MSMSSqlDB):
        ext = os.path.splitext(ms1_isos)[1]
        if ext == '.db':
            spectrum_db = MSMSSqlDB(ms1_isos)
        else:
            spectrum_db = parse_to_database(ms1_isos)
    else:
        spectrum_db = ms1_isos

    if isinstance(hypothesis, basestring):
        hypothesis = RecordDatabase(hypothesis, record_type=None)

    if output_path is None:
        output_path = hypothesis.connection_string
        if output_path != ":memory:":
            output_path = os.path.splitext(output_path) + '.composition_match.db'
    results = grouping_peaks.match_decon2ls_isos(spectrum_db, hypothesis, adducts, ms1_match_tolerance, output_path)
    return results


if __name__ == '__main__':
    taskmain(app.parse_args())
