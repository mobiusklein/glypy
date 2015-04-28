import re
from ..matching import NoShift, crossring_pattern, ppm_error, match_fragments
from collections import Counter

reducing_end_type = {"X", "Y", "Z"}
terminal_end_type = {"A", "B", "C"}


def find_matches(precursor, msms_db, shifts=None,
                 ms1_match_tolerance=1e-5,
                 ms2_match_tolerance=2e-5, ion_types="ABCXYZ"):
    '''
    Find all MS1 matches, find all MS2 matches in these matches, and merge the fragments found.
    '''
    shifts = shifts or [NoShift]
    results = []
    total_spectra = []
    precursor_ppm_errors = []
    scans_searched = set()
    i = 0
    ion_types = map(sorted, ion_types)
    precursor.fragments = [f for f in precursor.fragments if sorted(crossring_pattern.sub("", f.kind)) in (ion_types)]

    for shift in shifts:
        for row in msms_db.ppm_match_tolerance_search(precursor.intact_mass + shift.mass, ms1_match_tolerance):
            spectrum = msms_db.precursor_type.from_sql(row, msms_db)
            precursor_ppm_errors.append(ppm_error(precursor.mass() + shift.mass, spectrum.neutral_mass))
            scans_searched.update(spectrum.scan_ids)
            matches = match_fragments(precursor.fragments, spectrum.tandem_data,
                                      shifts=shifts, ms2_match_tolerance=ms2_match_tolerance)
            total_spectra.append(len(spectrum.tandem_data))
            results.append(matches)
            i += 1
    return results, total_spectra


def offset_frequency(matches):
    counter = Counter()
    for match in matches:
        counter[re.sub(r'[a-z0-9:]', '', match.match_key), match.charge] += 1

    return counter
