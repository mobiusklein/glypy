import re
import os
import logging
try:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                        datefmt="%H:%M:%S")
except:
    pass
from pygly2.algorithms import database

from ..matching import (DEFAULT_MS2_MATCH_TOLERANCE,
                        DEFAULT_MS1_MATCH_TOLERANCE,
                        MassShift,
                        NoShift,
                        crossring_pattern,
                        ppm_error, match_fragments)
from ..spectra import bupid_topdown_deconvoluter, spectra

from collections import Counter

logger = logging.getLogger(__name__)

reducing_end_type = {"X", "Y", "Z"}
terminal_end_type = {"A", "B", "C"}


def main(structure_database, observed_data,
         ms1_match_tolerance=DEFAULT_MS1_MATCH_TOLERANCE,
         ms2_match_tolerance=DEFAULT_MS2_MATCH_TOLERANCE,
         shifts=None,
         ion_types="ABCXYZ"):
    """Summary

    Parameters
    ----------
    structure_database : RecordDatabase
        Description
    observed_data : MSMSSqlDB
        Description
    ms1_match_tolerance : float, optional
        Description
    ms2_match_tolerance : float, optional
        Description
    shifts : list of MassShift, optional
        Description
    ion_types : iterable, optional
        Description

    Returns
    -------
    list of GlycanRecord
    """
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
        if len(results[1]) > 0:
            matches.append(results)
    return matches


def find_matches(precursor, msms_db, shifts=None,
                 ms1_match_tolerance=1e-5,
                 ms2_match_tolerance=2e-5, ion_types="ABCXYZ"):
    '''Find all MS1 matches, find all MS2 matches in these matches, and merge the fragments found.

    Parameters
    ----------
    precursor : TYPE
        Description
    msms_db : TYPE
        Description
    shifts : TYPE, optional
        Description
    ms1_match_tolerance : float, optional
        Description
    ms2_match_tolerance : float, optional
        Description
    ion_types : str, optional
        Description

    Returns
    -------
    tuple of int, list, list
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
    return precursor.id, results, total_spectra


def offset_frequency(matches):
    """Count the occurences of each fragment type

    Parameters
    ----------
    matches : list

    Returns
    -------
    Counter
    """
    counter = Counter()
    for match in matches:
        key, shift = match.match_key.split(":")
        counter[re.sub(r'[a-z0-9]', '', key) + shift, match.charge] += 1

    return counter
