import csv
import os

from pygly2.algorithms.database import RecordDatabase

from .spectrum_model import ObservedPrecursorSpectrum, Scan, MSMSSqlDB
from ..matching import MassShift, NoShift, make_struct, ppm_error


def tryfloat(obj):
    try:
        return float(obj)
    except:
        return obj


class MS1DeconRow(object):
    '''
    Describe a row of output from Decon2LS "DeconTools"
    '''
    def __init__(self, scan_number=0, charge=0, abundance=0.0, mz=0.0, fit=0.0,
                 average_mass_weight=0.0, monoisotopic_mass_weight=0.0, most_abundant_mass_weight=0.0,
                 full_width_half_height=0.0, signal_to_noise=0.0, monoisotopic_abundance=0.0,
                 monoisotopic_plus_2_abundance=0.0, flag=False, inference_score=0.0):
        self.scan_number = int(scan_number)
        self.charge = int(charge)
        self.abundance = abundance
        self.mz = mz
        self.fit = fit

        self.average_mass_weight = average_mass_weight
        self.monoisotopic_mass_weight = monoisotopic_mass_weight
        self.most_abundant_mass_weight = most_abundant_mass_weight

        self.full_width_half_height = full_width_half_height
        self.signal_to_noise = signal_to_noise

        self.monoisotopic_abundance = monoisotopic_abundance
        self.monoisotopic_plus_2_abundance = monoisotopic_plus_2_abundance

        self.flag = bool(flag)

        self.inference_score = inference_score

    @classmethod
    def from_csv(cls, stream, noise_threshold=0):
        reader = csv.reader(stream)
        reader.next()
        for i, line in enumerate(reader):
            row = cls(*map(tryfloat, line))
            row.precursor_id = i
            yield row

    @classmethod
    def from_observed_spectrum(cls, spectrum):
        return cls(**spectrum.other_data)

    def to_observed_spectrum(self):
        scan_ids = [self.scan_number]
        charge = self.charge
        neutral_mass = self.monoisotopic_mass_weight
        scan = Scan(self.scan_number, self.charge, self.mz)
        spectrum = ObservedPrecursorSpectrum([scan], scan_ids, charge, neutral_mass, [], **self.__dict__)
        spectrum.id = self.precursor_id
        return spectrum


class MS2DeconRow(object):
    def __init__(self, peak_index, scan_num, mz, intensity, full_width_half_height, signal_to_noise, feature_id):
        self.mz = mz
        self.intensity = intensity
        self.scan_num = scan_num
        self.full_width_half_height = full_width_half_height
        self.signal_to_noise = signal_to_noise
        self.feature_id = feature_id


class ResultsGroup(object):
    def __init__(self, ms1_score, observed_mass, composition, ppm_error, charge, scan_density,
                 average_a_plus_2_error, a_to_a_plus_2_ration, volume, signal_to_noise, centroid_error, centroid_scan,
                 max_scan_number, min_scan_number, calculated_mass, adduct_mass, num_adducts):
        pass


PrecursorMatch = make_struct("PrecursorMatch", ("match_key", "mass", "ppm_error",
                                                "intensity", "charge", "specturm_data", "peak_id"))


def match_ms1_hypothesis(ms1_isos, hypothesis, adducts=None, ms1_match_tolerance=1e-5):
    if adducts is None:
        adducts = [NoShift]
    if not isinstance(ms1_isos, MSMSSqlDB):
        spectrum_db = parse_to_database(ms1_isos)
    else:
        spectrum_db = ms1_isos

    if isinstance(hypothesis, basestring):
        hypothesis = RecordDatabase(hypothesis)

    for composition in hypothesis:
        matches = []
        scans_searched = set()
        charge_states = set()
        for adduct in adducts:
            for spectrum in spectrum_db.ppm_match_tolerance_search(
                    composition.intact_mass + adduct.mass, ms1_match_tolerance):
                match_ppm = ppm_error(composition.intact_mass + adduct.mass, spectrum.neutral_mass)
                match = PrecursorMatch(
                    str(composition.id) + ":" + adduct.name, spectrum.neutral_mass, match_ppm,
                    spectrum.get("abundance"), spectrum.charge, spectrum.other_data,
                    spectrum.id)
                matches.append(match)
                scans_searched.update(spectrum.scan_ids)
                charge_states.add(spectrum.charge)

        composition.matches = matches
        composition.scans_searched = scans_searched
        composition.charge_states = charge_states
        yield composition


def parse_to_database(ms1_isos, db_name=None):
    if db_name is None:
        db_name = os.path.splitext(ms1_isos)[0] + '.db'
    spectrum_db = MSMSSqlDB(db_name)
    spectrum_db.apply_schema()
    with open(ms1_isos) as isos_file:
        for row in MS1DeconRow.from_csv(isos_file):
            spec = row.to_observed_spectrum()
            spectrum_db.load_data((spec,), False)
    spectrum_db.commit()
    spectrum_db.apply_indices()
    return spectrum_db
