import csv
from .spectrum_model import ObservedPrecursorSpectrum, Scan


def tryfloat(obj):
    try:
        return float(obj)
    except:
        return obj


class DeconRow(object):
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

    @staticmethod
    def from_csv(stream, noise_threshold=0):
        reader = csv.reader(stream)
        reader.next()
        for line in reader:
            row = DeconRow(*map(tryfloat, line))
            yield row

    def to_observed_spectrum(self):
        scan_ids = [self.scan_number]
        charge = self.charge
        neutral_mass = self.monoisotopic_abundance
        scan = Scan(self.scan_number, self.charge, self.mz)
        spectrum = ObservedPrecursorSpectrum([scan], scan_ids, charge, neutral_mass, [], **self.__dict__)
        return spectrum


class ResultsGroup(object):
    def __init__(self, ms1_score, observed_mass, composition, ppm_error, charge, scan_density,
                 average_a_plus_2_error, a_to_a_plus_2_ration, volume, signal_to_noise, centroid_error, centroid_scan,
                 max_scan_number, min_scan_number, calculated_mass, adduct_mass, num_adducts):
        pass
