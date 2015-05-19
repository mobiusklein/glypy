import yaml
import itertools

from . import DeconIOBase
from .spectrum_model import ObservedPrecursorSpectrum
from .spectrum_model import ObservedTandemSpectrum
from .spectrum_model import Scan
from . import neutral_mass
from .constants import constants as ms_constants


class BUPIDYamlParser(DeconIOBase):

    def __init__(self, file_path=None):
        super(BUPIDYamlParser, self).__init__(file_path)

    def _load(self, file_path):
        self.file_path = file_path
        stream = open(file_path, 'r')
        try:
            loader = yaml.CLoader(stream)
        except:
            loader = yaml.Loader(stream)
        raw_data = (loader.get_data())
        self.data = dict()
        self._build_spectra(raw_data)

    def _build_spectra(self, raw_data):
        ion_id = 0
        for tandem_ms_ind, peak_data in enumerate(raw_data['peaks']):
            scan_id_range = [scan["id"] for scan in peak_data["scans"]]
            # Treat the first scan as representative
            precursor = peak_data["scans"][0]
            precursor_mz = precursor["mz"]
            precursor_charge = precursor["z"]
            if precursor_charge > ms_constants.MAX_PRECURSOR_CHARGE_STATE:
                continue
            precursor_neutral_mass = neutral_mass(
                precursor_mz, precursor_charge)
            tandem_data = [ObservedTandemSpectrum(*ion, precursor_id=tandem_ms_ind) for
                           ion in itertools.izip(
                               peak_data["mass"], peak_data["z"], peak_data["intensity"])
                           if (ion[1] <= precursor_charge)  # and
                           #   (ion[1] <= ms_constants.MAX_FRAGMENT_CHARGE_STATE) and
                           # mass_charge_ratio(ion[0], ion[1]) <= ms_constants.MACHINE_ACQUISITION_RANGE
                           ]
            for tand in tandem_data:
                tand.id = ion_id
                ion_id += 1

            observed_spectra = ObservedPrecursorSpectrum([Scan(s['id'], s['z'], s['mz']) for s in peak_data["scans"]],
                                                         scan_id_range,
                                                         precursor_charge,
                                                         precursor_neutral_mass,
                                                         tandem_data)
            observed_spectra.id = tandem_ms_ind
            self.data[tandem_ms_ind] = observed_spectra
