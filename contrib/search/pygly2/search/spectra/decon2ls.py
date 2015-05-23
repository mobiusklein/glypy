import csv
import os

from .spectrum_model import ObservedPrecursorSpectrum, Scan, MSMSSqlDB


def tryfloat(obj):
    try:
        return float(obj)
    except:
        return obj


class Decon2LSPeak(object):
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
    def from_csv(cls, stream, signal_to_noise_threshold=0):
        reader = csv.reader(stream)
        reader.next()
        for i, line in enumerate(reader):
            row = cls(*map(tryfloat, line))
            row.precursor_id = i
            if row.signal_to_noise > signal_to_noise_threshold:
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


class ScanInfo(object):
    @classmethod
    def sql_schema(cls):
        yield '''
        DROP TABLE IF EXISTS ScanInfo;
        CREATE TABLE ScanInfo(
            scan_id INT PRIMARY KEY,
            scan_time FLOAT,
            scan_type INT,
            base_peak_intensity FLOAT,
            total_ion_current FLOAT);
        '''

    def __init__(self, scan_id, scan_time, scan_type, base_peak_intensity, total_ion_current):
        self.scan_id = scan_id
        self.scan_time = scan_time
        self.scan_type = scan_type
        self.base_peak_intensity = base_peak_intensity
        self.total_ion_current = total_ion_current

    def to_sql(self):
        stmt = "INSERT OR REPLACE INTO ScanInfo \
        (scan_id, scan_time, scan_type, base_peak_intensity, total_ion_current)\
        VALUES ({scan_id}, {scan_time}, {scan_type}, {base_peak_intensity}, {total_ion_current})".format(
            **self.__dict__)
        yield stmt

    @classmethod
    def from_csv(cls, stream):
        reader = csv.DictReader(stream)
        for row in reader:
            inst = cls(row['scan_id'], row['scan_time'], row['scan_type'], row['bpi'], row['tic'])
            yield inst

    @classmethod
    def parse_to_sql(cls, stream, conn):
        cur = conn.cursor()
        for line in cls.sql_schema():
            cur.executescript(line)
        for block in cls.from_csv(stream):
            cur.execute(block.to_sql())


def parse_to_database(ms1_isos, db_name=None):
    if db_name is None:
        db_name = os.path.splitext(ms1_isos)[0] + '.db'
    spectrum_db = MSMSSqlDB(db_name)
    spectrum_db.apply_schema()
    with open(ms1_isos) as isos_file:
        for row in Decon2LSPeak.from_csv(isos_file):
            spec = row.to_observed_spectrum()
            spectrum_db.load_data((spec,), False)
    spectrum_db.commit()
    spectrum_db.apply_indices()
    return spectrum_db
