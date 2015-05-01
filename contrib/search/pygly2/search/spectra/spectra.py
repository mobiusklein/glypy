import itertools
import json
import re
import sqlite3
import logging
from copy import deepcopy

from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


from pygly2 import Composition

PROTON = Composition("H+").mass
db_logger = logging.getLogger(__name__)


def neutral_mass(mz, z):
    return (mz * z) - (z * PROTON)


def mass_charge_ratio(neutral_mass, z):
    return (neutral_mass + (z * PROTON)) / z


class MSMSSqlDB(object):
    def __init__(self, connection_string=":memory:", precursor_type=None):
        self.connection_string = connection_string
        self.connection = sqlite3.connect(connection_string)
        self.connection.row_factory = sqlite3.Row
        self.cursor = self.connection.cursor
        if precursor_type is None:
            precursor_type = ObservedPrecursorSpectrum
        self.precursor_type = precursor_type

    def init_schema(self):
        self.connection.executescript(self.precursor_type.sql_schema())
        self.connection.executescript(ObservedTandemSpectrum.sql_schema())
        self.connection.commit()

    def apply_indices(self):
        for ix_stmt in self.precursor_type.add_index():
            self.connection.executescript(ix_stmt)
        self.connection.commit()

    def load_data(self, precursors_list):
        for precursor in precursors_list:
            for stmt in precursor.to_sql():
                try:
                    self.connection.execute(stmt)
                except:
                    print(stmt)
                    raise
        self.connection.commit()

    def __getitem__(self, scan_id):
        results = []
        for row in self.execute('''SELECT *,
             ObservedPrecursorSpectrum.precursor_id, other_data FROM ObservedPrecursorSpectrum
             JOIN Scans ON ObservedPrecursorSpectrum.precursor_id =
             Scans.precursor_id WHERE scan_id = {0};'''.format(scan_id)):
            results.append(self.precursor_type.from_sql(row, self))
        return results

    def get_tandem_spectra(self, id):
        results = []
        for row in self.execute('''SELECT * FROM ObservedTandemSpectrum WHERE
                                   tandem_id = {}'''.format(id)):
            results.append(ObservedTandemSpectrum.from_sql(row, self))
        return results

    def __iter__(self):
        for row in self.execute("SELECT * FROM ObservedPrecursorSpectrum ORDER BY neutral_mass;"):
            yield self.precursor_type.from_sql(row, self)

    def execute(self, *args, **kwargs):
        return self.connection.execute(*args, **kwargs)

    def executemany(self, *args, **kwargs):
        return self.connection.executemany(*args, **kwargs)

    def executescript(self, *args, **kwargs):
        return self.connection.executescript(*args, **kwargs)

    def commit(self):
        self.connection.commit()

    def rollback(self):
        self.connection.rollback()

    def _find_boundaries(self, mass, tolerance):
        spread = mass * tolerance
        return (mass - spread, mass + spread)

    def ppm_match_tolerance_search(self, mass, tolerance, target_table="ObservedPrecursorSpectrum",
                                   precursor_id=None, mass_shift=0):
        boundaries = self._find_boundaries(mass + mass_shift, tolerance)
        results = self.execute("SELECT * FROM ObservedPrecursorSpectrum\
         WHERE neutral_mass BETWEEN %f AND %f;" % boundaries)
        for result in results:
            yield result


class Scan(object):
    def __init__(self, id, z, mz, **kwargs):
        self.scan_id = id
        self.z = z
        self.mz = mz

    def __repr__(self):
        rep = "Scan({scan_id} {mz}|{z})".format(**self.__dict__)
        return rep

    def to_json(self):
        return self.__dict__

    @classmethod
    def from_dict(cls, d):
        inst = cls(d['scan_id'], d["z"], d["mz"])
        return inst


class ObservedPrecursorSpectrum(object):
    __table_name__ = "ObservedPrecursorSpectrum"

    @classmethod
    def sql_schema(cls):
        return '''
    DROP TABLE IF EXISTS ObservedPrecursorSpectrum;
    CREATE TABLE ObservedPrecursorSpectrum (
    precursor_id INTEGER UNIQUE PRIMARY KEY NOT NULL,
    charge INTEGER,
    neutral_mass FLOAT,
    other_data TEXT
    );

    DROP TABLE IF EXISTS Scans;
    CREATE TABLE Scans(
    scan_id INTEGER UNIQUE PRIMARY KEY NOT NULL,
    mz FLOAT,
    z INTEGER,
    precursor_id INTEGER,
    FOREIGN KEY(precursor_id) REFERENCES ObservedPrecursorSpectrum(precursor_id)
    );
    '''

    def __init__(self, scan_data, scan_ids, z, neutral_mass, tandem_data, **data):
        self.charge = z
        self.scans = scan_data
        self.scan_ids = scan_ids
        self.neutral_mass = neutral_mass
        self.tandem_data = tandem_data
        self.other_data = data
        self.id = data.pop('id', None)

    def __repr__(self):
        return "<Scans: ({scans}), Neutral Mass: {neutral_mass}, Tandem Spectra: {num_tandem}>".format(
            num_tandem=len(self.tandem_data),
            **self.__dict__)

    def to_json(self):
        json_dict = deepcopy(self.__dict__)
        json_dict['tandem_data'] = [t.to_json()
                                    for t in json_dict['tandem_data']]
        return json_dict

    def to_sql(self):
        stmt = '''INSERT INTO {table} (precursor_id, neutral_mass, charge, other_data)
        VALUES ({id}, {neutral_mass}, {charge}, '{other_data}');'''.format(
            table="ObservedPrecursorSpectrum",
            id=self.id, neutral_mass=self.neutral_mass, charge=self.charge,
            other_data=json.dumps(self.other_data))
        yield stmt

        for tandem in self.tandem_data:
            yield (tandem.to_sql(self.id))

        for scan in self.scans:
            insert_stmt = '''INSERT INTO {table} (scan_id, mz, z, precursor_id)
            VALUES ({id}, {mz}, {z},{precursor_id});'''.format(
                table="Scans", precursor_id=self.id,  **scan)
            yield (insert_stmt)

    @classmethod
    def from_sql(cls, row, cursor):
        neutral_mass = row['neutral_mass']
        charge = row['charge']
        id = row['precursor_id']
        other_data = json.loads(row['other_data'])
        scans_rows = cursor.execute("SELECT * FROM Scans WHERE Scans.precursor_id={0};".format(id))
        scans = [Scan.from_dict(row) for row in scans_rows]
        scan_ids = [scan.scan_id for scan in scans]
        tandem_rows = cursor.execute("SELECT * FROM ObservedTandemSpectrum WHERE\
         ObservedTandemSpectrum.precursor_id={0};".format(id))
        tandem_spectra = [ObservedTandemSpectrum.from_sql(row, cursor) for row in tandem_rows]
        instance = cls(scans, scan_ids, charge, neutral_mass, tandem_spectra, **other_data)
        instance.id = id
        return instance

    @classmethod
    def add_index(cls):
        stmt = '''CREATE INDEX IF NOT EXISTS neutral_mass_index ON ObservedPrecursorSpectrum(neutral_mass DESC);
        '''
        yield stmt
        for stmt in ObservedTandemSpectrum.add_index():
            yield stmt


class ObservedTandemSpectrum(object):
    __table_name__ = "ObservedTandemSpectrum"

    @classmethod
    def sql_schema(cls):
        return '''
    DROP TABLE IF EXISTS ObservedTandemSpectrum;
    CREATE TABLE ObservedTandemSpectrum (
    tandem_id INTEGER UNIQUE PRIMARY KEY NOT NULL,
    intensity FLOAT,
    charge INTEGER,
    neutral_mass FLOAT,
    annotation TEXT,
    other_data TEXT,
    precursor_id INTEGER,
    FOREIGN KEY(precursor_id) REFERENCES ObservedPrecursorSpectrum(precursor_id));
    '''

    def __init__(self, mass, charge, intensity, id=None, annotation=None, **data):
        if annotation is None:
            annotation = data.pop("annotation", [])
        self.charge = charge
        self.mass = mass  # backwards compatibility
        self.neutral_mass = mass
        self.intensity = intensity
        self.id = data.pop("id", None) or id
        self.annotation = annotation
        self.other_data = data
        self.precursor_id = data.pop("precursor_id", None)

    def __repr__(self):
        return "<ObservedTandemSpectra {neutral_mass}, {charge}, {intensity}>".format(**self.__dict__)

    def to_json(self):
        return self.__dict__

    def to_sql(self, precursor_id=None):
        return '''INSERT INTO {table} (tandem_id, intensity, charge, neutral_mass, annotation, other_data, precursor_id)
        VALUES ({id}, {intensity}, {charge}, {neutral_mass}, '{annotation}', '{other_data}', {precursor_id});'''.format(
            table="ObservedTandemSpectrum",
            id=self.id, charge=self.charge, neutral_mass=self.neutral_mass, intensity=self.intensity,
            other_data=json.dumps(self.other_data), precursor_id=precursor_id, annotation=json.dumps(self.annotation)
        )

    def update(self, db_handle):
        return '''UPDATE {table} SET intensity={intensity}, charge={charge}, neutral_mass={neutral_mass},
                                     annotation={annotation}, other_data={other_data}
                                 WHERE tandem_id = {id};'''.format(
                                    table="ObservedTandemSpectrum",
                                    id=self.id, charge=self.charge, neutral_mass=self.neutral_mass,
                                    intensity=self.intensity, other_data=json.dumps(self.other_data),
                                    precursor_id=self.precursor_id, annotation=json.dumps(self.annotation)
        )

    @classmethod
    def from_sql(cls, row, cursor):
        mass = row['neutral_mass']
        charge = row['charge']
        other_data = json.loads(row['other_data'])
        annotation = json.loads(row["annotation"])
        intensity = row['intensity']
        id = row['tandem_id']
        instance = cls(mass, charge, intensity, id=id, annotation=annotation, **other_data)
        return instance

    @classmethod
    def add_index(cls):
        yield "CREATE INDEX IF NOT EXISTS parent_spectrum ON ObservedTandemSpectrum(precursor_id DESC);"


def extract_annotations(*spectra):
    return [dict(s.annotation) for s in spectra]


def matched_spectra(scan_ids, neutral_mass, ppm_error, match_key, **data):
    return dict(scan_ids=scan_ids,
                neutral_mass=neutral_mass,
                ppm_error=ppm_error,
                key=match_key,
                other_data=data)


def label_formatter(label):
    label = str(label).strip()
    if len(label) > 12:
        label = '\n'.join(
            s for s in re.split(r'([^\+-]+[\+-])', label) if len(s) != 0)

    return label, label.count("\n")


def plot_observed_spectra(spectra, title=None, label_formatter_fn=label_formatter,
                          include_target_ion_ladder=True, label_target_ion_ladder=False,
                          annotation_source=None):
    mzs = np.array(
        [mass_charge_ratio(spec.neutral_mass, spec.charge) for spec in spectra])
    intensities = np.array([spec.intensity for spec in spectra])
    charge = np.array([spec.charge for spec in spectra])
    fig = plt.figure()
    plt.xlabel("m/z")
    plt.ylabel("Relative Intensity")
    plt.bar(mzs, intensities, width=0.1, edgecolor='black')

    uppermost_text = 0
    label_map = extract_annotations(*spectra) if annotation_source is None else None
    if label_map is not None:
        track_colors = {}
        color_iter = itertools.cycle(colors)
        label_offset = itertools.cycle(label_shift)

        def get_color(key):
            try:
                return track_colors[key]
            except:
                track_colors[key] = color_iter.next()
                return track_colors[key]

        all_keys = {k for s in label_map for k in s.keys()}
        # Ensure that each key is associated with a color in case we have only
        # common ions
        map(get_color, all_keys)

        arrowprops = dict(width=0.01, headwidth=0.01, alpha=0.2)
        need_text_offset_reset = True
        for i, annots in enumerate(label_map):
            xy = (mzs[i], intensities[i])
            chg = charge[i]
            if need_text_offset_reset:
                offset_counter = label_offset.next()
                need_text_offset_reset = False

            all_same = len(set(map(str, annots.values()))) == 1 and len(
                annots) == len(all_keys) and (len(all_keys) != 1)
            if all_same:
                fmt_label, lines = label_formatter_fn(annots.values()[0])
                if chg > 1:
                    fmt_label = r"{1}$\tt{{^{{+{0}}}}}$".format(chg, fmt_label)
                plt.annotate(fmt_label, xy, arrowprops=arrowprops,
                             xytext=(-15, 10 + offset_counter + (lines * 12)),
                             textcoords='offset points',
                             color='black', fontweight=600,
                             va="top", ha="left")
                need_text_offset_reset = True
                offset_counter += (lines + 1) * 14
                uppermost_text = max(uppermost_text, offset_counter)
            else:
                for track, text in annots.items():
                    fmt_label, lines = label_formatter_fn(text)
                    if chg > 1:
                        fmt_label = r"{1}$\tt{{^{{+{0}}}}}$".format(chg, fmt_label)
                    plt.annotate(fmt_label, xy, arrowprops=arrowprops,
                                 xytext=(-15, 10 + offset_counter +
                                         (lines * 12)),
                                 textcoords='offset points',
                                 color=get_color(track),
                                 va="bottom", ha="left")
                    need_text_offset_reset = True
                    offset_counter += (lines + 1) * 14
                    uppermost_text = max(uppermost_text, offset_counter)
        # Generate the reference peaks expected from the matched precursors to
        # provide a visual guide for evaluating a match.
        # if include_target_ion_ladder:
        #     max_intensity = max(intensities)
        #     for target_seq in all_keys:
        #         # Parse the key which is the sequence into an object and generate the fragments
        #         # as a flat list.
        #         seq_obj = Sequence(target_seq)
        #         b_frags = list(
        #             itertools.chain.from_iterable(seq_obj.get_fragments('b')))
        #         y_frags = list(
        #             itertools.chain.from_iterable(seq_obj.get_fragments('y')))


        #         # No charge information in fragment, use 1?
        #         b_frag_mz = [mass_charge_ratio(f.mass, 1) for f in b_frags]
        #         y_frag_mz = [mass_charge_ratio(f.mass, 1) for f in y_frags]

        #         # Draw the bars for these reference ions. Make them visible but
        #         # don't wash out the actual observed ions
        #         plt.bar(b_frag_mz, [max_intensity] * len(b_frag_mz), width=0.01,
        #                 edgecolor=get_color(target_seq), alpha=0.1)
        #         plt.bar(y_frag_mz, [max_intensity] * len(y_frag_mz), width=0.01,
        #                 edgecolor=get_color(target_seq), alpha=0.1)

        #         if label_target_ion_ladder:
        #             layer = 1
        #             layer_size = 20
        #             label_offset = itertools.cycle(label_shift_small)
        #             for i, frag in enumerate(b_frags):
        #                 offset_counter = label_offset.next() + \
        #                     (layer * layer_size)
        #                 xy = (b_frag_mz[i], max_intensity)
        #                 label = frag.get_fragment_name()
        #                 fmt_label, lines = label_formatter_fn(label)
        #                 if len(fmt_label) > 4:
        #                     lines += 1
        #                 plt.annotate(fmt_label, xy, arrowprops=arrowprops,
        #                              xytext=(
        #                                  0, 0 + offset_counter + (lines * 12)),
        #                              textcoords='offset points',
        #                              color=get_color(target_seq),
        #                              va="bottom", ha="left")
        #                 uppermost_text = max(
        #                     uppermost_text, offset_counter + (lines * 12))

        #             for i, frag in enumerate(y_frags):
        #                 offset_counter = label_offset.next() + \
        #                     (layer * layer_size)
        #                 xy = (y_frag_mz[i], max_intensity)
        #                 label = frag.get_fragment_name()
        #                 fmt_label, lines = label_formatter_fn(label)

        #                 plt.annotate(fmt_label, xy, arrowprops=arrowprops,
        #                              xytext=(0, 0 + offset_counter +
        #                                      (lines * 12)),
        #                              textcoords='offset points',
        #                              color=get_color(target_seq),
        #                              va="bottom", ha="left")
        #                 uppermost_text = max(uppermost_text, offset_counter)

        #             layer += 1

        plt.figlegend(map(lambda x: mpatches.Patch(color=x), track_colors.values()),
                      track_colors.keys(), "upper right", fontsize=12)
        plt.ylim(0, (max(intensities) + uppermost_text) * 1.25)
        plt.xlim(0, (max(mzs)) * 1.15)
    fig.set_size_inches(12, 10)
    return fig

label_shift = [60, 0, 40, 0, 20]
label_shift_small = [60, 0, 40, 0, 20]
colors = ['red', 'blue', 'orange', 'green', 'teal', 'darkcyan']
