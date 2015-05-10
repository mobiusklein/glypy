import os
import json
import logging
from collections import Iterable

from pygly2.utils import opener, pickle

from .spectrum_model import (neutral_mass,
                      mass_charge_ratio,
                      ObservedPrecursorSpectrum,
                      ObservedTandemSpectrum,
                      IonMatchAnnotation,
                      MSMSSqlDB)

from .constants import constants

decon_io_logger = logging.getLogger("DeconIO")


class DeconIOBase(object):

    @staticmethod
    def prepare(data_dict):
        holder = DeconIOBase(None)
        holder.data = data_dict
        return holder

    def __init__(self, file_path=None):
        self.file_path = file_path
        if file_path is None:
            self.data = dict()
        else:
            self._load(file_path)

    def _load(self, file_path):
        self.data = pickle.load(opener(file_path))

    def __iter__(self):
        return iter(self.data.items())

    def subset(self, indices):
        if not isinstance(indices, Iterable):
            indices = [indices]
        return self.prepare({i: self.data[i] for i in indices})

    def __getitem__(self, ix):
        return self.data[ix]

    def get_db_filename(self, file_path=None):
        if file_path is None:
            file_path = self.file_path + '.db'
        return file_path

    def to_db(self, file_path=None, overwrite=True):
        if file_path is None:
            file_path = self.file_path + '.db'
        if file_path is None:
            file_path = ":memory:"
        exists = os.path.exists(file_path)
        decon_io_logger.debug("Deconvoluted Ions Database file exists? %s", exists)
        db = MSMSSqlDB(file_path)
        if (exists and overwrite) or not exists:
            decon_io_logger.debug("Initializing database")
            db.init_schema()
            db.load_data(self.data.values())
            db.apply_indices()
        return db

    def index_by_scan_ids(self):
        index = {}
        for i, precursor in self:
            for scan in precursor.scans:
                index[scan['id']] = precursor
        return self.prepare(index)

    def to_json(self, file_path):
        json.dump({k: v.to_json() for k, v in self.data.items()}, open(file_path, 'wb'))


def default_loader(file_path):
    try:
        return pickle.load(opener(file_path, "rb"))
    except Exception, e:
        decon_io_logger.warning("Failed to open %s for default_loader", file_path, exc_info=e)


def default_serializer(obj, file_path):
    try:
        pickle.dump(obj, opener(file_path, "wb"))
    except Exception, e:
        decon_io_logger.warning("Failed to open %s for default_serializer, for %s", file_path, obj, exc_info=e)
