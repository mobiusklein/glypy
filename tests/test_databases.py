import unittest
from pygly2.algorithms import database
from .common import load


class GlycanRecordTest(unittest.TestCase):

    def test_record_creation(self):
        rec = database.GlycanRecord(load("broad_n_glycan"))

        self.assertEqual(rec.structure.mass(), rec.mass())

        saccharides = {u'GlcNA': 6, u'Gal': 4, u'aMan': 2, u'Fuc': 1, u'Man': 1}
        self.assertEqual(saccharides, rec.monosaccharides)

        self.assertEqual(database.GlycanRecord.extract_composition(rec), '"Gal:4 aMan:2 Fuc:1 GlcNA:6 Man:1"')

    def test_to_sql(self):
        rec = database.GlycanRecord(load("broad_n_glycan"))
        db = database.RecordDatabase()
        db.create(load("broad_n_glycan"))
        self.assertEqual(db[1], rec)


class RecordDatabaseTest(unittest.TestCase):

    def test_load_data(self):
        rec = database.GlycanRecord(load("broad_n_glycan"))
        rec2 = database.GlycanRecord(load("complex_glycan"))
        db = database.RecordDatabase(records=[rec, rec2])
        self.assertTrue(db[1] == rec)
        self.assertTrue(db[2] == rec2)

    def test_ppm_search(self):
        rec = database.GlycanRecord(load("broad_n_glycan"))
        rec2 = database.GlycanRecord(load("complex_glycan"))
        db = database.RecordDatabase(records=[rec, rec2])
        self.assertEqual(rec, (db.ppm_match_tolerance_search(rec.mass(), 1e-5)).next())
