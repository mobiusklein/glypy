import sqlite3
from collections import Counter

from pygly2.utils import pickle
from pygly2.io.nomenclature import identity


class GlycanRecord(object):
    __table_name__ = "GlycanRecord"

    __table__ = '''
    drop table if exists {table_name};
    create table {table_name}(
        glycan_id integer unique primary key not null,
        mass float not null,
        structure text not null
        /*rest*/
    );
    '''.format(table_name=__table_name__)

    def __init__(self, structure, motifs=None, dbxref=None, aglycones=None, **kwargs):
        self.structure = structure
        self.motifs = motifs or {}
        self.dbxref = dbxref or {}
        self.aglycones = aglycones or {}
        self.id = kwargs.get('id')

    def mass(self, average=False, charge=0, mass_data=None):
        return self.structure.mass(average=average, charge=charge, mass_data=mass_data)

    @property
    def monosaccharides(self):
        return Counter(map(identity.identify, self.structure))

    def __repr__(self):
        rep = '<GlycanRecord {id} {mass}>\n{glycoct}'.format(
                id=(self.id or ''), mass=self.mass(), glycoct=self.structure.to_glycoct()
            )
        return rep
