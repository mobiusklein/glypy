import os
from glypy.io import glycomedb
from glypy.algorithms import database

print("Moving 'glycomedb.db' to 'old_glycomedb.db'")
os.rename("glycomedb.db", "old_glycomedb.db")
print("Donwloading all structures to 'glycomedb.db'")
db = glycomedb.download_all_structures("glycomedb.db")
db._patch_querymethods()
print("Extracting Human glycans to 'human_glycans.db'")
hdb = database.dbopen("human_glycans.db", record_type=db.record_type, flag='w')
hdb.load_data(db.query_by_taxon_id(9606), set_id=False)
hdb.apply_indices()
print("Extracting Human N-glycans to 'human_n_glycans.db'")
nhdb = database.dbopen("human_n_glycans.db", record_type=db.record_type, flag='w')
nhdb.load_data(hdb.from_sql(hdb.execute("SELECT * FROM {table_name} WHERE is_n_glycan=1;")), set_id=False)
nhdb.apply_indices()
