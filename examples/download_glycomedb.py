import os
import warnings
from glypy.io import glycomedb
from glypy.algorithms import database
import logging
logging.basicConfig(level=logging.DEBUG, filemode='w',
                    format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                    datefmt="%H:%M:%S")
warnings.simplefilter("error")
logging.info("Moving 'glycomedb.db' to 'old_glycomedb.db'")
try:
    os.remove("old_glycomedb.db")
except:
    pass
os.rename("glycomedb.db", "old_glycomedb.db")
logging.info("Downloading all structures to 'glycomedb.db'")
db = glycomedb.download_all_structures("glycomedb.db")
db._patch_querymethods()
logging.info("Extracting Human glycans to 'human_glycans.db'")
hdb = database.dbopen("human_glycans.db", record_type=db.record_type, flag='w')
hdb.load_data(db.query_by_taxon_id(9606), set_id=False)
hdb.apply_indices()
logging.info("Extracting Human N-glycans to 'human_n_glycans.db'")
nhdb = database.dbopen("human_n_glycans.db", record_type=db.record_type, flag='w')
nhdb.load_data(hdb.from_sql(hdb.execute("SELECT * FROM {table_name} WHERE is_n_glycan=1;")), set_id=False)
nhdb.apply_indices()
try:
    from glypy.search.hypothesis import ms2_fragment_database_hypothesis
    logging.info("De-duplicating human_n_glycans.db")
    keepers = ms2_fragment_database_hypothesis.duplicate_check(nhdb)
    ddnhdb = database.dbopen("human_n_glycans_deduplicated.db", record_type=db.record_type, flag='w')
    ddnhdb.load_data(keepers, set_id=False)
    ddnhdb.apply_indices()
except ImportError:
    pass
