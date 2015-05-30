import sqlitedict
import os
from glypy.io import glycoct
from glypy.io.nomenclature import identity
from collections import Counter


def prepare_file(file_path):
    key = int(os.path.splitext(os.path.basename(file_path))[0])
    glycan = glycoct.read(file_path).next()
    composition = Counter(map(identity.identify, glycan))
    return key, {"key": key, "composition": composition, "glycan": glycan}


def build_database(data_path, db_path=None):
    if db_path is None:
        db_path = "./glycomedb.db"
    files = os.listdir(data_path)
    db = sqlitedict.open(db_path)
    total = float(len(files))
    for i, f in enumerate(files):
        if i % 10 == 0:
            print("%f%%" % (i/total * 100))
        if f[-3:] != 'txt':
            continue
        try:
            key, record = prepare_file(f)
            print(record)
            db[key] = record
        except (glycoct.GlycoCTSectionUnsupported, glycoct.GlycoCTError), e:
            print(e)
    db.commit()
    db.close()

if __name__ == '__main__':
    import sys
    build_database(sys.argv[1], None)
