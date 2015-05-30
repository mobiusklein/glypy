import json
from collections import defaultdict

import glypy

def main():
    store = defaultdict(list)
    for name, glycoct in glypy.monosaccharides.items():
        store[str(glycoct)].append(name)
    json.dump(store.values(), open("monosaccharide_synonyms.json", 'wb'), indent=2)

def schematize():
    store = defaultdict(list)
    for name, glycoct in glypy.monosaccharides.items():
        store[str(glycoct)].append(name)

    json.dump([{"structure": glycoct, "names": names} for glycoct, names in store.items()], 
        open("monosaccharide_schematic.json", 'wb'), indent=2)



if __name__ == '__main__':
    main()
    schematize()
