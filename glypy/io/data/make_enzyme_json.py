import json
from glypy.enzyme import EnzymeDatabase, EnzymeCommissionNumber

source = EnzymeDatabase._build()
keeper = EnzymeDatabase()
for k in list(source.direct_store.keys()):
    k = EnzymeCommissionNumber.parse(k)
    if k[0] in (1, 4, 5, 6):
        continue
    if k[0] == 2 and k[1] != 4:
        continue
    if k[0] == 3 and k[1] != 2:
        continue
    keeper.add(source[k])

with open("enzyme.json", 'w') as fp:
    keeper._dump(fp)
