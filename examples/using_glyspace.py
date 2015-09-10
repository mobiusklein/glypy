from glypy.io.glyspace import GlySpaceRDF, NSGlyTouCan, NSGlycan, NSGlycoinfo, NSGlycomeDB, NSSKOS

client = GlySpaceRDF()

glycomedb_to_glyspace_map = {}

for glyspace, glycomedb in client.subject_objects(NSSKOS["exactMatch"]):
    glyspace = glyspace.split("/")[-1]
    glycomedb = glycomedb.split("/")[-1]
    glycomedb_to_glyspace_map[glycomedb] = glyspace

for k, v in glycomedb_to_glyspace_map.items():
    print "%s, %s" % (k, v)
