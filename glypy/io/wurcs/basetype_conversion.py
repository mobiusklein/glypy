from glypy.utils import invert_dict

# BaseTypeForRelativeConfiguration
# "3" and "4" in relative configuration are replace from
# "1" and "2" in "D" absolute configuration, respectively
base_type_to_descriptors = {
    "gro": "x",  # "xgro" has no relative configuration
    "thr": "34",
    "ery": "44",
    "ara": "344",
    "rib": "444",
    "lyx": "334",
    "xyl": "434",
    "all": "4444",
    "alt": "3444",
    "man": "3344",
    "glc": "4344",
    "gul": "4434",
    "ido": "3434",
    "tal": "3334",
    "gal": "4334"
}


descriptors_to_base_type = invert_dict(base_type_to_descriptors)
