from pygly2.composition import Composition
import warnings

do_warn = True


class CompositionIndex(dict):
    def __getitem__(self, key):
        try:
            composition_dict = super(CompositionIndex, self).__getitem__(key)
        except KeyError:
            composition_dict = {}
            if do_warn:
                warnings.warn("{key} could not be found. It may not have an explicit composition".format(key=key),
                              stacklevel=3)
        return Composition(composition_dict)

# Monosaccharide Residue Compositions
_monosaccharide_compositions = {
    "tri": {
        "C": 3,
        "H": 6,
        "O": 3
    },
    "tet": {
        "C": 4,
        "H": 8,
        "O": 4
    },
    "pen": {
        "C": 5,
        "H": 10,
        "O": 5
    },
    "hex": {
        "C": 6,
        "H": 12,
        "O": 6
    },
    "hep": {
        "C": 7,
        "H": 14,
        "O": 7
    },
    "oct": {
        "C": 8,
        "H": 16,
        "O": 8
    },
    "non": {
        "C": 9,
        "H": 18,
        "O": 9
    }
}

monosaccharide_composition = CompositionIndex(_monosaccharide_compositions)

#: Modification Compositions defined as dictionaries, specified as composition
#: deltas relative to the glycan structure
_modification_compositions = {
    "d": {
        "O": -1,
    },
    "a": {
        "H": -2,
        "O": 1
    },
    "aldi": {
        "H": 2
    },
    "keto": {},

    #: Special terms for experimenting with derivitization
    "_reserve": {},
    "_cleave": {"H": -1},
    "_reduce": {"H": 2}
}

modification_compositions = CompositionIndex(_modification_compositions)

# Specified as composition deltas relative to the glycan structure
_substituent_compositions = {
    "n_acetyl": {
        "C": 2,
        "H": 5,
        "N": 1,
        "O": 1,
    },
    "sulfate": {
        "S": 1,
        "O": 3,
        "H": 2
    },
    "amino": {
        "N": 1,
        "H": 3,
    },
    "n_glycolyl": {
        "C": 2,
        "H": 5,
        "O": 2,
        "N": 1
    },
    "phosphate": {
        'P': 1,
        'O': 3,
        'H': 3
    },
    "acetyl": {
        'H': 4,
        'C': 2,
        'O': 1
    },
    "methyl": {
        'H': 4,
        'C': 1
    },
    "n_sulfate": {
        "H": 3,
        'S': 1,
        'O': 3,
        'N': 1
    },
    "nitrate": {
        "N": 1,
        "O": 3
    },
    "iodo": {
        "I": 1,
        "H": 1,
    },
    "flouro": Composition("FH"),
    "imino": Composition("NH2"),
    "chloro": Composition("ClH"),
    "formyl": Composition("CHOH"),
    "bromo": Composition("BrH"),
    "pyruvate": Composition("COCOCH3H"),
    "n_alanine": {
      "C": 3,
      "H": 8,
      "N": 2,
      "O": 1
    },
    "n_dimethyl": Composition("N(CH3)2H"),
    "n_formyl": Composition("NHCHOH"),
    "n_methyl": Composition("NHCH3H"),
    "n_succinate": Composition("NCOCH2CH2COOHH"),
    "n_trifluoroacetyl": Composition("NHCOCF3H"),
    "thio": Composition("SHH")
}

substituent_compositions = CompositionIndex(_substituent_compositions)
