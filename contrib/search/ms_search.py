import collections
from pygly2 import Composition
from pygly2.algorithms import database

PROTON = Composition("H+").mass


def neutral_mass(mz, z):
    return (mz * z) - (z * PROTON)


def mass_charge_ratio(neutral_mass, z):
    return (neutral_mass + (z * PROTON)) / z


default_fragmentation_parameters = {
    "kind": "BYX",
    "max_cleavages": 2,
    "average": False,
    "charge": 0
}
