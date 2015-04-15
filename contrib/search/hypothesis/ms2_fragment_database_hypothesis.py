import os
from collections import defaultdict
from pygly2.algorithms import database
from .common_transforms import (monoisotopic_mass, reduced_mass,
                                permethelylated_mass,
                                deuteroreduced_permethylated_mass, derivatized_mass)

from pygly2.utils import identity

default_fragmentation_parameters = {
    "kind": "BYX",
    "max_cleavages": 2,
    "average": False,
    "charge": 0
}


def mass_transform(record, derivatize_fn=identity, adduct_mass=0, adduct_number=0):
    derivatize_fn(record.structure)
    mass = record.mass() + adduct_number * adduct_mass
    return mass


def extract_fragments(record, fragmentation_parameters=None):
    fragmentation_parameters = fragmentation_parameters or default_fragmentation_parameters
    container = defaultdict(dict)
    for frag in record.structure.fragments(**fragmentation_parameters):
        container[frag.kind][tuple(frag.included_nodes)] = frag
    return container


def prepare_database(in_database, out_database=None, mass_transform_parameters=None, fragmentation_parameters=None):
    if isinstance(in_database, str):
        in_database = database.RecordDatabase(in_database)
    if out_database is None:
        out_database_string = os.path.splitext(in_database.connection_string)[0]
        database.RecordDatabase(out_database_string, record_type=in_database.record_type)
    elif isinstance(out_database, str):
        out_database = database.RecordDatabase(out_database)
    for record in in_database:
        mass = mass_transform(record, **(mass_transform_parameters or {}))
        fragments = extract_fragments(record, fragmentation_parameters or default_fragmentation_parameters)
        record.fragments = fragments
        out_database.load_data([record], commit=False, mass_params={"override": mass})
    out_database.commit()
    return out_database
