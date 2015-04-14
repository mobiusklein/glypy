from collections import defaultdict
from math import fabs


from pygly2.structure import glycan
from pygly2.algorithms import database
from pygly2.composition import composition_transform
from pygly2.structure.monosaccharide import ReducedEnd


default_fragmentation_parameters = {
    "kind": "BYX",
    "max_cleavages": 2,
    "average": False,
    "charge": 0
}


def extract_fragments(record, fragmentation_parameters=None):
    fragmentation_parameters = fragmentation_parameters or default_fragmentation_parameters
    container = defaultdict(dict)
    for frag in record.structure.fragments(**fragmentation_parameters):
        container[frag.kind][tuple(frag.included_nodes)] = frag
    return container
