from math import fabs
from collections import defaultdict
from pygly2.utils import make_struct


FragmentMatch = make_struct("FragmentMatch", ["match_key", "mass", "ppm_error", "intensity", "charge"])


def ppm_error(x, y):
    return (x.mass - y.mass) / y.mass


def collect_similar_ions(fragments, tolerance=2e-5, redundant=True):
    groups = defaultdict(list)
    membership = dict()
    for index in fragments:
        for other in fragments:
            if other.name in membership and not redundant:
                continue
            if fabs(ppm_error(index, other)) < tolerance:
                groups[index.name].append(other)
                membership[other.name] = index.name
    return groups


def match_fragments(fragments, peak_list):
    matches = []
    for fragment in fragments:
        for peak in peak_list:
            match_error = fabs(ppm_error(fragment, peak))
            if match_error <= 2e-5:
                matches.append(FragmentMatch(fragment.name, peak.mass, match_error, peak.intensity, peak.charge))
    return matches


def merge_matches(matches):
    pass
