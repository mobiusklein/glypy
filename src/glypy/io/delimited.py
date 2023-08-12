import re
from typing import Dict, Deque

import glypy


DEFAULT_MONOSACCHARIDES = {
    "H": glypy.monosaccharides.Hex,
    "N": glypy.monosaccharides.HexNAc,
    "F": glypy.monosaccharides.Fuc,
    "A": glypy.monosaccharides.NeuAc,
    "G": glypy.monosaccharides.NeuGc,
}


def parse_delimited(text: str, monosaccharides: Dict[str, glypy.Monosaccharide]=None, start='\(', end='\)'):
    if monosaccharides is None:
        monosaccharides = DEFAULT_MONOSACCHARIDES
    segments = Deque()
    current_depth = 0
    root = current_segment = {"node": None, "children": []}
    for c in re.compile(fr"({start}|[^{start}{end}]+|{end})").findall(text):
        if c == "(":
            current_depth += 1
            segments.append(current_segment)
            current_segment = {"node": None, "children": []}
            segments[-1]["children"].append(current_segment)
        elif c == ")":
            current_depth -= 1
            assert current_depth >= 0
            current_segment = segments.pop()
        else:
            current_segment["node"] = monosaccharides[c].clone()
    root_spec = root["children"][0]
    segments.clear()
    segments.append(root_spec)
    while segments:
        node_spec = segments.pop()
        for j, child_spec in enumerate(node_spec["children"], 2):
            node_spec["node"].add_monosaccharide(child_spec["node"], -1)
            segments.append(child_spec)
    return glypy.Glycan(root_spec["node"]).canonicalize()
