'''Implements a reader and writer for the `WURCS 2.0` format.
'''

from glypy.io.file_utils import ParserInterface
from .parser import loads
from .node_type import NodeTypeSpec
from .writer import dumps
from .utils import WURCSError
from .carbon_descriptors import CarbonDescriptors


class WURCSParser(ParserInterface):
    def process_result(self, line):
        structure = loads(line)
        return structure


__all__ = [
    "WURCSParser", "loads", "dumps",
    "NodeTypeSpec", "CarbonDescriptors",
    "WURCSError",
]
