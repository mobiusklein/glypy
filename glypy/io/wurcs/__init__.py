from glypy.io.file_utils import ParserInterface
from .parser import loads
from .node_type import NodeTypeSpec
from .writer import dumps


class WURCSParser(ParserInterface):
    def process_result(self, line):
        structure = loads(line)
        return structure
