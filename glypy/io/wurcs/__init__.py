from glypy.io.file_utils import ParserInterface
from .parser import loads


class WURCSParser(ParserInterface):
    def process_result(self, line):
        structure = loads(line)
        return structure
