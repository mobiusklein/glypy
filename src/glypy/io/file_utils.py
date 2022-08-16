import re
from glypy.utils import opener, Enum, root, StringIO
from glypy.structure.glycan import NamedGlycan
from six import add_metaclass


class ParserState(Enum):
    defline = 1
    sequence = 2


class ParserError(ValueError):
    pass


class FormatRegisteringMeta(type):
    def __new__(cls, name, parents, attrs):
        new_type = type.__new__(cls, name, parents, attrs)
        if not hasattr(cls, "registry"):
            cls.registry = {}

        if 'format_name' in attrs:
            cls.registry[attrs["format_name"]] = new_type
        return new_type


class TextFileParserBase(object):

    def __init__(self, *args, **kwargs):
        self._iter = None

    def parse(self):  # pragma: no cover
        raise NotImplementedError()

    def __iter__(self):
        self._iter = self.parse()
        return self._iter

    def __next__(self):
        if self._iter is None:
            iter(self)
        return next(self._iter)

    next = __next__

    @classmethod
    def loads(cls, text):
        return cls(StringIO(text))


@add_metaclass(FormatRegisteringMeta)
class FastaLikeFileParser(TextFileParserBase):
    format_name = 'fasta'

    def __init__(self, file_spec, processor):
        super(FastaLikeFileParser, self).__init__()
        self.state = ParserState.defline
        self.handle = opener(file_spec)
        self.defline = None
        self.sequence_chunks = []
        self.processor = processor

    def parse(self):
        for line in self.handle:
            line = line.lstrip()
            if self.state is ParserState.defline:
                if line[0] == ">":
                    self.defline = re.sub(r"[\n\r]", "", line[1:])
                    self.state = ParserState.sequence
                else:
                    continue
            else:
                if not re.match(r"^(\s+|>)", line):
                    self.sequence_chunks.append(re.sub(r"[\n\r]", "", line))
                else:
                    if self.defline is not None:
                        yield self.pack()

                    self.sequence_chunks = []
                    self.defline = None
                    self.state = ParserState.defline
                    if line[0] == '>':
                        self.defline = re.sub(r"[\n\r]", "", line[1:])
                        self.state = ParserState.sequence

        if len(self.sequence_chunks) > 0:
            yield self.pack()

    def pack(self):
        d = self.processor(''.join(self.sequence_chunks))
        name = self.defline
        root_node = root(d)
        return NamedGlycan(name=name, root=root_node, index_method='dfs')


@add_metaclass(FormatRegisteringMeta)
class StructurePerLineParser(TextFileParserBase):
    format_name = 'line'

    def __init__(self, file_spec, processor):
        super(StructurePerLineParser, self).__init__()
        self.handle = opener(file_spec)
        self.processor = processor

    def parse(self):
        for line in self.handle:
            data = self.processor(line)
            yield data


class ParserInterface(TextFileParserBase):

    def __init__(self, file_spec, file_type='line'):
        line_fetcher = FormatRegisteringMeta.registry[file_type]
        self.reader = line_fetcher(file_spec, processor=self.process_result)
        self.file_type = file_type

        super(ParserInterface, self).__init__()

    def parse(self):
        return self.reader.parse()

    @classmethod
    def loads(cls, text, file_type='line'):
        return cls(StringIO(text), file_type=file_type)

    @classmethod
    def load(cls, fd, file_type='line'):
        reader = cls(fd, file_type=file_type)
        data = list(reader)
        if len(data) == 1:
            return data[0]
        else:
            return data
