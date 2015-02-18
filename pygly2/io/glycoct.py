import re
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

from ..utils import opener, StringIO
from ..utils.multimap import OrderedMultiMap
from ..structure import monosaccharide, substituent, link, constants, glycan
from .format_constants_map import anomer_map, superclass_map, link_replacement_composition_map

Glycan = glycan.Glycan
Monosaccharide = monosaccharide.Monosaccharide
Substituent = substituent.Substituent
Link = link.Link

START = "!START"
RES = "RES"
LIN = "LIN"
REP = "REP"
ALT = "ALT"
UND = "UND"

subsituent_start = "s"
base_start = "b"
repeat_start = "r"
alternative_start = "a"

res_pattern = re.compile(
    '''
    (?P<anomer>[abxo])?
    (?P<conf_stem>(?:-[dlx][a-z]+)+)?-?
    (?P<superclass>[A-Z]+)-?
    (?P<indices>[0-9x]+:[0-9x]+)
    (?P<modifications>(\|[0-9x]+:[0-9a-z]+)+)?
    ''', re.VERBOSE)

conf_stem_pattern = re.compile(r'(?P<config>[dlx])(?P<stem>[a-z]+)')

modification_pattern = re.compile(r"\|?(\d+):([^\|;\n]+)")

link_pattern = re.compile(
    r'''(?P<doc_index>\d+)?:
    (?P<parent_residue_index>\d+)
    (?P<parent_atom_replaced>[odhnx])
    \((?P<parent_attachment_position>-?\d+)[\+\-]
        (?P<child_attachment_position>-?\d+)\)
    (?P<child_residue_index>\d+)
    (?P<child_atom_replaced>[odhnx])
        ''', re.VERBOSE)


def try_int(v):
    try:
        return int(v)
    except:
        return v


class GlycoCTError(Exception):
    pass


class GlycoCTSectionUnsupported(GlycoCTError):
    pass


class GlycoCT(object):
    '''
    Simple State-Machine parser for condensed GlycoCT representations. Yields
    |Glycan| instances. 
    '''

    @classmethod
    def loads(cls, glycoct_str):
        '''Parse results from |str|'''
        rep = StringIO(glycoct_str)
        return cls(rep)

    def __init__(self, stream):
        '''
        Creates a parser of condensed GlycoCT. 

        Parameters
        ----------
        stream: basestring or file-like
            A path to a file or a file-like object to be processed
        '''
        self.graph = {}
        self.handle = opener(stream, "r")
        self.state = START
        self.counter = 0
        self.repeats = {}
        self.root = None
        self._iter = None

    def _read(self):
        for line in self.handle:
            for token in re.split(r"\s|;", line):
                logger.debug(token)
                if "" == token.strip():
                    continue
                yield token

    def _reset(self):
        self.graph = {}
        self.root = None
        #self.state = START

    def __iter__(self):
        '''
        Calls :meth:`parse` and stores it for reuse with :meth:`__next__`
        '''
        self._iter = self.parse()
        return self._iter

    def next(self):
        '''
        Calls :meth:`parse` if the internal iterator has not been instantiated

        '''
        if self._iter is None:
            iter(self)
        return self._iter.next()

    #: Alias for next. Supports Py3 Iterator interface
    __next__ = next

    def handle_residue_line(self, line):
        '''
        Handle a base line, creates an instance of |Monosaccharide|
        and adds it to :attr:`graph` at the given index.

        Called by :meth:`parse`
        '''
        _, ix, residue_str = re.split(r"^(\d+)b", line, maxsplit=1)
        residue_dict = res_pattern.search(residue_str).groupdict()

        mods = residue_dict.pop("modifications")
        modifications = OrderedMultiMap()
        if mods is not None:
            for p, mod in modification_pattern.findall(mods):
                modifications[try_int(p)] = constants.Modification[mod]

        residue_dict["modifications"] = modifications

        conf_stem = residue_dict.pop("conf_stem")
        if conf_stem is not None:
            config, stem = zip(*conf_stem_pattern.findall(conf_stem))
        else:
            config = ('x',)
            stem = ('x',)
        residue_dict['stem'] = stem
        residue_dict['configuration'] = config

        residue_dict["ring_start"], residue_dict["ring_end"] = list(map(
            try_int, residue_dict.pop("indices").split(":")))

        residue_dict['anomer'] = anomer_map[residue_dict['anomer']]
        residue_dict['superclass'] = superclass_map[residue_dict['superclass']]
        residue = monosaccharide.Monosaccharide(**residue_dict)
        self.graph[ix] = residue
        residue.id = int(ix)
        if self.root is None:
            self.root = residue

    def handle_residue_substituent(self, line):
        '''
        Handle a substituent line, creates an instance of |Substituent|
        and adds it to :attr:`graph` at the given index. The |Substituent| object is not yet linked
        to a |Monosaccharide| instance.

        Called by :meth:`parse`

        '''
        _, ix, subsituent_str = re.split(r"^(\d+)s:", line, maxsplit=1)
        sub = Substituent(subsituent_str.strip())
        self.graph[ix] = sub

    def handle_linkage(self, line):
        '''
        Handle a linkage line, creates an instance of |Link| and
        attaches it to the two referenced nodes in :attr:`graph`. The parent node is always
        an instance of |Monosaccharide|, and the child node
        may either be an instance of |Monosaccharide| or
        |Substituent| or |Monosaccharide|.

        Called by :meth:`parse`

        See also |Link| for more information on the impact of instantiating
        a |Link| object.
        '''
        link_dict = link_pattern.search(line)
        if link_dict is not None:
            link_dict = link_dict.groupdict()
        else:
            raise GlycoCTError("Could not interpret link", line)
        id = link_dict['doc_index']
        parent = self.graph[link_dict["parent_residue_index"]]
        parent_atom_replaced = link_replacement_composition_map[link_dict["parent_atom_replaced"]]
        parent_attachment_position = int(link_dict["parent_attachment_position"])

        child = self.graph[link_dict["child_residue_index"]]
        child_atom_replaced = link_replacement_composition_map[link_dict["child_atom_replaced"]]
        child_attachment_position = int(link_dict["child_attachment_position"])

        link_obj = Link(
            parent, child,
            parent_position=parent_attachment_position, child_position=child_attachment_position,
            parent_loss=parent_atom_replaced, child_loss=child_atom_replaced, id=id)

    def parse(self):
        '''
        Returns an iterator that yields each complete :class:`Glycan` instance
        from the underlying text stream.
        '''
        for line in self._read():
            if RES == line.strip():
                self.state = RES
                logger.debug("RES")
                if self.root is not None:
                    logger.debug("yielding root")
                    yield Glycan(self.root)
                    self._reset()
            elif LIN == line.strip():
                if self.state != RES:
                    raise GlycoCTError("LIN before RES")
                self.state = LIN

            elif REP == line.strip():
                raise GlycoCTSectionUnsupported(REP)
            elif ALT == line.strip():
                raise GlycoCTSectionUnsupported(ALT)
            elif UND == line.strip():
                raise GlycoCTSectionUnsupported(UND)

            elif re.search(r"^(\d+)b", line) and self.state == RES:
                logger.debug("handling residue")
                self.handle_residue_line(line)

            elif re.search(r"^(\d+)s:", line) and self.state == RES:
                logger.debug("handling subsituent")
                self.handle_residue_substituent(line)

            elif re.search(r"^(\d+):(\d+)", line) and self.state == LIN:
                logger.debug("handling linkage")
                self.handle_linkage(line)
            else:
                raise GlycoCTError("Unknown format error: {}".format(line))
        yield Glycan(self.root)


def read(stream):
    '''
    A convenience wrapper for :class:`GlycoCT`
    '''
    return GlycoCT(stream)


def loads(glycoct_str):
    '''
    A convenience wrapper for :meth:`GlycoCT.loads`
    '''

    return GlycoCT.loads(glycoct_str)
