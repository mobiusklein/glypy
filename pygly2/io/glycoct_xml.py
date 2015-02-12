import logging
from operator import itemgetter
from collections import defaultdict
from ..utils import opener, StringIO, ET
from ..utils.multimap import OrderedMultiMap
from ..structure import monosaccharide, substituent, link, glycan
from .format_constants_map import anomer_map, link_replacement_composition_map
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


basetype_unpacker = itemgetter("id", "anomer", "superclass", "ringStart", "ringEnd")

START = "!START"
RES = "RES"
LIN = "LIN"
REP = "REP"
ALT = "ALT"
UND = "UND"


def try_int(v):
    try:
        return int(v)
    except:
        return v


class GlycoCTXMLError(Exception):
    pass


class GlycoCTXMLSectionUnsupported(GlycoCTXMLError):
    pass


class GlycoCTXML(object):
    @classmethod
    def loads(cls, glycoct_str):
        '''Parse results from |str|'''
        rep = StringIO(glycoct_str)
        return cls(rep)

    def __init__(self, stream):
        self.graph = {}
        self.state = START
        self.handle = opener(stream, "r")
        self.counter = 0
        self.repeats = {}
        self.buffer = defaultdict(list)
        self.root = None
        self._iter = None

    def _reset(self):
        self.graph = {}
        self.buffer = defaultdict(list)
        self.root = None
        self.state = START

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

    def parse(self):
        for evt, entity in ET.iterparse(self.handle, ("start", "end")):
            if evt != "end":
                if entity.tag == "sugar":
                    self.state = START
                elif entity.tag == 'residues':
                    if self.state == START:
                        self.state = RES
                    else:
                        raise GlycoCTXMLError("<residues> not the first section encountered")
                elif entity.tag == 'linkages':
                    if self.state != RES:
                        raise GlycoCTXMLError("<linkages> not following <residues>")
                elif entity.tag in {'repeat'}:
                    raise GlycoCTXMLSectionUnsupported("<{}> section is not supported".format(entity.tag))
                continue
            if entity.tag == "stemtype":
                self.buffer['stem'].append(entity.attrib['type'][1:])
                self.buffer['configuration'].append(entity.attrib['type'][0])
            elif entity.tag == "modification":
                self.buffer['modification'].append(
                    (entity.attrib['type'], try_int(entity.attrib['pos_one'])))
            elif entity.tag == "basetype":
                id, anomer, superclass, ring_start, ring_end = basetype_unpacker(entity.attrib)
                anomer = anomer_map[anomer]
                id = int(id)
                modifications = OrderedMultiMap()
                mods = self.buffer.pop("modification", None)
                if mods is not None:
                    for mod, pos in mods:
                        modifications[pos] = mod
                residue = monosaccharide.Monosaccharide(
                    anomer=anomer, superclass=superclass, stem=self.buffer.pop('stem'),
                    configuration=self.buffer.pop("configuration"), ring_start=ring_start,
                    ring_end=ring_end, modifications=modifications, id=id)
                self.graph[id] = residue
                if self.root is None:
                    self.root = residue
            elif entity.tag == "substituent":
                substituent_obj = substituent.Substituent(entity.attrib['name'])
                self.graph[int(entity.attrib['id'])] = substituent_obj
            elif entity.tag == "connection":
                parent_id = try_int(entity.attrib['parent'])
                child_id = try_int(entity.attrib['child'])
                parent_node = self.graph[parent_id]
                child_node = self.graph[child_id]
                link.Link(
                    parent_node, child_node, parent_position=self.buffer.pop("parent_position"),
                    child_position=self.buffer.pop("child_position"),
                    parent_loss=self.buffer.pop('parent_loss'), child_loss=self.buffer.pop('child_loss'),
                    id=self.buffer.pop('id'))
            elif entity.tag == "parent":
                self.buffer['parent_position'] = int(entity.attrib['pos'])
            elif entity.tag == "child":
                self.buffer['child_position'] = int(entity.attrib['pos'])
            elif entity.tag == "linkage":
                self.buffer["id"] = int(entity.attrib['id'])
                self.buffer["parent_loss"] = link_replacement_composition_map[entity.attrib['parentType']]
                self.buffer["child_loss"] = link_replacement_composition_map[entity.attrib['childType']]
            elif entity.tag == 'sugar':
                if self.root is not None:
                    yield glycan.Glycan(self.root)
                self._reset()


def read(stream):
    '''
    A convenience wrapper for :class:`GlycoCTXML`
    '''
    return GlycoCTXML(stream)


def loads(glycoct_str):
    '''
    A convenience wrapper for :meth:`GlycoCTXML.loads`
    '''

    return GlycoCTXML.loads(glycoct_str)
