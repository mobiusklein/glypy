'''
A parser for :title-reference:`GlycoCT{XML}` format.

This module provides a simple parser for GlycoCT XML. It has fewer
features than :mod:`glypy.io.glycoct`, and does not support writing.

It is included for historical purposes.
'''
# pragma: no cover

from operator import itemgetter
from collections import defaultdict
from glypy.utils import opener, StringIO, ET
from glypy.utils.multimap import OrderedMultiMap
from glypy.structure import monosaccharide, substituent, link
from glypy.structure.glycan import Glycan
from .format_constants_map import (anomer_map, superclass_map,
                                   link_replacement_composition_map, modification_map)
from .tree_builder_utils import try_int
from .file_utils import ParserError

basetype_unpacker = itemgetter("id", "anomer", "superclass", "ringStart", "ringEnd")

START = "!START"
RES = "RES"
LIN = "LIN"
REP = "REP"
ALT = "ALT"
UND = "UND"


class GlycoCTXMLError(ParserError):
    pass


class GlycoCTXMLSectionUnsupported(GlycoCTXMLError):
    pass


class GlycoCTXML(object):
    '''Parse :title-reference:`GlycoCT{XML}` text data into |Glycan| objects.

    The parser implements the :class:`Iterator` interface, yielding successive glycans
    from a text stream separated by empty lines.
    '''
    @classmethod
    def loads(cls, glycoct_str):
        '''Parse results from |str|'''
        rep = StringIO(glycoct_str)
        return cls(rep)

    def __init__(self, stream, structure_class=Glycan):
        self.graph = {}
        self.state = START
        self.handle = opener(stream, "r")
        self.counter = 0
        self.repeats = {}
        self.buffer = defaultdict(list)
        self.root = None
        self._iter = None
        self.structure_class = structure_class

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
        return next(self._iter)

    __next__ = next

    def _make_iterator(self):
        for evt, entity in ET.iterparse(self.handle, ("start", "end")):
            entity.tag = entity.tag.split("}")[-1]
            yield evt, entity

    def parse(self):
        for evt, entity in self._make_iterator():
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
            else:
                if entity.tag == "stemtype":
                    self.buffer['stem'].append(entity.attrib['type'][1:])
                    self.buffer['configuration'].append(entity.attrib['type'][0])
                elif entity.tag == "modification":
                    self.buffer['modification'].append(
                        (entity.attrib['type'], try_int(entity.attrib['pos_one'])))
                elif entity.tag == "basetype":
                    id, anomer, superclass, ring_start, ring_end = basetype_unpacker(entity.attrib)
                    superclass = superclass_map[superclass.upper()]
                    anomer = anomer_map[anomer]
                    id = int(id)
                    modifications = OrderedMultiMap()
                    mods = self.buffer.pop("modification", None)
                    if mods is not None:
                        for mod, pos in mods:
                            modifications[pos] = modification_map[mod]
                    is_reduced = "aldi" in modifications[1]
                    if is_reduced:
                        modifications.pop(1, "aldi")

                    residue = monosaccharide.Monosaccharide(
                        anomer=anomer, superclass=superclass, stem=self.buffer.pop('stem'),
                        configuration=self.buffer.pop("configuration"), ring_start=ring_start,
                        ring_end=ring_end, modifications=modifications, reduced=is_reduced, id=id)
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
                        yield self.structure_class(self.root)
                    self._reset()
                entity.clear()


def read(stream):
    '''
    A convenience wrapper for :class:`GlycoCTXML`

    Parameters
    ----------
    stream : file-like
        The text stream to parse structures from

    Returns
    -------
    :class:`~.GlycoCTXML`
    '''
    return GlycoCTXML(stream)


def load(stream, structure_class=Glycan, allow_multiple=True):
    """Read all structures from the provided text stream.

    Parameters
    ----------
    stream : file-like
        The text stream to parse structures from
    structure_class : type, optional
        :class:`~.Glycan` subclass to use
    allow_multiple: bool, optional
        If :const:`False`, makes loading multiple structures an error.

    Returns
    -------
    :class:`~.Glycan` or :class:`list` of :class:`~.Glycan`
    """

    g = GlycoCTXML(stream, structure_class=structure_class)
    first = next(g)
    if not allow_multiple:
        return first
    second = None
    try:
        second = next(g)
        collection = [first, second]
        collection.extend(g)
        return collection
    except StopIteration:
        return first


def loads(text, allow_multiple=True):
    """Read all structures from the provided text string.

    Parameters
    ----------
    text : str
        The text to parse structures from

    Returns
    -------
    :class:`~.Glycan` or :class:`list` of :class:`~.Glycan`
    """

    g = GlycoCTXML.loads(text)
    first = next(g)
    if not allow_multiple:
        return first
    second = None
    try:
        second = next(g)
        collection = [first, second]
        collection.extend(g)
        return collection
    except StopIteration:
        return first
