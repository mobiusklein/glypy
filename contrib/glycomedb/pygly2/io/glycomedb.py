from collections import namedtuple

Taxon = namedtuple("Taxon", ("tax_id", "name", 'entries'))
Aglyca = namedtuple("Aglyca", ("name", "reducing", 'entries'))
Entry = namedtuple("Entry", ("database", "id"))
Motif = namedtuple("Motif", ("name", "id", "motif_class"))

import requests
from lxml import etree
from pygly2.io import glycoct

get_url_template = "http://www.glycome-db.org/database/showStructure.action?glycomeId={id}"
xpath = ".//condenced"

def get(id):
    '''
    Get the structure for `id` from :title-reference:`GlycomeDB`.
    '''
    r = requests.get(get_url_template.format(id=id))
    r.raise_for_status()
    tree = etree.fromstring(r.content)
    condensed = tree.find(xpath).text
    return glycoct.loads(condensed).next()


#: GlycomeDB supplies a detailed schema link which allows `lxml` to easily pull out more than just the 
#: GlycoCT string. To download a more informative record, use :func:`get_record`

def get_record(id):
    '''
    Get the complete record for `id` from :title-reference:`GlycomeDB`. 

    Returns
    -------
    |Glycan|
    '''
    r = requests.get(get_url_template.format(id=id))
    r.raise_for_status()
    tree = etree.fromstring(r.content)
    return GlycomeDBRecord(tree, id)


search_url = "http://www.glycome-db.org/database/searchSubStructure.action"

def search_like(glycan_obj):
    post_data = {
        "sequencetype": "glycoct_condenced",
        "sequence": glycan_obj.to_glycoct(),
        "msView": "yes"
    }
    r = requests.post(search_url, data=post_data)
    return r


def make_entries(annotation):
    '''
    Helper method to create :class:`Entry` objects from <entry></entry> tags in :class:`Taxon` and :class:`Aglyca`.
    '''
    res = [Entry(node.attrib['database'], node.attrib['id']) for node in annotation.findall(".//entry")]
    return res

class GlycomeDBRecord(object):
    '''
    Describes a complete entry from :title-reference:`GlycomeDB`.

    Attributes
    ----------
    glycan: |Glycan| instance built from the structure stored in this record
    taxa: |list| of :class:`Taxon` objects
    aglycon: |list| of :class:`Aglyca` objects
    motifs: |list| of :class:`Motif` objects
    '''
    def __init__(self, xml_tree, id):
        self.glycan = glycoct.loads(xml_tree.find(xpath).text).next()
        self.taxa = [Taxon(t.attrib['ncbi'], t.attrib['name'], make_entries(t)) for t in xml_tree.findall(".//taxon")]
        self.aglycon = [Aglyca(t.attrib['name'], t.attrib['reducing'], make_entries(t)) for t in xml_tree.findall(".//aglyca")]
        self.motifs = [Motif(t.attrib['name'], t.attrib['id'], t.attrib['class']) for t in xml_tree.findall(".//motif")]
        self.id = id


