import requests
from lxml import etree
from pygly2.io import glycoct
from collections import namedtuple

Taxon = namedtuple("Taxon", ("tax_id", "name", 'entries'))
Aglyca = namedtuple("Aglyca", ("name", "reducing", 'entries'))
Entry = namedtuple("Entry", ("database", "id"))
Motif = namedtuple("Motif", ("name", "id", "motif_class"))


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


#: GlycomeDB supplies a detailed schema link which allows `lxml` to easily pull out
#: more than just the GlycoCT string. To download a more informative record, use :func:`get_record`

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


def search_substructure(glycan_obj):
    post_data = {
        "sequencetype": "glycoct_condenced",
        "sequence": glycan_obj.to_glycoct(),
        "msView": "yes"
    }
    session = requests.Session()
    r = session.post(_search_url, data=post_data)
    r.raise_for_status()
    tree = etree.fromstring(r.content)
    last_page = int([p.attrib for p in tree.findall('.//page')][-1]["number"])
    matches = [GlycomeDBSearchMatch(**s.attrib) for s in tree.findall(".//structure")]
    for m in matches:
        yield m
    for page_id in range(1, last_page):
        get_more = session.post(_get_more_results_url, data={"page": page_id})
        get_more.raise_for_status()
        tree = etree.fromstring(get_more.content)
        matches = [GlycomeDBSearchMatch(**s.attrib) for s in tree.findall(".//structure")]
        if len(matches) == 0:
            raise StopIteration(get_more.content)
        for m in matches:
            yield m


_search_url = "http://www.glycome-db.org/database/searchSubStructure.action"
_get_more_results_url = "http://www.glycome-db.org/database/showMultiStructure.action"


def search_minimum_common_substructure(glycan_obj, minimum_residues=1):
    post_data = {
        "sequencetype": "glycoct_condenced",
        "sequence": glycan_obj.to_glycoct(),
        "msView": "yes",
        "minResidueCount": minimum_residues
    }
    session = requests.Session()
    r = session.post(_mcs_search_url, data=post_data)
    r.raise_for_status()
    tree = etree.fromstring(r.content)
    last_page = int([p.attrib for p in tree.findall('.//page')][-1]["number"])
    matches = [GlycomeDBSearchMatch(**s.attrib) for s in tree.findall(".//structure")]
    for m in matches:
        yield m
    for page_id in xrange(1, last_page):
        get_more = session.post(_get_more_results_url, data={"page": page_id})
        get_more.raise_for_status()
        tree = etree.fromstring(get_more.content)
        matches = [GlycomeDBSearchMatch(**s.attrib) for s in tree.findall(".//structure")]
        if len(matches) == 0:
            raise StopIteration(get_more.content)
        for m in matches:
            yield m

_mcs_search_url = "http://www.glycome-db.org/database/searchMCS.action"


class GlycomeDBSearchMatch(object):
    def __init__(self, id, score, cross_reference, species_number, **kwargs):
        self.id = id
        self.score = score
        self.num_crossreferences = cross_reference
        self.num_species = species_number

    def get(self):
        return get(self.id)

    def get_record(self):
        return get_record(self.id)


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


