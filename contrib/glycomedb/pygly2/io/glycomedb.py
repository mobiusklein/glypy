import requests
from lxml import etree

from pygly2.io import glycoct
from pygly2.algorithms.database import (Taxon, Aglyca, Motif,
                                        DatabaseEntry, GlycanRecord,
                                        RecordDatabase)

cache = None


def set_cache(path):
    global cache
    cache = RecordDatabase(path)


def add_cache(record):
    if cache is None:
        return
    for stmt in record.to_sql():
        try:
            cache.connection.execute(stmt)
        except Exception, e:
            print(e, record.id)


def check_cache(key):
    if cache is None:
        return
    return key in cache


get_url_template = "http://www.glycome-db.org/database/showStructure.action?glycomeId={id}"
xpath = ".//condenced"


def get(id):
    '''
    Get the structure for `id` from :title-reference:`GlycomeDB`.
    '''
    if check_cache(id):
        return cache[id].structure
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
    if check_cache(id):
        return cache[id]
    r = requests.get(get_url_template.format(id=id))
    r.raise_for_status()
    tree = etree.fromstring(r.content)
    return glycan_record_from_xml(tree, id)


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
    try:
        last_page = int([p.attrib for p in tree.findall('.//page')][-1]["number"])
    except:
        last_page = 1
    matches = [GlycomeDBSearchMatch(**s.attrib) for s in tree.findall(".//structure")]
    for m in matches:
        yield m
    for page_id in range(2, last_page):
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
    try:
        last_page = int([p.attrib for p in tree.findall('.//page')][-1]["number"])
    except:
        last_page = 1
    matches = [GlycomeDBSearchMatch(**s.attrib) for s in tree.findall(".//structure")]
    for m in matches:
        yield m
    for page_id in xrange(2, last_page):
        get_more = session.post(_get_more_results_url, data={"page": page_id})
        get_more.raise_for_status()
        tree = etree.fromstring(get_more.content)
        matches = [GlycomeDBSearchMatch(**s.attrib) for s in tree.findall(".//structure")]
        if len(matches) == 0:
            raise StopIteration(get_more.content)
        for m in matches:
            yield m

_mcs_search_url = "http://www.glycome-db.org/database/searchMCS.action"


def search_by_species(tax_id):
    post_data = {
        "species": tax_id,
        "type": "ncbi",
        "subTree": "yes",
    }
    session = requests.Session()
    r = session.post(_taxa_search_url, data=post_data)
    r.raise_for_status()
    tree = etree.fromstring(r.content)
    try:
        last_page = int([p.attrib for p in tree.findall('.//page')][-1]["number"])
    except:
        last_page = 1
    matches = [GlycomeDBSearchMatch(**s.attrib) for s in tree.findall(".//structure")]
    for m in matches:
        yield m
    for page_id in xrange(2, last_page):
        get_more = session.post(_get_more_results_url, data={"page": page_id})
        get_more.raise_for_status()
        tree = etree.fromstring(get_more.content)
        matches = [GlycomeDBSearchMatch(**s.attrib) for s in tree.findall(".//structure")]
        if len(matches) == 0:
            raise StopIteration(get_more.content)
        for m in matches:
            yield m

_taxa_search_url = "http://www.glycome-db.org/database/searchBySpecies.action"


class GlycomeDBSearchMatch(object):
    '''
    A search match which carries information about the scored match and the
    id number for retrieving the |Glycan| or |GlycanRecord|.
    '''
    def __init__(self, id, score, cross_reference, species_number, **kwargs):
        self.id = id
        self.score = score
        self.num_crossreferences = cross_reference
        self.num_species = species_number

    def get(self):
        """Fetch the |Glycan| referenced by this search result

        Returns
        -------
        Glycan:
            The |Glycan| referenced
        """
        return get(self.id)

    def get_record(self):
        """Fetch the |GlycanRecord| referenced by this search result

        Returns
        -------
        GlycanRecord:
            The |GlycanRecord| referenced by this search result
        """
        return get_record(self.id)


def make_entries(annotation):
    '''
    Helper method to create :class:`Entry` objects from <entry></entry> tags in :class:`Taxon` and :class:`Aglyca`.
    '''
    res = [DatabaseEntry(node.attrib['database'], node.attrib['id']) for node in annotation.findall(".//entry")]
    return res


def glycan_record_from_xml(xml_tree, id):
    '''
    Converts an XML document and the associated database into an instance of
    `GlycanRecord`.

    Parameters
    ----------
    xml_tree: lxml.etree
        XML document to consume
    id:
        GlycomeDB id number to assign this record

    Returns
    -------
    GlycanRecord:
        Constructed record
    '''
    structure = glycoct.loads(xml_tree.find(xpath).text).next()
    taxa = [Taxon(t.attrib['ncbi'], t.attrib['name'], make_entries(t)) for t in xml_tree.findall(".//taxon")]
    aglycon = [Aglyca(t.attrib['name'].replace(
        "'", "`"), t.attrib['reducing'], make_entries(t)) for t in xml_tree.findall(".//aglyca")]
    motifs = [Motif(t.attrib['name'], t.attrib['id'], t.attrib['class']) for t in xml_tree.findall(".//motif")]
    dbxref = [e for c in [t.entries for t in taxa] + [t.entries for t in aglycon] for e in c]
    dbxref.append(DatabaseEntry("GlycomeDB", id))
    record = GlycanRecord(structure, motifs=motifs, dbxref=dbxref, aglycones=aglycon, taxa=taxa, id=id)
    add_cache(record)
    return record
