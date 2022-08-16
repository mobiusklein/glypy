# pragma: no cover
import gzip
import logging
import requests
from lxml import etree

from glypy.utils import StringIO
from glypy.io import glycoct
from glypy.algorithms.database import (Taxon, Aglyca, Motif,
                                       DatabaseEntry, GlycanRecord,
                                       GlycanRecordWithTaxon,
                                       RecordDatabase)

try:
    xrange
except NameError:
    xrange = range


import warnings
warnings.warn("Glycome-DB is no longer available. All Methods Will Fail.", DeprecationWarning)


logger = logging.getLogger(__name__)

cache = None


def download_all_structures(db_path, record_type=GlycanRecordWithTaxon):  # pragma: no cover
    response = requests.get(u'http://www.glycome-db.org/http-services/getStructureDump.action?user=eurocarbdb')
    response.raise_for_status()
    handle = gzip.GzipFile(fileobj=StringIO(response.content))
    xml = etree.parse(handle)
    db = RecordDatabase(db_path, record_type=record_type)
    misses = []
    i = 0
    for structure in xml.iterfind(".//structure"):
        try:
            glycomedb_id = int(structure.attrib['id'])
            i += 1
            glycoct_str = structure.find("sequence").text
            taxa = [Taxon(t.attrib['ncbi'], None, None) for t in structure.iterfind(".//taxon")]
            glycan = glycoct.loads(glycoct_str)
            if (glycoct.loads(str(glycan)).mass() - glycan.mass()) > 0.00001:
                raise Exception("Mass did not match on reparse")
            record = record_type(glycan, taxa=taxa, id=glycomedb_id)
            db.load_data(record, commit=False, set_id=False)
            if i % 1000 == 0:
                print(i, "Records parsed.")
        except Exception as e:
            misses.append((glycomedb_id, e))
            print(glycomedb_id, e)
    db.set_metadata("misses", misses)
    db.commit()
    return db


def set_cache(path):
    '''Set the path to the current cache for use with this client.
    This will create an instance of :class:`RecordDatabase` using the
    path provided.

    Any record lookups will first check the cache for that identifier, and
    if found will return from there. Any records not found will be fetched from
    the remote database and saved in the cache before returning.
    '''
    global cache
    logger.info("Setting glycomedb client cache to %r", path)
    cache = RecordDatabase(path)


def add_cache(record):
    if cache is None:
        return
    for stmt in record.to_sql():
        try:
            cache.connection.execute(stmt)
        except Exception as e:
            logger.error("An error occurred while adding %r", record, exc_info=e)


def check_cache(key):
    if cache is None:
        return False
    res = int(key) in cache
    if not res:
        logger.debug("Cache miss %r", key)
    else:
        logger.debug("Cache hit %r", key)
    return res


get_url_template = "http://www.glycome-db.org/database/showStructure.action?glycomeId={id}"
xpath = ".//condenced"


def _get_raw(id):
    r = requests.get(get_url_template.format(id=id))
    r.raise_for_status()
    tree = etree.fromstring(r.content)
    condensed = tree.find(xpath).text
    return condensed


def get(id):
    '''
    Get the structure for `id` from :title-reference:`GlycomeDB`.

    GlycomeDB supplies a detailed schema link which allows `lxml` to easily pull out
    more than just the GlycoCT string. To download a more informative record, use :func:`get_record`

    Parameters
    ----------
    id: str or int

    Returns
    -------
    Glycan
    '''
    if check_cache(id):
        return cache[id].structure
    r = requests.get(get_url_template.format(id=id))
    r.raise_for_status()
    tree = etree.fromstring(r.content)
    condensed = tree.find(xpath).text
    return glycoct.loads(condensed)


#: GlycomeDB supplies a detailed schema link which allows `lxml` to easily pull out
#: more than just the GlycoCT string. To download a more informative record, use :func:`get_record`

def get_record(id):
    '''
    Get the complete record for `id` from :title-reference:`GlycomeDB`.

    Returns
    -------
    GlycanRecord
    '''
    if check_cache(id):
        return cache[id]
    r = requests.get(get_url_template.format(id=id))
    r.raise_for_status()
    tree = etree.fromstring(r.content)
    return glycan_record_from_xml(tree, id)


def search_substructure(glycan_obj):  # pragma: no cover
    '''
    Select all structures which have contain `glycan_obj`.

    Parameters
    ----------
    glycan_obj: Glycan
        The glycan structure to search for

    Yields
    ------
    GlycomeDBSearchMatch
    '''
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


def search_minimum_common_substructure(glycan_obj, minimum_residues=1):  # pragma: no cover
    '''
    Select all structures which have a minimum common substructure of `minimum_residues`
    with `glycan_obj`.

    Parameters
    ----------
    glycan_obj: Glycan
        The glycan structure to search for
    minimum_residues: int
        The minimum number of residues to match

    Yields
    ------
    GlycomeDBSearchMatch
    '''
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
    '''
    Select all structures which are associated with the provided taxonomy id, or
    one of its children.

    Parameters
    ----------
    tax_id: int or str
        The taxonomy id to search for

    Yields
    ------
    GlycomeDBSearchMatch
    '''
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
            session.close()
            # raise StopIteration(get_more.content)
            return
        for m in matches:
            yield m
    session.close()


_taxa_search_url = "http://www.glycome-db.org/database/searchBySpecies.action"


def download_results(results, to=None):  # pragma: no cover
    if to is not None:
        if isinstance(to, str):
            set_cache(to)
        elif isinstance(to, RecordDatabase):
            global cache
            cache = to
    for r in results:
        print(r.id)
        try:
            r.get_record()
        except glycoct.GlycoCTError as e:
            logger.exception("Exception %r", e, exc_info=e)


class GlycomeDBSearchMatch(object):
    '''
    A search match which carries information about the scored match and the
    id number for retrieving the |Glycan| or |GlycanRecord|.
    '''
    def __init__(self, id, score, cross_reference, species_number, **kwargs):
        self.id = int(id)
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
    structure = glycoct.loads(xml_tree.find(xpath).text)
    taxa = [Taxon(t.attrib['ncbi'], t.attrib['name'], make_entries(t)) for t in xml_tree.findall(".//taxon")]
    aglycon = [Aglyca(t.attrib['name'].replace(
        "'", "`"), t.attrib['reducing'], make_entries(t)) for t in xml_tree.findall(".//aglyca")]
    motifs = [Motif(t.attrib['name'], t.attrib['id'], t.attrib['class']) for t in xml_tree.findall(".//motif")]
    dbxref = [e for c in [t.entries for t in taxa] + [t.entries for t in aglycon] for e in c]
    dbxref.append(DatabaseEntry("GlycomeDB", id))
    record = GlycanRecord(structure, motifs=motifs, dbxref=dbxref, aglycones=aglycon, taxa=taxa, id=id)
    record.id = id
    add_cache(record)
    return record

if __name__ == "__main__":
    import sys
    download_all_structures(sys.argv[1])
