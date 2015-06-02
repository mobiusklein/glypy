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

# http://glytoucan.org/glyspace/documentation/apidoc.html
# http://code.glytoucan.org/system/glyspace

logger = logging.getLogger(__name__)

cache = None


def set_cache(path, record_type=GlycanRecord):
    global cache
    logger.info("Setting glyspace client cache to %r", path)
    cache = RecordDatabase(path, record_type=record_type)


def add_cache(record):
    if cache is None:
        return
    for stmt in record.to_sql():
        try:
            cache.connection.execute(stmt)
        except Exception, e:
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


get_url_template = "http://glytoucan.org/glyspace/service/glycans/{accession}.json"
get_rdf_url_template = "http://glytoucan.org/glyspace/service/glycans/{accession}/rdf"


def get(accession):
    '''
    Get the structure for `accession` from :title-reference:`GlySpace`.
    '''
    if check_cache(accession):
        return cache[accession].structure
    r = requests.get(get_url_template.format(accession=accession))
    r.raise_for_status()
    condensed = r.json["structure"]
    return glycoct.loads(condensed).next()


#: GlySpace supplies a detailed schema link which allows `lxml` to easily pull out
#: more than just the GlycoCT string. To download a more informative record, use :func:`get_record`

def get_record(accession):
    '''
    Get the complete record for `accession` from :title-reference:`GlySpace`.

    Returns
    -------
    |Glycan|
    '''
    if check_cache(accession):
        return cache[accession]
    r = requests.get(get_url_template.format(accession=accession))
    r.raise_for_status()
    return glycan_record_from_json(r.json())


def glycan_record_from_json(response, record_type=GlycanRecord):
    structure = glycoct.loads(response['structure']).next()
    accession = response['accessionNumber']
    record_id = response['glycanId']
    motifs = response['motifs']
    dbxref = [DatabaseEntry("glySpace", record_id), DatabaseEntry("glySpace", accession)]
    return GlycanRecord(structure, id=record_id, dbxref=dbxref)
