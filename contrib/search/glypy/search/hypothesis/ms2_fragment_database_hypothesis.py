import os
import logging
from glypy.algorithms import database
from .common_transforms import monoisotopic_mass
from glypy.utils import identity, pickle

from multiprocessing import Pool
from functools import partial

logger = logging.getLogger("ms2_fragment_database_hypothesis")


default_fragmentation_parameters = {
    "kind": "ABCYXZ",
    "max_cleavages": 1,
    "average": False,
    "charge": 0
}


default_mass_transform_parameters = {
    "derivatize_fn": monoisotopic_mass,
    "adduct_mass": 0,
    "adduct_number": 0
}


def mass_transform(record, derivatize_fn=identity, adduct_mass=0, adduct_number=0):
    """Apply a derivatization to `record`, and add any adducts to the intact mass

    Parameters
    ----------
    record : GlycanRecord
        Record to mutate
    derivatize_fn : function, optional
        Function to apply derivatizion and any other mass modifications present on
        the entire structure and its fragments
    adduct_mass : float, optional
        The mass shift induced by the adduct attached to the intact mass
    adduct_number : int, optional
        The number of adducts to include

    Returns
    -------
    float : intact mass of the derivatized structure
    """
    derivatize_fn(record)
    record.adduct_mass = adduct_mass
    record.adduct_number = adduct_number
    mass = record.mass() + adduct_number * adduct_mass
    return mass


def extract_fragments(record, fragmentation_parameters=None):
    """Generate all fragments for `record` using `fragmentation_parameters`

    Parameters
    ----------
    record : GlycanRecord
        The structure record to generate fragments with
    fragmentation_parameters : dict, optional
        A mapping passed as **kwargs to :meth:`Glycan.fragments`. Defaults to :data:`default_fragmentation_parameters`

    Returns
    -------
    list : list of :class:`Fragment` objects
    """
    default_fragmentation_parameters = {
        "kind": "ABCYXZ",
        "max_cleavages": 1,
        "average": False,
        "charge": 0
    }
    fragmentation_parameters = fragmentation_parameters or default_fragmentation_parameters
    fragments = {}
    record.structure = record.structure.clone().reindex()
    for frag in (record.structure.fragments("ABCXYZ", 1)):
        fragments[frag.name] = frag
    for frag in record.structure.fragments("BCYZ", 3):
        fragments[frag.name] = frag
    return list(fragments.values())


def record_handle(record, mass_transform_parameters, fragmentation_parameters):
    """Summary
    
    Parameters
    ----------
    record : TYPE
        Description
    mass_transform_parameters : TYPE
        Description
    fragmentation_parameters : TYPE
        Description
    
    Returns
    -------
    TYPE : Description
    """
    mass = mass_transform(record, **(mass_transform_parameters or {}))
    try:
        fragments = extract_fragments(record, fragmentation_parameters)
    except Exception, e:
        print "An exception, {} for record {}".format(e, record.id)
        # raise e
        return None
    record.fragments = fragments
    record.intact_mass = mass
    return record


def _chunk_iter(iter, size=50):
    results = []
    for entry in iter:
        results.append(entry)
        if len(results) == size:
            yield results
            results = []
    yield results


def _strip_bound_db_gen(record_iter):
    for record in record_iter:
        record._bound_db = 0
        yield record


def prepare_database(in_database, out_database=None, mass_transform_parameters=None,
                     fragmentation_parameters=None, n_processes=1):
    """Translate a database of structures into a new set of experimental structures, including
    theoretical fragments and modified masses.

    Parameters
    ----------
    in_database : str or RecordDatabase
        Database to translate. If passed a string, it will be treated as a connection string
    out_database : str or RecordDatabase, optional
        Storage for the resulting translated records. If passed a string, it will be treated
        as a connection string. If no value is given, a connection string will be interpolated
        from `in_database`
    mass_transform_parameters : dict, optional
        A dictionary of parameters passed by **kwargs to :func:`mass_transform`
    fragmentation_parameters : dict, optional
        A dictionary of parameters passed by **kwargs to :func:`extract_fragments`
    n_processes : int, optional
        The number of processes to run on. Defaults to 1.

    Returns
    -------
    out_database : RecordDatabase
    """
    if isinstance(in_database, str):
        in_database = database.RecordDatabase(in_database)
    if out_database is None:
        out_database_string = os.path.splitext(in_database.connection_string)[0] + ".out.db"
        out_database = database.RecordDatabase(out_database_string, record_type=in_database.record_type, flag='w')
    elif isinstance(out_database, str):
        out_database = database.RecordDatabase(out_database, record_type=in_database.record_type, flag='w')
    logger.info("%d records in input database", len(in_database))
    if n_processes == 1:
        for i, record in enumerate(in_database):
            mass = mass_transform(record, **(mass_transform_parameters or {}))
            fragments = extract_fragments(record, fragmentation_parameters)
            record.fragments = fragments
            record.intact_mass = mass
            out_database.load_data([record], commit=False, mass_params={"override": mass})
            logger.info("%d records processed", i)
            print(i)
    else:
        logger.info("Using pool with %d workers", n_processes)
        worker_pool = Pool(n_processes, maxtasksperchild=3)
        taskfn = partial(record_handle,
                         mass_transform_parameters=mass_transform_parameters,
                         fragmentation_parameters=fragmentation_parameters)
        do_work = True
        size = 30
        g = _chunk_iter(_strip_bound_db_gen(in_database), size)
        while do_work:
            job = g.next()
            if len(job) < size:
                do_work = False
            for i, record in enumerate(worker_pool.imap_unordered(taskfn, job, chunksize=1)):
                if record is None:
                    continue
                mass = record.intact_mass
                out_database.load_data([record], set_id=False, commit=False, mass_params={"override": mass})
                logger.info("Finished %d, %d records processed", record.id, i)
                # print(record.id, i)
            out_database.commit()
            logger.info("There are %d records in the output database", len(out_database))
            del job
    out_database.commit()
    out_database.apply_indices()
    return out_database


def duplicate_check(db):
    '''Check the passed iterable of |GlycanRecord| objects for topological duplicates.
    
    Parameters
    ----------
    db : TYPE
        Description
    '''
    blacklist = set()
    keepers = []
    for rec in db:
        print rec.id
        if rec.id in blacklist:
            continue
        keepers.append(rec)
        anomers = []
        for node in rec.structure:
            anomers.append(node.anomer)
            node.anomer = None
        m = rec.mass()
        isobarics = list(db.ppm_match_tolerance_search(m, 1e-5))
        print len(isobarics)
        for iso in isobarics:
            if iso.id == rec.id:
                continue

            iso_anomers = []
            for node in iso.structure:
                iso_anomers.append(node.anomer)
                node.anomer = None

            if iso.structure.topological_equality(rec.structure):
                blacklist.add(iso.id)

            for an, node in zip(iso_anomers, iso.structure):
                node.anomer = an

        for an, node in zip(anomers, rec.structure):
            node.anomer = an

    return keepers
