import os
import sqlite3
import logging
import functools
try:
    from collections import Counter
    from collections.abc import Iterable, Callable
except ImportError:
    from collections import Counter, Iterable, Callable


import glypy

from glypy.utils import pickle, classproperty, make_struct
from glypy.io.nomenclature.identity import naive_name_monosaccharide
from glypy.algorithms import subtree_search

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

Taxon = make_struct("Taxon", ("tax_id", "name", 'entries'))
Aglyca = make_struct("Aglyca", ("name", "reducing", 'entries'))
DatabaseEntry = make_struct("DatabaseEntry", ("database", "id"))
Motif = make_struct("Motif", ("name", "id", "motif_class"))


def _resolve_column_data_mro(cls):
    '''
    Given a class with :attr:`__column_data_map` mangled attributes
    from its class hierarchy, extract in descending order each
    `dict`, overwriting old settings as it approaches the most
    recently descended class.

    Parameters
    ----------
    cls: type
        The type to attempt to extract column_data mappings for
        along the MRO

    Returns
    -------
    dict:
        The column_data mapping describing the entire class hierarchy along
        `cls`'s MRO.
    '''
    attrs = []
    for attr in dir(cls):
        if "column_data_map" in attr:
            attrs.append(attr)
    mapping = {attr[1:].split("__column_data_map")[0]: attr for attr in attrs}
    meta_map = {}
    for typ in cls.mro()[::-1][1:]:
        column_data = mapping.get(typ.__name__, "__column_data_map")
        meta_map.update(getattr(cls, column_data, {}))
    return meta_map


def _extract_querymethods(cls):
    methods = {}
    for name, value in cls.__dict__.items():
        if isinstance(value, QueryMethod):
            methods[name] = value
    return methods


def _resolve_querymethods_mro(cls):
    methods = {}
    for typ in cls.__mro__[::-1]:
        methods.update(_extract_querymethods(typ))
    return methods


def column_data(name, dtype, transform):
    '''
    Decorator for adding a new column to the SQL table mapped to a record class

    Parameters
    ----------
    name: str
        Name of the new column
    dtype: str
        The SQL data type to encode the column as
    transform: function
        The function to extract the value of the column from a record

    Returns
    -------
    function:
        Decorator that will call :meth:`.add_column_data` with `name`, `dtype` and
        `transform` on the decorated class after instantiation
    '''
    def func(cls):
        cls.add_column_data(name, dtype, transform)
        return cls
    return func


class QueryMethod(object):
    def __init__(self, func=None):
        self.func = func

    def bind(self, obj):
        """Binds database parameters to the wrapped function

        Parameters
        ----------
        obj: RecordDatabase

        Returns
        -------
        function: The partially parameterized function bound to `obj`
        """
        @functools.wraps(self.func)
        def forward_binding(*args, **kwargs):
            return self.func(obj.record_type, obj, *args, **kwargs)
        return forward_binding

    def __call__(self, record_type, conn, *args, **kwargs):
            return self.func(record_type, conn, *args, **kwargs)


def querymethod(func):
    """Decorator for creating patching methods onto databases

    Parameters
    ----------
    func : method
        Method to be bound

    Returns
    -------
    QueryMethod
    """
    return QueryMethod(func)


class GlycanRecordBase(object):
    '''
    Defines the base class for SQL serialize-able records carrying
    glycan structure data and column_data. Includes tools for extending the
    SQL schema describing the structure to make new information query-able.

    .. warning::
        This is a base class. All user-defined record classes should be based on :class:`.GlycanRecord`

    Intended for storage in a :class:`RecordDatabase` instance, a wrapper around
    an Sqlite3 database file

    Attributes
    ----------
    structure: |Glycan|
        The |Glycan| structure being described
    motifs: |dict|
        A collection of motif names to motif locations
    dbxref: |dict|
        A collection of database names to database IDs for this structure
    aglycones: |dict|
        A collection of aglycone names to structures
    taxa: |dict|
        A collection of taxon ids to taxon objects
    id: |int|
        A unique numeric identifier used as a primary key

    The basic table schema includes a primary key, `glycan_id`, mapping to :attr:`id`.
    Additionally, it includes the mass calculated at the time of serialization under
    the `mass` column, and the entire record structure is pickled and stored under the
    `structure` column.

    The translation to SQL values is carried out by :meth:`.to_sql`, and is restored from
    a query row by :meth:`.from_sql`.
    '''

    #: Default table name used
    __table_name = "GlycanRecord"

    #: The default table schema. Additional
    #: items are added on, replacing /*rest*/
    __table_schema__ = '''
    DROP TABLE IF EXISTS {table_name};
    CREATE TABLE {table_name}(
        glycan_id INTEGER UNIQUE PRIMARY KEY NOT NULL,
        mass float NOT NULL,
        structure TEXT NOT NULL/*rest*/
    );
    '''

    #: The storage for base-class specific
    #: column_data mappings. Add column_data here to
    #: include in all GlycanRecordBase'd objects
    __column_data_map = {}

    @classmethod
    def add_column_data(cls, name, dtype, transform):
        '''
        Function-based approach to modifying the class-specific
        :attr:`.__column_data_map` attribute. Must use :func:`getattr`
        to prevent class-private name mangling

        Parameters
        ----------
        name: str
            Name of the new column_data field
        dtype: str
            The SQL data type to encode the column as
        transform: function
            The function to extract the value of the column_data from a record
        '''
        try:
            getattr(cls, "_{}__column_data_map".format(cls.__name__))[name] = (dtype, transform)
        except:
            setattr(cls, "_{}__column_data_map".format(cls.__name__), {name: (dtype, transform)})

    @classmethod
    def sql_schema(cls, *args, **kwargs):
        '''
        Generate the base table schema for the main GlycanRecord table.

        Parameters
        ----------
        inherits: dict
            Metadata extensions to include in the table schema

        Yields
        ------
        str:
            An SQL script block describing the GlycanRecord table. Descendants
            may yield additional statements.
        '''
        ext_stmts = []
        meta_map = dict(cls.__column_data_map)
        meta_map.update(kwargs.get('inherits', {}))
        for name, type_transformer in meta_map.items():
            template = "{name} {type}\n".format(name=name, type=type_transformer[0])
            ext_stmts.append(template)
        ext_def = ", ".join(ext_stmts)
        if len(ext_def) > 0:
            ext_def = ", " + ext_def
        schema = cls.__table_schema__.format(table_name=cls.table_name)
        schema = schema.replace("/*rest*/", ext_def)
        yield schema

    @classmethod
    def add_index(cls, *args, **kwargs):
        '''
        Generate the base table's indices for fast search

        Yields
        ------
        str:
            The SQL script block describing the mass_index of the GlycanRecord table
        '''
        yield '''CREATE INDEX IF NOT EXISTS mass_index ON {table_name}(mass DESC);'''.format(
            table_name=cls.table_name)
        for index in cls._collect_ext_indices(kwargs.get("inherits", {})):
            yield index

    def __init__(self, structure, motifs=None, dbxref=None, aglycones=None, taxa=None, **kwargs):
        self.structure = structure
        self.motifs = motifs or []
        self.dbxref = dbxref or []
        self.aglycones = aglycones or []
        self.taxa = taxa or []
        self.id = kwargs.get('id')
        self._bound_db = kwargs.get("bound_db")

    @classmethod
    def replicate(cls, record):
        instance = cls(record.structure)
        for name in dir(record):
            if "_" == name[0] or isinstance(getattr(record, name), Callable):
                continue
            try:
                setattr(instance, name, getattr(record, name))
            except:
                pass
        return instance

    def mass(self, average=False, charge=0, mass_data=None, override=None):
        '''
        Calculates the mass of :attr:`structure`. If `override` is not |None|,
        return this instead.

        See Also
        --------
        :meth:`glypy.structure.glycan.Glycan.mass`
        '''
        if override is not None:
            return override
        return self.structure.mass(average=average, charge=charge, mass_data=mass_data)

    def __repr__(self):  # pragma: no cover
        rep = '<{type} {id} {mass}>\n{glycoct}'.format(
                id=(self.id or ''), mass=self.mass(), glycoct=self.structure.serialize("glycoct"),
                type=self.__class__.__name__
            )
        return rep

    def __getstate__(self):
        state = dict(self.__dict__)
        state.pop("_bound_db", None)
        return state

    def to_sql(self, id=None, mass_params=None, inherits=None):
        '''
        Translates the :class:`GlycanRecord` instance into SQL.

        Parameters
        ----------
        id: int
            The primary key to use, overwriting :attr:`id` if present. Optional
        mass_params: tuple
            Parameters to pass to :meth:`.mass`. The output is stored
            in the SQL record as the `mass` value
        inherits: dict
            Mapping of inherited column_data properties to include in the record

        Yields
        ------
        str:
            The SQL insert statement adding this record to the database
        '''
        inherits = dict(inherits or {})
        inherits.update(inherits)

        template = '''INSERT INTO {table_name} (glycan_id, mass, structure /*rest*/)
         VALUES ({id}, {mass}, "{structure}" /*values*/);'''
        ext_names = ', '.join(inherits)
        if len(ext_names) > 0:
            ext_names = ', ' + ext_names
        ext_values = ', '.join(["{}".format(v) for k, v in self._collect_ext_data().items()])
        if len(ext_values) > 0:
            ext_values = ', ' + ext_values
        if id is not None:
            self.id = id

        template = template.replace("/*rest*/", ext_names).replace("/*values*/", ext_values)
        values = {}
        values['id'] = self.id
        values['mass'] = self.mass(**(mass_params or {}))
        _bound_db = getattr(self, "_bound_db", None)
        self._bound_db = None
        values['structure'] = pickle.dumps(self)
        self._bound_db = _bound_db
        values['table_name'] = self.__table_name
        yield template.format(**values)

    def to_update_sql(self, mass_params=None, inherits=None, *args, **kwargs):
        '''
        Generates SQL for use with ``UPDATE {table_name} set ... where glycan_id = {id};``.

        Called by :meth:`update`
        '''

        inherits = dict(inherits or {})
        inherits.update(inherits)

        template = '''UPDATE {table_name} SET mass = {mass},
         structure = "{structure}" /*rest*/ WHERE glycan_id = {id};'''

        ext_names = list(inherits)
        ext_values = ["{}".format(v) for k, v in self._collect_ext_data().items()]
        ext_parts = ', '.join(["{} = {}".format(name, value) for name, value in zip(ext_names, ext_values)])
        if len(ext_parts) > 0:
            ext_parts = ", " + ext_parts

        template = template.replace("/*rest*/", ext_parts)
        values = {}
        values['id'] = self.id
        values['mass'] = self.mass(**(mass_params or {}))
        values['structure'] = pickle.dumps(self)
        values['table_name'] = self.__table_name

        yield template.format(**values)

    def update(self, mass_params=None, inherits=None, commit=True, *args, **kwargs):
        """Execute SQL ``UPDATE`` instructions, writing this object's values back to the
        database it was last extracted from.

        Parameters
        ----------
        mass_params: dict
            Parameters passed to :meth:`mass`

        Raises
        ------
        ValueError:
            If the record is not bound to a database

        """
        if self._bound_db is None:
            raise ValueError("Cannot commit an unbound record")
        cur = self._bound_db.cursor()
        for stmt in self.to_update_sql(mass_params=mass_params, inherits=inherits):
            cur.execute(stmt)
        if commit:
            cur.connection.commit()

    @classproperty
    def table_name(cls):
        '''
        A property to get the class-specific table name. Must be used
        to prevent :class:`.RecordDatabase` from mangling references to
        its :attr:`record_type`'s internal table name field on request.
        '''
        return cls.__table_name

    @table_name.setter
    def table_name(cls, value):
        '''
        Setter for :attr:`.table_name`
        '''
        cls.__table_name = value

    @classmethod
    def from_sql(cls, row, *args, **kwargs):
        '''
        Translate a Row object from sqlite3 into a GlycanRecord object

        Parameters
        ----------
        row: sqlite3.Row
            A dict-like object containing the pickled value of the record in
            the `structure` field

        Returns
        -------
        GlycanRecord:
            The unpickled :class:`GlycanRecord` object. Sub-classes may perform
            more complex operations like decompressing or joining other tables in
            the database.
        '''
        record = pickle.loads(bytes(row["structure"]))
        record._bound_db = kwargs.get("database")
        return record

    def _collect_ext_data(self, inherits=None):
        '''
        Apply each column_data mapping transform sequentially, storing each result
        in a |dict| object for generating SQL.
        '''
        meta_map = dict(inherits or {})
        meta_map.update(dict(self.__column_data_map))
        data = {}
        for name, type_transformer in meta_map.items():
            ext_type, transformer = type_transformer
            data[name] = transformer(self)
        return data

    @classmethod
    def _collect_ext_indices(cls, inherits=None):
        meta_map = dict(inherits or {})
        meta_map.update(dict(cls.__column_data_map))
        for name in meta_map:
            yield "CREATE INDEX IF NOT EXISTS {column_name}_{table_name}_index ON {table_name}({column_name});".format(
                column_name=name, table_name=cls.table_name)

    def __eq__(self, other):
        '''
        Equality testing is done between :attr:`structure`
        '''
        try:
            return self.structure == other.structure
        except Exception as e:  # pragma: no cover
            logger.exception(exc_info=e)
            return False

    def __ne__(self, other):  # pragma: no cover
        return not self == other

    def __root__(self):  # pragma: no cover
        return self.structure.root

    def __tree__(self):  # pragma: no cover
        return self.structure


def extract_composition(record, max_size=120):
    '''
    Given a :class:`GlycanRecord`, translate its :attr:`.monosaccharides` property
    into a string suitable for denormalized querying.

    Transforms the resulting, e.g. Counter({u'GlcNA': 6, u'Gal': 4, u'aMan': 2, u'Fuc': 1, u'Man': 1})
    into the string "Gal:4 aMan:2 Fuc:1 GlcNA:6 Man:1" which could be partially matched in
    queries using SQL's LIKE operator.

    Parameters
    ----------
    record: GlycanRecord
        The record to serialize :attr:`.monosaccharides` for
    max_size: int
        The maximum size of the resulting string allowed under the target
        SQL schema

    Returns
    -------
    str:
        The string representation of `record.monosaccharides`.
    '''
    composition_list = ["{}:{}".format(name, count) for name, count in sorted(record.monosaccharides.items())]
    if sum(map(len, composition_list)) + len(composition_list) > max_size:
        raise ValueError(
            "The resulting composition string is larger than {} characters.".format(max_size))
    return '\"' + ' '.join(composition_list) + '\"'


def _query_composition(prefix=None, **kwargs):
    if prefix is None:
        prefix = ''
    else:
        prefix += '.'
    col_name = prefix + "composition"
    composition_list = ["{} LIKE '%{}:{}%'".format(col_name, name, count) for name, count in kwargs.items()]
    return ' AND '.join(composition_list)


def is_n_glycan(record):
    '''
    A predicate testing if the :title:`N-linked Glycan` core motif is present in `record.structure`.

    Returns
    -------
    int:
        0 if |False|, 1 if |True|. Sqlite doesn't have a dedicated |bool| data type.

    See Also
    --------
    :func:`glypy.algorithms.subtree_search.subtree_of`
    '''
    return int(subtree_search.subtree_of(
        n_glycan_core, record.structure, exact=False) == 1)
n_glycan_core = glypy.glycans["N-Linked Core"]


@column_data("composition", "VARCHAR(120)", extract_composition)
@column_data("is_n_glycan", "BOOLEAN", is_n_glycan)
@column_data("glycoct", "TEXT", lambda x: "\"{}\"".format(str(x.structure)))
class GlycanRecord(GlycanRecordBase):
    '''
    An extension of :class:`GlycanRecordBase` to add additional features and better support for extension
    by both metaprogramming and inheritance.
    '''

    __column_data_map = {}

    @querymethod
    def query_like_composition(cls, conn, record=None, prefix=None):
        stmt = "SELECT * FROM {table_name} WHERE " + _query_composition(prefix, **record.monosaccharides) + ";"
        for result in conn.from_sql(conn.execute(stmt)):
            yield result

    @classmethod
    def add_index(cls, *args, **kwargs):
        '''
        Generate the base table's indices for fast search

        Yields
        ------
        str:
            The SQL script block describing the mass_index of the GlycanRecord table
        '''
        for stmt in super(GlycanRecord, cls).add_index(*args, **kwargs):
            yield stmt
        for stmt in cls._collect_ext_indices(kwargs.get("inherits", {})):
            yield stmt

    @property
    def monosaccharides(self):
        '''
        Returns a mapping of the counts of monosaccharides found in :attr:`structure`. Generic names are found
        using :func:`naive_name_monosaccharide`.

        .. note::
            This property is mapped to the the database column ``composition`` by :func:`extract_composition`.


        See Also
        --------
        :func:`extract_composition`
        :func:`naive_name_monosaccharide`
        '''
        res = Counter(map(naive_name_monosaccharide, self.structure))
        res.pop(None, None)
        return res

    @property
    def composition(self):
        return ''.join("{}{}".format(name, number)
                       for name, number in sorted(self.monosaccharides.items()))

    @property
    def is_n_glycan(self):
        '''
        Returns |True| if :attr:`structure` has the :title:`N-linked Glycan` core motif.

        .. note:: This property is mapped to the the database column ``is_n_glycan`` by :func:`is_n_glycan`
        '''
        return bool(is_n_glycan(self))

    @classmethod
    def sql_schema(cls, *args, **kwargs):
        meta_map = dict(cls.__column_data_map)
        meta_map.update(kwargs.pop("inherits", {}))
        return super(GlycanRecord, cls).sql_schema(inherits=_resolve_column_data_mro(cls))

    def to_sql(self, *args, **kwargs):
        inherits = _resolve_column_data_mro(self.__class__)
        return super(GlycanRecord, self).to_sql(*args, inherits=inherits, **kwargs)

    def to_update_sql(self, *args, **kwargs):
        inherits = _resolve_column_data_mro(self.__class__)
        inherits = inherits or _resolve_column_data_mro(self.__class__)
        kwargs['inherits'] = inherits
        return super(GlycanRecord, self).to_update_sql(*args, **kwargs)

    def update(self, mass_params=None, inherits=None, commit=True, *args, **kwargs):
        inherits = inherits or _resolve_column_data_mro(self.__class__)
        kwargs['inherits'] = inherits
        super(GlycanRecord, self).update(mass_params=mass_params, *args, **kwargs)

    def _collect_ext_data(self):
        inherits = _resolve_column_data_mro(self.__class__)
        data = super(GlycanRecord, self)._collect_ext_data(inherits=inherits)
        return data

    @classmethod
    def _collect_ext_indices(cls, inherits=None):
        inherits = _resolve_column_data_mro(cls)
        data = super(GlycanRecord, cls)._collect_ext_indices(inherits=inherits)
        return data


class GlycanRecordWithTaxon(GlycanRecord):
    __taxa_table_schema__ = '''
    DROP TABLE IF EXISTS RecordTaxonomy;
    CREATE TABLE RecordTaxonomy(
        glycan_id INTEGER NOT NULL,
        taxon_id INTEGER NOT NULL
    );
    '''

    @classmethod
    def sql_schema(cls, *args, **kwargs):
        for line in super(GlycanRecordWithTaxon, cls).sql_schema(*args, **kwargs):
            yield line
        yield cls.__taxa_table_schema__

    @classmethod
    def add_index(cls, *args, **kwargs):
        for line in super(GlycanRecordWithTaxon, cls).add_index(*args, **kwargs):
            yield line
        yield "CREATE INDEX IF NOT EXISTS TaxonomyIndex ON RecordTaxonomy(taxon_id);"
        yield "CREATE INDEX IF NOT EXISTS TaxonomyIndex2 ON RecordTaxonomy(glycan_id);"

    def to_sql(self, *args, **kwargs):
        for line in super(GlycanRecordWithTaxon, self).to_sql(*args, **kwargs):
            yield line
        for taxon in self.taxa:
            yield "INSERT OR REPLACE INTO RecordTaxonomy (glycan_id, taxon_id) VALUES ({}, {});".format(
                self.id, int(taxon.tax_id))

    @querymethod
    def query_by_taxon_id(cls, conn, taxon_ids):
        # Passed an iterable of taxa to search
        try:
            taxon_ids = tuple(taxon_ids)
            taxon_ids_token = str(taxon_ids)
            if len(taxon_ids) == 1:
                taxon_ids_token = taxon_ids_token.replace(',', '')

            return conn.from_sql(conn.execute(
                """SELECT {table_name}.* FROM {table_name} JOIN RecordTaxonomy taxa ON
                   taxa.glycan_id = {table_name}.glycan_id WHERE taxa.taxon_id IN [!SUB!];
                """.replace('[!SUB!]', taxon_ids_token)))
        # Passed a single taxon
        except:
            return conn.from_sql(conn.execute(
                """SELECT {table_name}.* FROM {table_name} JOIN RecordTaxonomy taxa ON
                   taxa.glycan_id = {table_name}.glycan_id WHERE taxa.taxon_id = ?;""",
                (int(taxon_ids),)))


class RecordDatabase(object):
    '''
    A wrapper around an Sqlite3 database for storing and searching GlycanRecord
    objects.

    This class defines a handful general data access methods as well as the ability
    to directly write SQL queries against the database.

    Calls :meth:`.apply_schema`. If ``records`` is not |None|, calls :meth:`.apply_indices`.

    ``record_type``, the stored class type is used for inferring the table schema,
    table name, and :meth:`GlycanRecord.from_sql` function.

    If ``records`` is not provided, no records are added. If records are provided, they are inserted
    with :meth:`.load_data` and afterwards :meth:`.apply_indices` is called.

    Attributes
    ----------
    connection_string: |str|
        The path to the Sqlite database file, or the ":memory:" special
        keyword defining the database to be held directly in memory.
    record_type: :class:`type`
        The class type of the records assumed to be stored in this database. The stored
        class type is used for inferring the table schema, table name, and :meth:`GlycanRecord.from_sql`
        function.

    Parameters
    ----------
    connection_string: str
        The path to the Sqlite database file, or the ":memory:" special
        keyword defining the database to be held directly in memory.
    record_type: type
        The class type of the records assumed to be stored in this database. Defaults to :class:`GlycanRecord`
    records: list
        A list of `record_type` records to insert immediately on table creation.
    '''
    def __init__(self, connection_string=":memory:", record_type=GlycanRecord, records=None, flag='c'):

        created_new = False
        if connection_string == ":memory:" or not os.path.exists(connection_string):
            created_new = True

        if flag == "c":
            pass
            # If 'c', create if not present, but don't clear old data
        elif flag == "w":
            created_new = True
            # If 'w', clear the table before taking any operations

        self.connection_string = connection_string
        self.connection = sqlite3.connect(connection_string)
        self.connection.row_factory = sqlite3.Row
        self.cursor = self.connection.cursor

        # Check to see if the record type matches what is already
        # stored in the database.
        try:
            if record_type is None:
                self.record_type = GlycanRecord
            else:
                self.record_type = record_type
            rec = next(iter(self))
            if not type(rec) == record_type:
                if record_type is not None:
                    logger.warn("Record class {0} is not the same as {1}. Using {0}".format(type(rec), record_type))
                    print("Record class {0} is not the same as {1}. Using {0}".format(type(rec), record_type))
                record_type = type(rec)
        except (StopIteration, sqlite3.OperationalError):
            # If record_type is None and no types were inferred, use GlycanRecord
            if record_type is None:
                record_type = GlycanRecord

        self.record_type = record_type
        self._id = 0

        if records is not None:
            self.apply_schema()
            self.load_data(records)
            self.apply_indices()
        elif created_new:
            self.apply_schema()
        else:
            self._id = len(self)
        try:
            self._patch_querymethods()
        except Exception as e:
            logger.error(exc_info=e)

    def apply_schema(self):
        '''
        Executes each SQL block yielded by :attr:`.record_type`'s :meth:`.sql_schema` class method.
        Commits all pending changes.

        Called during initialization if the database is newly created.

        .. danger::
            The SQL table definition statements generated may drop existing tables. Calling
            this function on an already populated table can cause data loss.
        '''
        self.executescript('''
        DROP TABLE IF EXISTS metadata;
        CREATE TABLE metadata(
            name VARCHAR(40) PRIMARY KEY,
            content TEXT
        );''')
        self.executescript('\n'.join(self.record_type.sql_schema()))
        self.commit()
        self._id = 0

    def apply_indices(self):
        '''
        Executes each SQL block yielded by :attr:`.record_type`'s :meth:`.add_index` class method.
        Commits all pending changes.

        May be called during initialization if data was added.
        '''
        for ix_stmt in self.record_type.add_index():

            self.execute(ix_stmt)
        self.commit()

    def load_data(self, record_list, commit=True, set_id=True, cast=True, **kwargs):
        '''
        Given an iterable of :attr:`.record_type` objects,
        assign each a primary key value and insert them into the
        database.

        Forwards all ``**kwargs`` to :meth:`to_sql` calls.

        Parameters
        ----------
        record_list: GlycanRecord or iterable of GlycanRecords
        commit: bool
            Whether or not to commit all changes to the database
        set_id: bool
        cast: bool
        '''
        if not isinstance(record_list, Iterable):
            record_list = [record_list]
        for record in record_list:
            if set_id:
                self._id += 1
                record.id = self._id
            if cast and not isinstance(record, self.record_type):
                record = self.record_type.replicate(record)
            for stmt in record.to_sql(**kwargs):
                try:
                    self.connection.execute(stmt)
                except:
                    print(stmt)
                    raise
        if commit:
            self.commit()

    def get_metadata(self, key=None):
        """Retrieve a value from the key-value store
        in the database's metadata table.

        If `key` is |None| all of the keys will be retrieved
        and returned as a |dict|.

        Parameters
        ----------
        key : str, optional
            The key value to retrieve.

        Returns
        -------
        any or |dict|
        """
        if key is None:
            return {
                row[0]: pickle.loads(bytes(row[1])) for row in self.execute("SELECT name, content FROM metadata")
            }
        res = [pickle.loads(bytes(row[0])) for row in self.execute("SELECT content FROM metadata WHERE name = ?", (key,))]
        if len(res) > 0:
            return res[0]
        else:
            raise KeyError("Key {} not found".format(key))

    def set_metadata(self, key, value):
        """Set a key-value pair in the database's metadata table

        Parameters
        ----------
        key : str
            Key naming value
        value : any
            Value to store in the metadata table
        """
        self.execute("INSERT OR REPLACE INTO metadata (name, content) VALUES (?, ?);", (key, pickle.dumps(value)))
        self.commit()

    def bind(self, record):
        record._bound_db = self

    def _patch_querymethods(self):
        """Patch on additional methods defined on :attr:`record_type`
        """
        for name, value in _resolve_querymethods_mro(self.record_type).items():
            if isinstance(value, QueryMethod):
                setattr(self, name, value.bind(self))

    def __len__(self):
        """The number of records in the database

        Returns
        -------
        int
        """
        res = (self.execute("SELECT count(glycan_id) FROM {table_name};").next())["count(glycan_id)"]
        return res or 0

    def create(self, structure, *args, **kwargs):
        '''
        A convenience function for creating a database entry for |Glycan| `structure`. Passes
        along all arguments to :attr:`record_type` initialization methods.

        Parameters
        ----------
        structure: :class:`Glycan`
        commit: bool
            If |True|, commit changes after adding this record. Defaults to |True|.
        '''
        commit = kwargs.pop("commit", True)
        record = self.record_type(structure=structure, *args, **kwargs)
        self.load_data([record], commit=commit)
        return record

    def __getitem__(self, keys):
        '''
        Look up records in the database by primary key. Also accepts :class:`slice` objects.

        Returns
        -------
        :class:`.record_type` or :class:`list` of :class:`.record_type`
        '''
        results = []
        key_type = int
        if isinstance(keys, slice):
            key_type = slice
        if isinstance(keys, (tuple, list, set)):
            key_type = tuple
        else:
            keys = int(keys)
        if key_type is int:
            results = list(self.from_sql(
                self.execute("SELECT * FROM {table_name} WHERE glycan_id = ?", (keys,))))
        elif key_type is slice:
            begin = keys.start or 1
            end = keys.stop or len(self)

            results = list(self.from_sql(
                self.execute(
                    "SELECT * FROM {table_name} WHERE glycan_id BETWEEN ? AND ?", begin, end)))
        elif key_type is tuple:
            group = tuple(map(int, keys))
            if len(group) == 1:
                group = str(group).replace(',', '')
            results = list(self.from_sql(
                self.execute(
                    "SELECT * FROM {table_name} WHERE glycan_id IN {}".format(
                        group, table_name=self.record_type.table_name))))

        if len(results) == 0 and key_type is int:
            raise IndexError("No record found for %r" % keys)
        list(map(self.bind, results))
        if len(results) > 1:
            return results
        elif key_type is int:
            return results[0]
        else:
            return results

    def __iter__(self):
        '''
        Iterate sequentially over each entry in the database.
        '''
        for row in self.from_sql(self.execute(self.stn("select * from {table_name};"))):
            yield row

    def __contains__(self, key):
        try:
            self[int(key)]
            return True
        except IndexError:
            return False

    def __repr__(self):  # pragma: no cover
        rep = "<RecordDatabase {} records>".format(len(self))
        return rep

    def sub_table_name(self, string, key='table_name'):
        '''
        A convenience function called to substitute in the primary table name
        into raw SQL queries. By default it looks for the token {table_name}.
        '''
        return string.format(**{key: self.record_type.table_name})

    stn = sub_table_name

    def execute(self, query, *args, **kwargs):
        '''
        A wrapper around :meth:`sqlite3.Connection.execute`. Will format
        the query string to substitute in the main table name if the {table_name}
        token is present
        '''
        return self.connection.execute(self.stn(query), *args, **kwargs)

    def executemany(self, query, param_iter, *args, **kwargs):
        '''
        A wrapper around :meth:`sqlite3.Connection.executemany`. Will format
        the query string to substitute in the main table name if the {table_name}
        token is present.
        '''
        return self.connection.executemany(self.stn(query), param_iter, *args, **kwargs)

    def executescript(self, script, *args, **kwargs):
        '''
        A wrapper around :meth:`sqlite3.Connection.executescript`.
        '''
        return self.connection.executescript(script, *args, **kwargs)

    def commit(self):
        '''
        A wrapper around :meth:`sqlite3.Connection.commit`. Writes pending changes to the database.
        '''
        self.connection.commit()

    def rollback(self):
        '''
        A wrapper around :meth:`sqlite3.Connection.rollback`. Reverses the last set of changes to the database.
        '''
        self.connection.rollback()

    def close(self):
        self.connection.close()

    def _find_boundaries(self, mass, tolerance):
        spread = mass * tolerance
        return (mass - spread, mass + spread)

    def ppm_match_tolerance_search(self, mass, tolerance, mass_shift=0):
        '''
        Rapidly search the database for entries with a recorded mass within
        ``tolerance`` parts per million mass error of ``mass``.

        :math:`[mass - (tolerance * mass), mass + (tolerance * mass)]`

        '''
        lower, upper = self._find_boundaries(mass + mass_shift, tolerance)
        results = self.execute("SELECT * FROM {table_name}\
         WHERE mass BETWEEN ? AND ?;", (lower, upper))
        for result in results:
            yield self.record_type.from_sql(result, database=self)

    def from_sql(self, rows, from_sql_fn=None):
        """Convenience function to convert `rows` into objects through `from_sql_fn`,
        by default, :meth:`self.record_type.from_sql`

        Parameters
        ----------
        rows : sqlite3.Row or an iterable of sqlite3.Row
            Collection of objects to convert
        from_sql_fn : function, optional
            Function to perform the conversion. Defaults to :meth:`self.record_type.from_sql`

        Yields
        -------
        Type returned by `from_sql_fn`
        """
        if from_sql_fn is None:
            from_sql_fn = self.record_type.from_sql
        if isinstance(rows, sqlite3.Row):
            rows = [rows]
        for row in rows:
            record = from_sql_fn(row, self)
            self.bind(record)
            yield record

#: Open a database
dbopen = RecordDatabase
