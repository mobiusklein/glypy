import sqlite3
from collections import Counter, Iterable
from functools import partial

from pygly2.utils import pickle, classproperty, make_struct
from pygly2.io.nomenclature import identity


Taxon = make_struct("Taxon", ("tax_id", "name", 'entries'))
Aglyca = make_struct("Aglyca", ("name", "reducing", 'entries'))
DatabaseEntry = make_struct("DatabaseEntry", ("database", "id"))
Motif = make_struct("Motif", ("name", "id", "motif_class"))


def _resolve_metadata_mro(cls):
    '''
    Given a class with :attr:`__metadata_map` mangled attributes
    from its class hierarchy, extract in descending order each
    `dict`, overwriting old settings as it approaches the most
    recently descended class.

    Parameters
    ----------
    cls: type
        The type to attempt to extract metadata mappings for
        along the MRO

    Returns
    -------
    dict:
        The metadata mapping describing the entire class hierarchy along
        `cls`'s MRO.
    '''
    attrs = []
    for attr in dir(cls):
        if "metadata_map" in attr:
            attrs.append(attr)
    mapping = {attr[1:].split("__metadata_map")[0]: attr for attr in attrs}
    meta_map = {}
    for typ in cls.mro()[::-1][1:]:
        metadata = mapping.get(typ.__name__, "__metadata_map")
        meta_map.update(getattr(cls, metadata, {}))
    return meta_map


def metadata(name, dtype, transform):
    '''
    Decorator for adding metadata to a record class

    Parameters
    ----------
    name: str
        Name of the new metadata field
    dtype: str
        The SQL data type to encode the column as
    transform: function
        The function to extract the value of the metadata from a record

    Returns
    -------
    function:
        Decorator that will call :meth:`.add_metadata` with `name`, `dtype` and
        `transform` on the decorated class after instantiation
    '''
    def func(cls):
        cls.add_metadata(name, dtype, transform)
        return cls
    return func


class QueryMethod(object):
    def __init__(self, func):
        self.func = func

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)


def querymethod(func):
    return classmethod(QueryMethod(func))


class RecordMethodsMeta(type):
    def __new__(cls, name, bases, body):
        querymethods = {}
        for fname, func in body.items():
            if isinstance(func, classmethod) and isinstance(func.__func__, QueryMethod):
                querymethods[fname] = func
        body['_query_methods'.format(name)] = querymethods
        typ = type.__new__(cls, name, bases, body)
        return typ


class GlycanRecordBase(object):
    __metaclass__ = RecordMethodsMeta
    '''
    Defines the base class for SQL serialize-able records carrying
    glycan structure data and metadata. Includes tools for extending the
    SQL schema describing the structure to make new information query-able.

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
    drop table if exists {table_name};
    create table {table_name}(
        glycan_id integer unique primary key not null,
        mass float not null,
        structure text not null/*rest*/
    );
    '''

    #: The storage for base-class specific
    #: metadata mappings. Add metadata here to
    #: include in all GlycanRecordBase'd objects
    __metadata_map = {}

    __query_methods = {}

    @classmethod
    def add_metadata(cls, name, dtype, transform):
        '''
        Function-based approach to modifying the class-specific
        :attr:`.__metadata_map` attribute.

        Parameters
        ----------
        name: str
            Name of the new metadata field
        dtype: str
            The SQL data type to encode the column as
        transform: function
            The function to extract the value of the metadata from a record
        '''
        cls.__metadata_map[name] = (dtype, transform)

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
        meta_map = dict(cls.__metadata_map)
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
        yield '''create index if not exists mass_index on {table_name}(mass desc);'''.format(
            table_name=cls.table_name)

    def __init__(self, structure, motifs=None, dbxref=None, aglycones=None, taxa=None, **kwargs):
        self.structure = structure
        self.motifs = motifs or []
        self.dbxref = dbxref or []
        self.aglycones = aglycones or []
        self.taxa = taxa or []
        self.id = kwargs.get('id')

    def mass(self, average=False, charge=0, mass_data=None):
        '''
        Calculates the mass of :attr:`structure`.

        See Also
        --------
        :meth:`pygly2.structure.glycan.Glycan.mass`
        '''
        return self.structure.mass(average=average, charge=charge, mass_data=mass_data)

    def __repr__(self):  # pragma: no cover
        rep = '<{type} {id} {mass}>\n{glycoct}'.format(
                id=(self.id or ''), mass=self.mass(), glycoct=self.structure.to_glycoct(),
                type=self.__class__.__name__
            )
        return rep

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
            Mapping of inherited metadata properties to include in the record

        Yields
        ------
        str:
            The SQL insert statement adding this record to the database
        '''
        inherits = dict(inherits or {})
        inherits.update(inherits)

        template = '''insert into {table_name} (glycan_id, mass, structure /*rest*/)
         values ({id}, {mass}, "{structure}" /*values*/);'''
        ext_names = ', '.join(inherits)
        if len(ext_names) > 0:
            ext_names = ', ' + ext_names
        ext_values = ', '.join(["{}".format(v) for k, v in self.collect_ext_data().items()])
        if len(ext_values) > 0:
            ext_values = ', ' + ext_values
        if id is not None:
            self.id = id

        template = template.replace("/*rest*/", ext_names).replace("/*values*/", ext_values)
        values = {}
        values['id'] = self.id
        values['mass'] = self.structure.mass(**(mass_params or {}))
        values['structure'] = pickle.dumps(self)
        values['table_name'] = self.__table_name
        yield template.format(**values)

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

    @classproperty
    def query_methods(cls):
        return cls._query_methods

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
        return pickle.loads(str(row["structure"]))

    def collect_ext_data(self, inherits=None):
        '''
        Apply each metadata mapping transform sequentially, storing each result
        in a |dict| object, returning the collection of extension data.

        '''
        meta_map = dict(inherits or {})
        meta_map.update(dict(self.__metadata_map))
        data = {}
        for name, type_transformer in meta_map.items():
            ext_type, transformer = type_transformer
            data[name] = transformer(self)
        return data

    def __eq__(self, other):
        return self.structure == other.structure

    def __ne__(self, other):
        return not self == other


def extract_composition(record, max_size=80):
    '''
    Given a :class:`.GlycanRecord`, translate its :attr:`.monosaccharides` property
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
    composition_list = ["{}:{}".format(name, count) for name, count in record.monosaccharides.items()]
    if sum(map(len, composition_list)) + len(composition_list) > max_size:
        raise ValueError(
            "The resulting composition string is larger than {} characters.".format(max_size))
    return '\"' + ' '.join(composition_list) + '\"'


def query_composition(prefix=None, **kwargs):
    if prefix is None:
        prefix = ''
    else:
        prefix += '.'
    col_name = prefix + "composition"
    composition_list = ["{} like %{}:{}%".format(col_name, name, count) for name, count in kwargs.items()]
    return ' and '.join(composition_list)


@metadata("composition", "varchar(80)", extract_composition)
class GlycanRecord(GlycanRecordBase):

    __metadata_map = {}

    extract_composition = staticmethod(extract_composition)
    query_composition = staticmethod(query_composition)

    @querymethod
    def find_like_composition(cls, cursor, select_stmt="select * from {table_name}", prefix=None, record=None):
        stmt = select_stmt + " " + cls.query_composition(prefix, **record.monosaccharides) + ";"
        return cursor.execute(stmt)

    @property
    def monosaccharides(self):
        return Counter(map(identity.identify, self.structure))

    @classmethod
    def sql_schema(cls, *args, **kwargs):
        meta_map = dict(cls.__metadata_map)
        meta_map.update(kwargs.pop("inherits", {}))
        return super(GlycanRecord, cls).sql_schema(inherits=_resolve_metadata_mro(cls))

    def to_sql(self, *args, **kwargs):
        inherits = _resolve_metadata_mro(self.__class__)
        return super(GlycanRecord, self).to_sql(*args, inherits=inherits, **kwargs)

    def collect_ext_data(self):
        inherits = _resolve_metadata_mro(self.__class__)
        data = super(GlycanRecord, self).collect_ext_data(inherits=inherits)
        return data


def make_rectype(recname="GlycanRecordType", **kwargs):
    '''
    A helper method for creating new specializations of the
    GlycanRecord class which only add new metadata mappings
    that the user does not want to add to all instances of the
    GlycanRecord type.

    Parameters
    ----------
    recname: str
        The name of the new subclass. Defaults to "GlycanRecordType".
        Does not need to be unique.
    **kwargs: dict
        Passing arbitrary key-value pairs of names and tuples of data type strings
        and transform functions to be added to the metadata mapping
    '''
    rectype = type(recname, (GlycanRecord,), {"_{}__metadata_map".format(recname): kwargs})
    return rectype


class RecordDatabase(object):
    '''
    A wrapper around an Sqlite3 database for storing and searching GlycanRecord
    objects.

    This class defines a handful general data access methods as well as the ability
    to directly write SQL queries against the database.

    Attributes
    ----------
    connection_string: |str|
        The path to the Sqlite database file, or the ":memory:" special
        keyword defining the database to be held directly in memory.
    record_type: :class:`type`
        The class type of the records assumed to be stored in this database. The stored
        class type is used for inferring the table schema, table name, and :meth:`GlycanRecord.from_sql`
        function.
    '''
    def __init__(self, connection_string=":memory:", record_type=GlycanRecord, records=None):
        '''
        Calls :meth:`.apply_schema`. If `records` is not |None|, calls :meth:`.apply_indices`.

        Parameters
        ----------
        connection_string: str
            The path to the Sqlite database file, or the ":memory:" special
            keyword defining the database to be held directly in memory.
        record_type: :class:`type`
            The class type of the records assumed to be stored in this database. The stored
            class type is used for inferring the table schema, table name, and :meth:`GlycanRecord.from_sql`
            function. Defaults to :class:`GlycanRecord`
        records: list
            A list of `record_type` records to insert immediately on table creation. If not provided,
            defaults to |None| and no records are added. If records are provided, they are inserted
            with :meth:`.load_data` and afterwards :meth:`.apply_indices` is called.
        '''
        self.connection_string = connection_string
        self.connection = sqlite3.connect(connection_string)
        self.connection.row_factory = sqlite3.Row
        self.cursor = self.connection.cursor
        self.record_type = record_type
        self._id = 0
        self.apply_schema()
        if records is not None:
            self.load_data(records)
            self.apply_indices()

    def apply_schema(self):
        '''
        Executes each SQL block yielded by :attr:`.record_type`'s :meth:`.sql_schema` class method.
        Commits all pending changes.

        Called during initialization.
        '''
        self.connection.executescript('\n'.join(self.record_type.sql_schema()))
        self.connection.commit()

    def apply_indices(self):
        '''
        Executes each SQL block yielded by :attr:`.record_type`'s :meth:`.add_index` class method.
        Commits all pending changes.

        May be called during initialization if data was added.
        '''
        for ix_stmt in self.record_type.add_index():
            self.connection.executescript(ix_stmt)
        self.connection.commit()

    def load_data(self, record_list):
        '''
        Given an iterable of :attr:`.record_type` objects,
        assign each a primary key value and insert them into the
        database.

        Commits all pending changes after all data is loaded.
        '''
        if not isinstance(record_list, Iterable):
            record_list = [record_list]
        for record in record_list:
            self._id += 1
            record.id = self._id
            for stmt in record.to_sql():
                try:
                    self.connection.execute(stmt)
                except:
                    print(stmt)
                    raise
        self.commit()

    def create(self, structure, *args, **kwargs):
        record = self.record_type(structure=structure, *args, **kwargs)
        self.load_data([record])

    def __getitem__(self, keys):
        results = []
        if isinstance(keys, int):
            keys = str(tuple([keys])).replace(',', '')
        elif isinstance(keys, slice):
            keys = tuple(range(keys.start or 1, keys.stop or self._id, keys.step or 1))
            if len(keys) == 1:
                keys = str(keys).replace(",", '')
        else:
            keys = tuple(keys)
        stmt = "select * from {table_name} where glycan_id in {keys};".format(
            table_name=self.record_type.table_name, keys=keys)
        for record in self.execute(stmt):
            results.append(self.record_type.from_sql(record, self))
        return results if len(results) > 1 else results[0]

    def __iter__(self):
        for row in self.execute(self.stn("select * from {table_name};")):
            yield self.record_type.from_sql(row, self)

    def sub_table_name(self, string, key='table_name'):
        return string.format(**{key: self.record_type.table_name})

    stn = sub_table_name

    def execute(self, query, *args, **kwargs):
        return self.connection.execute(self.stn(query), *args, **kwargs)

    def executemany(self, query_iter, *args, **kwargs):
        return self.connection.executemany(query_iter, *args, **kwargs)

    def executescript(self, script, *args, **kwargs):
        return self.connection.executescript(script, *args, **kwargs)

    def commit(self):
        self.connection.commit()

    def rollback(self):
        self.connection.rollback()

    def _find_boundaries(self, mass, tolerance):
        spread = mass * tolerance
        return (mass - spread, mass + spread)

    def ppm_match_tolerance_search(self, mass, tolerance, target_table=None,
                                   mass_shift=0):
        target_table = target_table or self.record_type.table_name
        lower, upper = self._find_boundaries(mass + mass_shift, tolerance)
        print(lower, upper)
        results = self.execute("select * from {table_name}\
         where mass between {lower} and {upper};".format(
            lower=lower, upper=upper, table_name=target_table))
        for result in results:
            yield self.record_type.from_sql(result)
