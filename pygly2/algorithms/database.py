import os
import sqlite3
from collections import Counter, Iterable

import pygly2
from pygly2.utils import pickle, classproperty, make_struct
from pygly2.io.nomenclature import identity
from pygly2.algorithms import subtree_search


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
    '''
    Defines the base class for SQL serialize-able records carrying
    glycan structure data and metadata. Includes tools for extending the
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
    __metaclass__ = RecordMethodsMeta
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
        :attr:`.__metadata_map` attribute. Must use :func:`getattr`
        to prevent class-private name mangling

        Parameters
        ----------
        name: str
            Name of the new metadata field
        dtype: str
            The SQL data type to encode the column as
        transform: function
            The function to extract the value of the metadata from a record
        '''
        getattr(cls, "_{}__metadata_map".format(cls.__name__))[name] = (dtype, transform)

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
        self._bound_db = kwargs.get("bound_db")

    def mass(self, average=False, charge=0, mass_data=None, override=None):
        '''
        Calculates the mass of :attr:`structure`. If `override` is not |None|,
        return this instead.

        See Also
        --------
        :meth:`pygly2.structure.glycan.Glycan.mass`
        '''
        if override is not None:
            return override
        return self.structure.mass(average=average, charge=charge, mass_data=mass_data)

    def __repr__(self):  # pragma: no cover
        rep = '<{type} {id} {mass}>\n{glycoct}'.format(
                id=(self.id or ''), mass=self.mass(), glycoct=self.structure.to_glycoct(),
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
        ext_values = ', '.join(["{}".format(v) for k, v in self._collect_ext_data().items()])
        if len(ext_values) > 0:
            ext_values = ', ' + ext_values
        if id is not None:
            self.id = id

        template = template.replace("/*rest*/", ext_names).replace("/*values*/", ext_values)
        values = {}
        values['id'] = self.id
        values['mass'] = self.structure.mass(**(mass_params or {}))
        _bound_db = self._bound_db
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

        template = '''update {table_name} set mass = {mass}, structure = "{structure}" /*rest*/ where glycan_id = {id};'''

        ext_names = list(inherits)
        ext_values = ["{}".format(v) for k, v in self._collect_ext_data().items()]
        ext_parts = ', '.join(["{} = {}".format(name, value) for name, value in zip(ext_names, ext_values)])
        if len(ext_parts) > 0:
            ext_parts = ", " + ext_parts

        template = template.replace("/*rest*/", ext_parts)
        values = {}
        values['id'] = self.id
        values['mass'] = self.structure.mass(**(mass_params or {}))
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
        record = pickle.loads(str(row["structure"]))
        record._bound_db = kwargs.get("database")
        return record

    def _collect_ext_data(self, inherits=None):
        '''
        Apply each metadata mapping transform sequentially, storing each result
        in a |dict| object for generating SQL.
        '''
        meta_map = dict(inherits or {})
        meta_map.update(dict(self.__metadata_map))
        data = {}
        for name, type_transformer in meta_map.items():
            ext_type, transformer = type_transformer
            data[name] = transformer(self)
        return data

    def __eq__(self, other):
        '''
        Equality testing is done between :attr:`structure`
        '''
        try:
            return self.structure == other.structure
        except:
            return False

    def __ne__(self, other):
        return not self == other


def extract_composition(record, max_size=120):
    '''
    Given a :class:`GlycanRecord`, translate its :attr:`.monosaccharides` property
    into a string suitable for denormalized querying.

    Transforms the resulting, e.g. Counter({u'GlcNA': 6, u'Gal': 4, u'aMan': 2, u'Fuc': 1, u'Man': 1})
    into the string "Gal:4 aMan:2 Fuc:1 GlcNA:6 Man:1" which could be partially matched in
    queries using SQL's LIKE operator.

    Parameters
    ----------
    record: :class:`GlycanRecord`
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


def naive_name_monosaccharide(monosaccharide):
    '''
    Generate a generic name for `monosaccharide`, based loosely on IUPAC
    naming schema without including information about linkage.

    The tendency for monosaccharides of superclass > 7 to have special names,
    which will be used preferentially if possible.

    Parameters
    ----------
    monosaccharide: Monosaccharide

    Returns
    -------
    str:
        A simple name based on `SuperClass`, modifications, and substituents.

    See Also
    --------
    :func:`pygly2.io.nomenclature.identity.identify`

    '''
    try:
        c = monosaccharide.clone()
        if monosaccharide.superclass.value > 7:
            return identity.identify(monosaccharide, tolerance=1)
        c.anomer = None
        return identity.identify(c)
    except identity.IdentifyException:
        try:
            c.stem = None
            c.configuration = None
            return identity.identify(c)
        except identity.IdentifyException:
            return "".join(mod.name for mod in c.modifications.values() if mod.name != 'aldi') +\
             c.superclass.name.title() + ''.join(
                ["".join(map(str.title, subst.name.split("_")))[:3] for p, subst in c.substituents()])


def is_n_glycan(record):
    '''
    A predicate testing if the :title:`N-linked Glycan` core motif is present in `record.structure`.

    Returns
    -------
    int:
        0 if |False|, 1 if |True|. Sqlite doesn't have a dedicated |bool| data type.

    See Also
    --------
    :func:`pygly2.algorithms.subtree_search.subtree_of`
    '''
    return int(subtree_search.subtree_of(
        n_glycan_core, record.structure, exact=False) == 1)
n_glycan_core = pygly2.glycans["N-Linked Core"]


@metadata("composition", "varchar(120)", extract_composition)
@metadata("is_n_glycan", "boolean", is_n_glycan)
class GlycanRecord(GlycanRecordBase):
    '''
    An extension of :class:`GlycanRecordBase` to add additional features and better support for extension
    by both metaprogramming and inheritance.
    '''

    __metadata_map = {}

    extract_composition = staticmethod(extract_composition)
    query_composition = staticmethod(query_composition)

    @querymethod
    def find_like_composition(cls, cursor, select_stmt="select * from {table_name}", prefix=None, record=None):
        stmt = select_stmt + " " + cls.query_composition(prefix, **record.monosaccharides) + ";"
        return cursor.execute(stmt)

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
        return Counter(map(naive_name_monosaccharide, self.structure))

    @property
    def is_n_glycan(self):
        '''
        Returns |True| if :attr:`structure` has the :title:`N-linked Glycan` core motif.

        .. note:: This property is mapped to the the database column ``is_n_glycan`` by :func:`is_n_glycan`
        '''
        return bool(is_n_glycan(self))

    @classmethod
    def sql_schema(cls, *args, **kwargs):
        meta_map = dict(cls.__metadata_map)
        meta_map.update(kwargs.pop("inherits", {}))
        return super(GlycanRecord, cls).sql_schema(inherits=_resolve_metadata_mro(cls))

    def to_sql(self, *args, **kwargs):
        inherits = _resolve_metadata_mro(self.__class__)
        return super(GlycanRecord, self).to_sql(*args, inherits=inherits, **kwargs)

    def to_update_sql(self, *args, **kwargs):
        inherits = _resolve_metadata_mro(self.__class__)
        inherits = inherits or _resolve_metadata_mro(self.__class__)
        kwargs['inherits'] = inherits
        return super(GlycanRecord, self).to_update_sql(*args, **kwargs)

    def update(self, mass_params=None, inherits=None, commit=True, *args, **kwargs):
        inherits = inherits or _resolve_metadata_mro(self.__class__)
        kwargs['inherits'] = inherits
        super(GlycanRecord, self).update(mass_params=mass_params, *args, **kwargs)

    def _collect_ext_data(self):
        inherits = _resolve_metadata_mro(self.__class__)
        data = super(GlycanRecord, self)._collect_ext_data(inherits=inherits)
        return data


def include_fragments(record, kind="BY", average=False, charge=0, mass_data=None):
    record.fragments = list(record.structure.fragments(
        kind=kind, average=average, charge=charge, mass_data=mass_data))


def make_rectype(recname="GlycanRecordType", **kwargs):
    '''
    Programmatically create a new ``type`` based on :class:`GlycanRecord` at
    run time.

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
    connection_string: :class:`str`
        The path to the Sqlite database file, or the ":memory:" special
        keyword defining the database to be held directly in memory.
    record_type: :class:`type`
        The class type of the records assumed to be stored in this database. Defaults to :class:`GlycanRecord`
    records: :class:`list`
        A list of `record_type` records to insert immediately on table creation.
    '''
    def __init__(self, connection_string=":memory:", record_type=GlycanRecord, records=None):

        created_new = False
        if connection_string == ":memory:" or not os.path.exists(connection_string):
            created_new = True

        self.connection_string = connection_string
        self.connection = sqlite3.connect(connection_string)
        self.connection.row_factory = sqlite3.Row
        self.cursor = self.connection.cursor
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

    def apply_schema(self):
        '''
        Executes each SQL block yielded by :attr:`.record_type`'s :meth:`.sql_schema` class method.
        Commits all pending changes.

        Called during initialization if the database is newly created.

        .. danger::
            The SQL table definition statements generated may drop existing tables. Calling
            this function on an already populated table can cause data loss.
        '''
        self.connection.executescript('\n'.join(self.record_type.sql_schema()))
        self.connection.commit()
        self._id = 0

    def apply_indices(self):
        '''
        Executes each SQL block yielded by :attr:`.record_type`'s :meth:`.add_index` class method.
        Commits all pending changes.

        May be called during initialization if data was added.
        '''
        for ix_stmt in self.record_type.add_index():
            self.connection.executescript(ix_stmt)
        self.connection.commit()

    def load_data(self, record_list, commit=True, **kwargs):
        '''
        Given an iterable of :attr:`.record_type` objects,
        assign each a primary key value and insert them into the
        database.

        Forwards all ``**kwargs`` to :meth:`to_sql` calls.

        Commits all pending changes after all data is loaded.
        '''
        if not isinstance(record_list, Iterable):
            record_list = [record_list]
        for record in record_list:
            self._id += 1
            record.id = self._id
            for stmt in record.to_sql(**kwargs):
                try:
                    self.connection.execute(stmt)
                except:
                    print(stmt)
                    raise
        if commit:
            self.commit()

    def __len__(self):
        res = (self.execute("select count(glycan_id) from {table_name};").next())["count(glycan_id)"]
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

    # TODO: Generate more efficient SQL for simple lookups.
    def __getitem__(self, keys):
        '''
        Look up records in the database by primary key. Also accepts :class:`slice` objects.

        Returns
        -------
        :class:`.record_type` or :class:`list` of :class:`.record_type`
        '''
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
            results.append(self.record_type.from_sql(record, database=self))
        return results if len(results) > 1 else results[0]

    def __iter__(self):
        '''
        Iterate sequentially over each entry in the database.
        '''
        for row in self.execute(self.stn("select * from {table_name};")):
            yield self.record_type.from_sql(row, database=self)

    def __contains__(self, key):
        try:
            self[key]
            return True
        except IndexError:
            return False

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

        `sqlite3.execute <https://docs.python.org/2/library/sqlite3.html#sqlite3.Connection.execute>`_
        '''
        return self.connection.execute(self.stn(query), *args, **kwargs)

    def executemany(self, query, param_iter, *args, **kwargs):
        '''
        A wrapper around :meth:`sqlite3.Connection.executemany`. Will format
        the query string to substitute in the main table name if the {table_name}
        token is present.

        `sqlite3.executemany <https://docs.python.org/2/library/sqlite3.html#sqlite3.Connection.executemany>`_
        '''
        return self.connection.executemany(self.stn(query), param_iter, *args, **kwargs)

    def executescript(self, script, *args, **kwargs):
        '''
        A wrapper around :meth:`sqlite3.Connection.executescript`.

        `sqlite3.executescript <https://docs.python.org/2/library/sqlite3.html#sqlite3.Connection.executescript>`_
        '''
        return self.connection.executescript(script, *args, **kwargs)

    def commit(self):
        '''
        A wrapper around :meth:`sqlite3.Connection.commit`. Writes pending changes to the database.

        `sqlite3.commit <https://docs.python.org/2/library/sqlite3.html#sqlite3.Connection.commit>`_
        '''
        self.connection.commit()

    def rollback(self):
        '''
        A wrapper around :meth:`sqlite3.Connection.rollback`. Reverses the last set of changes to the database.

        `sqlite3.rollback <https://docs.python.org/2/library/sqlite3.html#sqlite3.Connection.rollback>`_
        '''
        self.connection.rollback()

    def _find_boundaries(self, mass, tolerance):
        spread = mass * tolerance
        return (mass - spread, mass + spread)

    def ppm_match_tolerance_search(self, mass, tolerance, target_table=None, mass_shift=0):
        '''
        Rapidly search the database for entries with a recorded mass within ``tolerance`` ppm of ``mass``.
        '''
        target_table = target_table or self.record_type.table_name
        lower, upper = self._find_boundaries(mass + mass_shift, tolerance)
        results = self.execute("select * from {table_name}\
         where mass between {lower} and {upper};".format(
            lower=lower, upper=upper, table_name=target_table))
        for result in results:
            yield self.record_type.from_sql(result)
