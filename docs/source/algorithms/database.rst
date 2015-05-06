========================
Embedded Glycan Database
========================

We often want to work with large amounts of data, sometimes larger than is convenient to store in main memory. Python ships with builtin support for :emphasis:`sqlite`, available in the standard library. To ensure these features are available to all users, a simple *Object Relational Mapping* for the :mod:`sqlite3` module. The ORM system is built on top of a single class, :class:`GlycanRecord` which can be used to represent a collection of data that is partially exposed to `sqlite3`'s search engine through mapping functions.

For a more sophisticated ORM, please see http://www.sqlalchemy.org/

.. autoclass:: pygly2.algorithms.database.GlycanRecord
    :members:
    :inherited-members:
    :exclude-members: __str__, __getattribute__, __format__, __delattr__, __hash__, __new__, __reduce__, __reduce_ex__, __subclasshook__, __setattr__, __sizeof__, 


.. automodule:: pygly2.algorithms.database
    :members:
    :exclude-members: GlycanRecord, GlycanRecordBase

