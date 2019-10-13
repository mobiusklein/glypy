IUPAC Three Letter Code
=======================

A reader and writer for the ubiquitious trivial IUPAC notation for glycan structures. This notation
will be familiar to biologists and biochemists, and the parser attempts to recognize common shorthand
and special cases that may be inconsistent with the grammmar used to describe the format.

.. currentmodule:: glypy.io.iupac


High Level Functions
--------------------

.. autofunction:: dumps
.. autofunction:: loads

.. autoexception:: IUPACError

There are multiple dialects of IUPAC, affecting how monosaccharides are written and how linkages are
denoted. The ``dialect`` parameter in these functions control which dialect is expected or generated.

The `simple` dialect: ``Neu5Gc(a2-3)Gal(b1-4)[Fuc(a1-3)]Glc2NAc``
The `extended` dialect: ``a-D-Neup5Gc-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1-3)]b-D-Glcp2NAc``


Example
-------

.. code-block:: python

    >>> from glypy.io import iupac
    >>> iupac.loads('a-L-Fucp-(1-6)-[a-D-Neup5Ac-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1-3)]b-D-Glcp2NAc-(1-6)-[a-D-Neup5Gc-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1-3)]b-D-Glcp2NAc-(1-2)]a-D-Manp-(1-6)-[b-D-Glcp2NAc-(1-4)][b-D-Galp2NAc-(1-4)-b-D-Glcp2NAc-(1-4)-[a-D-Neup5Gc-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1-3)]b-D-Glcp2NAc-(1-2)]a-D-Manp-(1-3)]b-D-Manp-(1-4)-b-D-Glcp2NAc-(1-4)]?-D-Glcp2NAc')
    >>>


.. plot::

    from glypy.io import iupac
    from glypy import plot

    glycan = iupac.loads('a-L-Fucp-(1-6)-[a-D-Neup5Ac-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1-3)]b-D-Glcp2NAc-(1-6)-[a-D-Neup5Gc-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1-3)]b-D-Glcp2NAc-(1-2)]a-D-Manp-(1-6)-[b-D-Glcp2NAc-(1-4)][b-D-Galp2NAc-(1-4)-b-D-Glcp2NAc-(1-4)-[a-D-Neup5Gc-(2-3)-b-D-Galp-(1-4)-[a-L-Fucp-(1-3)]b-D-Glcp2NAc-(1-2)]a-D-Manp-(1-3)]b-D-Manp-(1-4)-b-D-Glcp2NAc-(1-4)]?-D-Glcp2NAc')

    dt, ax = plot.plot(glycan, label=True)
    ax.figure.set_figwidth(8)
    ax.figure.set_figheight(4)
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo / 2, hi / 1.2)

Object-Oriented Interface
-------------------------

    The high-level API is implemented using a set of cooperating types for converting between
    :mod:`glypy`'s objects and IUPAC text formats. All of these objects provide these behaviors
    by calling them on the appropriate argument


Serialization
~~~~~~~~~~~~~

.. autoclass:: SubstituentSerializer
.. autoclass:: ModificationSerializer
.. autoclass:: MonosaccharideSerializer
.. autoclass:: LinkageSerializer
.. autoclass:: GlycanSerializer


Deserialization
~~~~~~~~~~~~~~~

.. autoclass:: SubstituentDeserializer
.. autoclass:: ModificationDeserializer
.. autoclass:: MonosaccharideDeserializer
.. autoclass:: GlycanDeserializer


Derivatization
~~~~~~~~~~~~~~

.. autoclass:: DerivatizationAwareMonosaccharideSerializer
.. autoclass:: DerivatizationAwareMonosaccharideDeserializer


Simplified Format
~~~~~~~~~~~~~~~~~


.. autoclass:: SimpleMonosaccharideSerializer
.. autoclass:: SimpleLinkageSerializer

.. autoclass:: SimpleMonosaccharideDeserializer
.. autoclass:: SimpleLinkageDeserializer

