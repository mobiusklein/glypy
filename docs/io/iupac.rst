IUPAC Three Letter Code
=======================


.. currentmodule:: glypy.io.iupac


.. automodule:: glypy.io.iupac

    .. autofunction:: dumps
    .. autofunction:: loads
    
    .. autoexception:: IUPACError

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

    .. autoclass:: SubstituentSerializer
    .. autoclass:: ModificationSerializer
    .. autoclass:: MonosaccharideSerializer
    .. autoclass:: GlycanSerializer

    .. autoclass:: DerivatizationAwareMonosaccharideSerializer

    .. autoclass:: SubstituentDeserializer
    .. autoclass:: ModificationDeserializer
    .. autoclass:: MonosaccharideDeserializer
    .. autoclass:: GlycanDeserializer

    .. autoclass:: DerivatizationAwareMonosaccharideDeserializer
    .. autoclass:: CondensedMonosaccharideDeserializer

