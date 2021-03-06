GlycoCT
=======

.. currentmodule:: glypy.io.glycoct

.. automodule:: glypy.io.glycoct
    :no-members:

    High Level Functions
    --------------------

    .. autofunction:: dump

    .. autofunction:: load

    .. autofunction:: dumps

    .. autofunction:: loads

    .. autoexception:: GlycoCTError

    Examples
    --------

    .. code-block:: python

        >>> from glypy.io import glycoct
        >>> glycoct.loads("""RES
        1b:x-dglc-HEX-1:5
        2s:n-acetyl
        3b:b-dglc-HEX-1:5
        4s:n-acetyl
        5b:b-dman-HEX-1:5
        6b:a-dman-HEX-1:5
        7b:b-dglc-HEX-1:5
        8s:n-acetyl
        9b:a-lgal-HEX-1:5|6:d
        10b:b-dgal-HEX-1:5
        11b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
        12s:n-glycolyl
        13b:b-dglc-HEX-1:5
        14s:n-acetyl
        15b:b-dgal-HEX-1:5
        16s:n-acetyl
        17b:b-dglc-HEX-1:5
        18s:n-acetyl
        19b:a-dman-HEX-1:5
        20b:b-dglc-HEX-1:5
        21s:n-acetyl
        22b:a-lgal-HEX-1:5|6:d
        23b:b-dgal-HEX-1:5
        24b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
        25s:n-glycolyl
        26b:b-dglc-HEX-1:5
        27s:n-acetyl
        28b:a-lgal-HEX-1:5|6:d
        29b:b-dgal-HEX-1:5
        30b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
        31s:n-acetyl
        32b:a-lgal-HEX-1:5|6:d
        LIN
        1:1d(2+1)2n
        2:1o(4+1)3d
        3:3d(2+1)4n
        4:3o(4+1)5d
        5:5o(3+1)6d
        6:6o(2+1)7d
        7:7d(2+1)8n
        8:7o(3+1)9d
        9:7o(4+1)10d
        10:10o(3+2)11d
        11:11d(5+1)12n
        12:6o(4+1)13d
        13:13d(2+1)14n
        14:13o(4+1)15d
        15:15d(2+1)16n
        16:5o(4+1)17d
        17:17d(2+1)18n
        18:5o(6+1)19d
        19:19o(2+1)20d
        20:20d(2+1)21n
        21:20o(3+1)22d
        22:20o(4+1)23d
        23:23o(3+2)24d
        24:24d(5+1)25n
        25:19o(6+1)26d
        26:26d(2+1)27n
        27:26o(3+1)28d
        28:26o(4+1)29d
        29:29o(3+2)30d
        30:30d(5+1)31n
        31:1o(6+1)32d
        """)
        >>>

    .. plot::

        from glypy.io import glycoct
        from glypy import plot

        print("Plotting...")

        glycan = glycoct.loads("""RES
            1b:x-dglc-HEX-1:5
            2s:n-acetyl
            3b:b-dglc-HEX-1:5
            4s:n-acetyl
            5b:b-dman-HEX-1:5
            6b:a-dman-HEX-1:5
            7b:b-dglc-HEX-1:5
            8s:n-acetyl
            9b:a-lgal-HEX-1:5|6:d
            10b:b-dgal-HEX-1:5
            11b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
            12s:n-glycolyl
            13b:b-dglc-HEX-1:5
            14s:n-acetyl
            15b:b-dgal-HEX-1:5
            16s:n-acetyl
            17b:b-dglc-HEX-1:5
            18s:n-acetyl
            19b:a-dman-HEX-1:5
            20b:b-dglc-HEX-1:5
            21s:n-acetyl
            22b:a-lgal-HEX-1:5|6:d
            23b:b-dgal-HEX-1:5
            24b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
            25s:n-glycolyl
            26b:b-dglc-HEX-1:5
            27s:n-acetyl
            28b:a-lgal-HEX-1:5|6:d
            29b:b-dgal-HEX-1:5
            30b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
            31s:n-acetyl
            32b:a-lgal-HEX-1:5|6:d
            LIN
            1:1d(2+1)2n
            2:1o(4+1)3d
            3:3d(2+1)4n
            4:3o(4+1)5d
            5:5o(3+1)6d
            6:6o(2+1)7d
            7:7d(2+1)8n
            8:7o(3+1)9d
            9:7o(4+1)10d
            10:10o(3+2)11d
            11:11d(5+1)12n
            12:6o(4+1)13d
            13:13d(2+1)14n
            14:13o(4+1)15d
            15:15d(2+1)16n
            16:5o(4+1)17d
            17:17d(2+1)18n
            18:5o(6+1)19d
            19:19o(2+1)20d
            20:20d(2+1)21n
            21:20o(3+1)22d
            22:20o(4+1)23d
            23:23o(3+2)24d
            24:24d(5+1)25n
            25:19o(6+1)26d
            26:26d(2+1)27n
            27:26o(3+1)28d
            28:26o(4+1)29d
            29:29o(3+2)30d
            30:30d(5+1)31n
            31:1o(6+1)32d""")
        dt, ax = plot.plot(glycan, label=True)
        ax.figure.set_figwidth(8)
        ax.figure.set_figheight(4)
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo / 2, hi / 1.2)

    Object-Oriented Interface
    -------------------------

    .. autoclass:: GlycoCTReader

    .. autoclass:: GlycoCTWriter


    Implementation Details
    ----------------------

    .. autoclass:: RepeatedGlycoCTSubgraph

    .. autoclass:: UndeterminedGlycoCTSubgraph
