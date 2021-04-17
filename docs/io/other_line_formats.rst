Other Line Formats
------------------

Glycomics is full of other notations for glycan sequences and compositions, with varying
nomenclatures and syntx for everything. Below is a list of modules that parse at least the
common cases for these formats.

.. automodule:: glypy.io.byonic
    :members: loads, dumps

.. automodule:: glypy.io.glyconnect
    :members: loads, dumps

.. automodule:: glypy.io.gws
    :members: loads, GWSError
    :exclude-members: parse_gws

.. automodule:: glypy.io.cfg
    :members: loads, CFGError
