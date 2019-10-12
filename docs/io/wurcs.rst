WURCS
=====


.. automodule:: glypy.io.wurcs
    :exclude-members: CarbonDescriptors, NodeTypeSpec, dumps, loads

    High Level Functions
    --------------------

     .. autofunction:: dumps

     .. autofunction:: loads


    File Parser
    -----------

    .. autoclass:: WURCSParser


    Implementation Details
    ----------------------

    .. autoclass:: NodeTypeSpec

    .. autoclass:: CarbonDescriptors

    .. autoexception:: WURCSError


    Low-Level Parser and Writer Implementations
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    .. autoclass:: glypy.io.wurcs.parser.WURCSParser

    .. autoclass:: glypy.io.wurcs.writer.WURCSWriter