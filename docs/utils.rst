Utilities and commonly reused generics
--------------------------------------

:mod:`glypy` reuses several common structures and functions throughout the
library. Some of these features are intended to be internal only. Those that
affect the object interfaces that users see are described here.


User-Facing Utilities
=====================

.. currentmodule:: glypy.utils.base

.. autofunction:: root

.. autofunction:: tree

Error Types
~~~~~~~~~~~

.. autoexception:: RootProtocolNotSupportedError
.. autoexception:: TreeProtocolNotSupportedError


Enum Type Implementation
========================
.. automodule:: glypy.utils.enum
    :exclude-members: _EnumMeta, _EnumValue

    .. autoclass:: EnumMeta

    .. autoclass:: EnumValue


Multimap Implementation
=======================

:class:`~.Monosaccharide` and :class:`~.Substituent` objects'
store position-specific infomration about links and modifiers
using a :class:`Mapping`-like object that allows a key to be
used to designate multiple values. These are usually
:class:`~.OrderedMultiMap` objects, which remember the order in
which keys were added to them. These types are implemented in
the :mod:`glypy.utils.multimap` module.

.. automodule:: glypy.utils.multimap

    .. autoclass:: MultiMap

    .. autoclass:: OrderedMultiMap

