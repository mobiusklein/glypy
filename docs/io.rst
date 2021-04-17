Nomenclature and Serialization Formats
======================================

There are many ways of representing monosaccharides, substituents and glycan
structures in text. :mod:`glypy.io` includes modules for reading and writing
several of these formats.

Three of these formats, :title-reference:`IUPAC`, :title-reference:`WURCS`, and
:title-reference:`LinearCode` put an entire structure in a single line, while
:title-reference:`GlycoCT` uses multi-line blocks to denote different parts of
a structure. The :mod:`~.glypy.io.glycoct`, :mod:`~.glypy.io.iupac`,
:mod:`~.glypy.io.linear_code`, and :mod:`~.glypy.io.wurcs` modules all provide
``loads`` and ``dumps`` functions, similar to other Python serialization interfaces
for converting objects to and from strings.

.. note::
    :title-reference:`GlycoCT` and :title-reference:`WURCS` support more complex
    representations, including both structures and compositions, than :title-reference:`IUPAC`,
    but all three can represent essentially any monosaccharide. :title-reference:`LinearCode`
    can only represent a limited number monosaccharides, not including generic cases with
    unknown ring stereochemistry.

.. warning::
    :title-reference:`IUPAC` and :title-reference:`LinearCode` have to do complex
    heuristic reasoning to decode modified versions of monosaccharides with special
    names (e.g. ``Neu5Ac``, ``Fuc``) that imply modifications or substituents, limiting
    their performance. :title-reference:`GlycoCT` and :title-reference:`WURCS` do not
    require nearly as much introspection, making them considerably faster.


.. toctree::
    :maxdepth: 2

    io/glycoct
    io/glycoct_xml
    io/wurcs
    io/iupac
    io/linear_code
    io/other_line_formats
    io/nomenclature/identity

