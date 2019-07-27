import zlib

try:
    from collections.abc import MutableSet
except ImportError:
    from collections import MutableSet

from glypy.io.glycoct import (loads, dumps)


class DistinctGlycanSet(MutableSet):
    """Store a distinct set of unique :class:`~.Glycan` objects
    space-efficiently.

    Implements the :class:`MutableSet` interface.

    This type stores structures by serializing the :class:`~.Glycan`
    into :title-reference:`GlycoCT{condensed}` text, and then compresses
    the text using :func:`zlib.compress`.

    Attributes
    ----------
    raw_data_buffer: :class:`set`
        The container for managing the unique glycans
    """

    def __init__(self, structures=None):
        if structures is None:
            structures = []
        self.raw_data_buffer = set()

        if isinstance(structures, DistinctGlycanSet):
            self.update(structures)
        else:
            for structure in structures:
                self.add(structure)

    def add(self, structure):
        """Add `structure` to the set

        Parameters
        ----------
        structure: :class:`~.Glycan`
            The structure to add to the set

        Returns
        -------
        :class:`bytes`
            The compressed encoding of `structure`
        """
        key = self._transform_text(
            self._structure_to_text(structure))
        if key in self.raw_data_buffer:
            return key
        self.raw_data_buffer.add(key)
        return key

    def discard(self, structure):
        """Remove `structure` from the set

        Parameters
        ----------
        structure: :class:`~.Glycan`
            The structure to remove
        """
        key = self._transform_text(
            self._structure_to_text(structure))
        self.raw_data_buffer.discard(key)

    def _structure_to_text(self, structure):
        return dumps(structure).encode('utf-8')

    def _text_to_structure(self, text):
        return loads(text.decode('utf-8'))

    def _transform_text(self, text):
        return zlib.compress(text)

    def _untransform_text(self, compressed):
        return zlib.decompress(compressed)

    def encode(self, structure):
        """Encode `structure` into compressed bytes

        Parameters
        ----------
        structure: :class:`~.Glycan`
            The structure to encode

        Returns
        -------
        :class:`bytes`
            The compressed bytes encoding `structure`
        """
        return self._transform_text(self._structure_to_text(structure))

    def add_encoded(self, encoded):
        """Add a pre-encoded structure to the set

        Parameters
        ----------
        encoded: :class:`bytes`
            An encoded structure
        """
        self.raw_data_buffer.add(encoded)

    def discard_encoded(self, encoded):
        """Remove a pre-encoded structure from the set

        Parameters
        ----------
        encoded: :class:`bytes`
            An encoded structure
        """
        self.raw_data_buffer.discard(encoded)

    def has_encoded(self, encoded):
        return encoded in self.raw_data_buffer

    def pop(self):
        text = self._untransform_text(self.raw_data_buffer.pop())
        return self._text_to_structure(text)

    def __len__(self):
        return len(self.raw_data_buffer)

    def __iter__(self):
        for compressed in self.raw_data_buffer:
            yield self._text_to_structure(
                self._untransform_text(compressed))

    def __contains__(self, structure):
        text = self._structure_to_text(structure)
        compressed = self._transform_text(text)
        return compressed in self.raw_data_buffer

    @classmethod
    def from_buffer_slice(cls, buffer_slice):
        inst = cls()
        inst.raw_data_buffer.update(buffer_slice)
        return inst

    def update(self, other):
        if isinstance(other, DistinctGlycanSet):
            self.raw_data_buffer.update(other.raw_data_buffer)
        else:
            for x in other:
                self.add(x)

    def partition(self, nchunks=2):
        buffer_sequence = list(self.raw_data_buffer)
        n = len(buffer_sequence)
        chunk_size = max(int(n / nchunks), 1)
        i = 0
        chunks = []
        while i < n:
            chunk = buffer_sequence[i:i + chunk_size]
            i += chunk_size
            chunks.append(self.from_buffer_slice(chunk))
        return chunks

    def remove_all(self, other):
        if not isinstance(other, DistinctGlycanSet):
            other = DistinctGlycanSet(other)
        self.raw_data_buffer -= other.raw_data_buffer

    def __sub__(self, other):
        return self.from_buffer_slice(
            self.raw_data_buffer - other.raw_data_buffer)

    def __isub__(self, other):
        self.raw_data_buffer -= other.raw_data_buffer
        return self

    def __and__(self, other):
        return self.from_buffer_slice(
            self.raw_data_buffer & other.raw_data_buffer)

    def __iand__(self, other):
        self.raw_data_buffer &= other.raw_data_buffer
        return self

    def __or__(self, other):
        return self.from_buffer_slice(
            self.raw_data_buffer | other.raw_data_buffer)

    def __ior__(self, other):
        self.raw_data_buffer.update(other.raw_data_buffer)
        return self
