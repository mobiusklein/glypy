import unittest

from glypy.io import linear_code
from glypy.tests.common import load


class LinearCodeTests(unittest.TestCase):

    def test_translate(self):
        broad = load("broad_n_glycan")
        dup = linear_code.loads(linear_code.dumps(broad))
        self.assertEqual(broad, dup)

        # linear code doesn't know about modifications or
        # ring shape
        sulfated = load("sulfated_glycan")
        sulfated.reducing_end = None
        sulfated.root.ring_start = 1
        sulfated.root.ring_end = 5
        dup2 = linear_code.loads(linear_code.dumps(sulfated))
        self.assertEqual(dup2, sulfated)

        sulfated = load("sulfated_glycan")
        dup = linear_code.loads(linear_code.dumps(sulfated))
        self.assertNotEqual(sulfated, dup)

    def test_example(self):
        seq = 'NNa3Ab3(NNa6)AN'
        structure = linear_code.loads(seq)
        self.assertEqual(len(structure), 4)


if __name__ == '__main__':
    unittest.main()
