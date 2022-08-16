import unittest
from .common import monosaccharides, link, glycoct, load
from glypy.composition.composition_transform import derivatize
from glypy.composition import Composition


class LinkTests(unittest.TestCase):

    def test_link_equality(self):
        parent = monosaccharides['Hex']
        child = monosaccharides['Hex']
        other = monosaccharides['Hex']
        link_1 = link.Link(parent, child, parent_position=3,
                           child_position=3, parent_loss='H', child_loss='OH')
        link_2 = link.Link(child, other, parent_position=6,
                           child_position=3, parent_loss='H', child_loss='OH')
        self.assertEqual(link_1, link_1)
        self.assertNotEqual(link_1, link_2)
        self.assertFalse(link_1 is None)

    def test_loss_composition(self):
        parent = monosaccharides['Hex']
        child = monosaccharides['Hex']

        link_1 = link.Link(parent, child, parent_position=3,
                           child_position=3, parent_loss='H', child_loss='OH')

        self.assertEqual(link_1.parent_loss, Composition(formula="H"))
        self.assertEqual(link_1.child_loss, Composition(formula="OH"))

    def test_break_and_reconnect(self):
        parent = monosaccharides['Hex']
        child = monosaccharides['Hex']

        link_1 = link.Link(parent, child, parent_position=3,
                           child_position=3, parent_loss='H', child_loss='OH')
        link_1.break_link(refund=True)
        self.assertTrue(len(parent.links[3]) == 0)
        self.assertTrue(len(child.links[3]) == 0)

        self.assertFalse(link_1.is_attached())
        self.assertFalse(link_1.is_attached(deep=True))

        link_1._reconnect(refund=True)
        self.assertTrue(len(parent.links[3]) == 1)
        self.assertTrue(len(child.links[3]) == 1)

        self.assertTrue(link_1.is_attached())
        self.assertTrue(link_1.is_attached(deep=True))

    def test_traversal(self):
        parent = monosaccharides['Hex']
        child = monosaccharides['Hex']

        link_1 = link.Link(parent, child, parent_position=3,
                           child_position=3, parent_loss='H', child_loss='OH')
        self.assertTrue(link_1.is_parent(parent))
        self.assertTrue(link_1.is_child(child))
        self.assertEqual(link_1.to(parent), child)
        self.assertEqual(link_1.to(child), parent)
        self.assertRaises(KeyError, lambda: link_1.to(1))


class AmbiguousLinkTests(unittest.TestCase):

    def test_link_site_no_overlap(self):

        ambig = '''
                RES
                1b:x-dglc-HEX-x:x
                2s:n-acetyl
                3b:b-dglc-HEX-1:5
                4s:n-acetyl
                5b:b-dman-HEX-1:5
                6b:a-dman-HEX-1:5
                7b:b-dglc-HEX-1:5
                8s:n-acetyl
                9b:b-dgal-HEX-1:5
                10b:a-dman-HEX-1:5
                11b:a-dman-HEX-1:5
                LIN
                1:1d(2+1)2n
                2:1o(4+1)3d
                3:3d(2+1)4n
                4:3o(4+1)5d
                5:5o(3|6+1)6d
                6:6o(2+1)7d
                7:7d(2+1)8n
                8:7o(4+1)9d
                9:5o(3|6+1)10d
                10:10o(3|6+1)11d
                '''
        exact = '''
                RES
                1b:x-dglc-HEX-1:5
                2s:n-acetyl
                3b:b-dglc-HEX-1:5
                4s:n-acetyl
                5b:b-dman-HEX-1:5
                6b:a-dman-HEX-1:5
                7b:b-dglc-HEX-1:5
                8s:n-acetyl
                9b:a-dman-HEX-1:5
                10b:a-dman-HEX-1:5
                11b:a-dman-HEX-1:5
                LIN
                1:1d(2+1)2n
                2:1o(4+1)3d
                3:3d(2+1)4n
                4:3o(4+1)5d
                5:5o(3+1)6d
                6:6o(2+1)7d
                7:7d(2+1)8n
                8:5o(6+1)9d
                9:9o(3+1)10d
                10:9o(6+1)11d
                '''
        ambig = glycoct.loads(ambig)
        exact = glycoct.loads(exact)
        ambig.set_reducing_end(True)
        exact.set_reducing_end(True)
        self.assertEqual(ambig.total_composition(), exact.total_composition())
        self.assertEqual(ambig, ambig.clone())
        self.assertEqual(derivatize(ambig.clone(), 'methyl').total_composition(),
                         derivatize(exact.clone(), 'methyl').total_composition())
