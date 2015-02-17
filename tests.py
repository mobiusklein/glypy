import unittest
import json
import logging
import cPickle
#logging.basicConfig(level="DEBUG")

import pygly2

from pygly2.structure import monosaccharide, constants, substituent, glycan, link, named_structures, structure_composition
from pygly2.io import glycoct
from pygly2.utils import enum, StringIO, identity, nullop, multimap
from pygly2.composition import Composition, composition_transform, composition


monosaccharides = named_structures.monosaccharides


class GlycoCTParserTests(unittest.TestCase):
    _file_path = "./test_data/glycoct.txt"

    def test_parse_file(self):
        for g in glycoct.read(self._file_path):
            self.assertTrue(isinstance(g, glycan.Glycan))


# monosaccharide_masses = json.load(open("./test_data/monosaccharide_masses.json"))
monosaccharide_structures = json.load(open("./pygly2/structure/data/monosaccharides.json"))

wiki_masses = {
    "Iduronic Acid": 194.04,
    "Bacillosamine": 162.10,
    "Allose": 180.06,
}

#: Not common in any way other than
#: reused in many tests
common_glycan = '''
RES
1b:b-dglc-HEX-1:5
2b:b-dgal-HEX-1:5
3b:b-dglc-HEX-1:5
4s:n-acetyl
5b:a-lgal-HEX-1:5|6:d
6b:b-dgal-HEX-1:5
7b:b-dglc-HEX-1:5
8s:n-acetyl
9b:a-lgal-HEX-1:5|6:d
10b:b-dgal-HEX-1:5
LIN
1:1o(4+1)2d
2:2o(3+1)3d
3:3d(2+1)4n
4:3o(3+1)5d
5:3o(4+1)6d
6:6o(3+1)7d
7:7d(2+1)8n
8:7o(3+1)9d
9:7o(4+1)10d'''

branchy_glycan = '''
RES
1b:x-dglc-HEX-x:x
2s:n-acetyl
3b:b-dman-HEX-1:5
4b:a-dman-HEX-1:5
5b:b-dglc-HEX-1:5
6s:n-acetyl
7b:b-dgal-HEX-1:5
8b:b-dglc-HEX-1:5
9s:n-acetyl
10b:b-dgal-HEX-1:5
11b:a-dman-HEX-1:5
12b:b-dglc-HEX-1:5
13s:n-acetyl
14b:b-dgal-HEX-1:5
LIN
1:1d(2+1)2n
2:1o(4+1)3d
3:3o(3+1)11d
4:3o(6+1)4d
5:4o(2+1)8d
6:4o(6+1)5d
7:5d(2+1)6n
8:5o(4+1)7d
9:8d(2+1)9n
10:8o(4+1)10d
11:11o(2+1)12d
12:12d(2+1)13n
13:12o(4+1)14d'''


class NamedStructureTests(unittest.TestCase):
    def test_accessors(self):
        self.assertEqual(named_structures.monosaccharides.Fucose,
                         named_structures.monosaccharides["Fucose"])
        self.assertNotEqual(named_structures.monosaccharides.Man,
                            named_structures.monosaccharides.Hex)
        self.assertEqual(named_structures.glycans["N-Linked Core"],
                         named_structures.glycans["N-Linked Core"])


class StructureCompositionTests(unittest.TestCase):
    def test_missing_composition(self):
        import warnings
        warnings.filterwarnings('ignore')
        self.assertEqual({}, structure_composition.monosaccharide_composition['hexadodecose'])
        warnings.filterwarnings('always')


class MonosaccharideTests(unittest.TestCase):
    _file_path = "./test_data/glycoct.txt"
    glycan = iter(glycoct.read(_file_path)).next()

    def test_from_glycoct(self):
        s = self.glycan.root.to_glycoct()
        b = StringIO(s)
        g = iter(glycoct.read(b)).next()
        self.assertEqual(g.root.to_glycoct(), s)

    def test_named_structure_masses(self):
        for name, mass in wiki_masses.items():
            structure = named_structures.monosaccharides[name]
            self.assertAlmostEqual(mass, structure.mass(), 2)

    def test_named_structure_glycoct(self):
        for name, glycoct_str in monosaccharide_structures.items():
            structure = named_structures.monosaccharides[name]
            test, ref = (structure.to_glycoct(), glycoct_str)
            i = 0
            for j, k in zip(test, ref):
                if j != k:
                    test_loc = test.replace('\n', ' ')[i-10:i+10]
                    ref_loc = ref.replace('\n', ' ')[i-10:i+10]
                    raise AssertionError("{j} != {k} at {i} in {name}\n{test_loc}\n{ref_loc}".format(**locals()))
                i += 1

    def test_ring_limit_modification(self):
        structure = named_structures.monosaccharides['Hex']
        self.assertRaises(IndexError, lambda: structure.add_modification('d', 8))

    def test_occupancy_limit_modification(self):
        structure = named_structures.monosaccharides['Hex']
        structure.add_modification('d', 4)
        self.assertRaises(ValueError, lambda: structure.add_modification('d', 4))

    def test_ring_limit_substituent(self):
        structure = named_structures.monosaccharides['Hex']
        self.assertRaises(IndexError, lambda: structure.add_substituent('methyl', 8))

    def test_occupancy_limit_substituent(self):
        structure = named_structures.monosaccharides['Hex']
        structure.add_substituent('methyl', 4)
        self.assertRaises(ValueError, lambda: structure.add_substituent('methyl', 4))

    def test_add_remove_modifcations(self):
        for name, mass in wiki_masses.items():
            structure = named_structures.monosaccharides[name]
            ref = structure.clone()
            self.assertAlmostEqual(ref.mass(), structure.mass())
            open_sites, unknowns = structure.open_attachment_sites()
            n_sites = len(open_sites)
            comp_delta = n_sites * structure_composition.modification_compositions[constants.Modification.aldi]
            for site in open_sites:
                structure.add_modification(constants.Modification.aldi, site)
            self.assertEqual(structure.total_composition(), ref.total_composition() + comp_delta)
            self.assertEqual([], structure.open_attachment_sites()[0])
            for site in open_sites:
                structure.drop_modification(site, constants.Modification.aldi)
            self.assertEqual(structure.total_composition(), ref.total_composition())

    def test_add_remove_substituents(self):
        for name, mass in wiki_masses.items():
            structure = named_structures.monosaccharides[name]
            ref = structure.clone()
            self.assertAlmostEqual(ref.mass(), structure.mass())
            open_sites, unknowns = structure.open_attachment_sites()
            n_sites = len(open_sites)
            mass_delta = substituent.Substituent('methyl').mass() * n_sites - Composition("H2").mass * n_sites
            ping = True
            for site in open_sites:
                if ping:
                    structure.add_substituent(substituent.Substituent('methyl'), position=site)
                    ping = False
                else:
                    structure.add_substituent('methyl', position=site)
                    ping = True
            self.assertAlmostEqual(structure.mass(), ref.mass() + mass_delta)
            for site in open_sites:
                structure.drop_substituent(site, substituent.Substituent('methyl'))
            self.assertAlmostEqual(structure.mass(), ref.mass())

    def test_validate_drop_substituents(self):
        structure = named_structures.monosaccharides['Hex']
        self.assertRaises(IndexError, lambda: structure.drop_substituent(99))
        self.assertRaises(IndexError, lambda: structure.drop_substituent(4))
        self.assertRaises(IndexError, lambda: structure.add_substituent(
            "n-acetyl", 3, max_occupancy=4).add_substituent(
            "methyl", 3, max_occupancy=4).drop_substituent(3, "n-glycolyl"))

    def test_add_remove_monosaccharides(self):
        for name, mass in wiki_masses.items():
            structure = named_structures.monosaccharides[name]
            ref = structure.clone()
            self.assertAlmostEqual(ref.mass(), glycan.Glycan(structure).mass())
            open_sites, unknowns = structure.open_attachment_sites()
            n_sites = len(open_sites)
            mass_delta = named_structures.monosaccharides["Hex"].mass() * n_sites - Composition("H2O").mass * n_sites
            for site in open_sites:
                structure.add_monosaccharide(named_structures.monosaccharides["Hex"], position=site, child_position=3)
            self.assertAlmostEqual(glycan.Glycan(structure).mass(), ref.mass() + mass_delta)
            for site in open_sites:
                structure.drop_monosaccharide(site)
            self.assertAlmostEqual(glycan.Glycan(structure).mass(), ref.mass())

    def test_validate_drop_monosacharide(self):
        structure = named_structures.monosaccharides['Hex']
        self.assertRaises(IndexError, lambda: structure.drop_monosaccharide(99))
        self.assertRaises(ValueError, lambda: structure.drop_monosaccharide(4))
        self.assertRaises(ValueError, lambda: structure.add_monosaccharide(
            named_structures.monosaccharides["Hex"], 3, max_occupancy=4).add_monosaccharide(
            named_structures.monosaccharides["Hex"], 3, max_occupancy=4).drop_monosaccharide(3))

    def test_validate_enums(self):
        structure = named_structures.monosaccharides['Hex']
        def t():
            structure.anomer = 5
        self.assertRaises(Exception, t)
        def t():
            structure.stem = "gibberish"
        self.assertRaises(Exception, t)
        def t():
            structure.superclass = "monose"
        self.assertRaises(Exception, t)
        def t():
            structure.configuration = 'not-real'
        self.assertRaises(Exception, t)

    def test_validate_reducing_end(self):
        structure = named_structures.monosaccharides['Hex']
        structure.reducing_end = 1
        structure.reducing_end = 3

        def t():
            structure.reducing_end = -1
        self.assertRaises(ValueError, t)

        def t():
            structure.reducing_end = 99
        self.assertRaises(IndexError, t)


class GlycanTests(unittest.TestCase):
    _file_path = "./test_data/glycoct.txt"

    def test_from_glycoct(self):
        for structure in glycoct.read(self._file_path):
            self.assertEqual(structure, glycoct.loads(structure.to_glycoct()).next())

    def test_fragments_preserve(self):
        for structure in glycoct.read(self._file_path):
            dup = structure.clone()
            self.assertEqual(structure, dup)
            list(dup.fragments('ZCBY', 3))

            self.assertEqual(structure, dup)

    def test_branch_counts(self):
        structure = glycoct.loads(branchy_glycan).next()
        self.assertEqual(structure.count_branches(), 3)

    def test_fragments_mass(self):
        structure = glycoct.loads(common_glycan).next()
        frags = list(structure.fragments('ZCBY', 1))
        import json
        frags_ref = json.load(open('test_data/fragments-example.json'))
        container = multimap.MultiMap()
        for frag in frags_ref:
            container[frag['kind']] = frag['mass']

        def almost_equal(a, b):
            e = 0.0001
            return (a - e) <= b <= (a + e)

        for frag in frags:
            structure.name_fragment(frag)
            kind = frag[0]
            mass = frag[-1]
            candidates = container[kind] + container[kind[::-1]]
            if len(candidates) == 0:
                raise AssertionError("No candidates found for {}".format(frag))
            res = (any(almost_equal(mass, x) for x in candidates))
            if not res:
                raise AssertionError("{} found no matches in {}".format(frag, candidates))

    def test_disjoint_subtrees(self):
        structure = glycoct.loads(common_glycan).next()
        f = open('test_data/test_disjoint_subtrees.pkl', 'rb')
        for tree in structure.fragments("BY", 3, structures=True):
            ref = cPickle.load(f)
            self.assertEqual(tree.link_ids, ref.link_ids)
            self.assertEqual(tree.parent_tree, ref.parent_tree)

    def test_reducing_end(self):
        structure = glycoct.loads(common_glycan).next()
        self.assertEqual(structure.reducing_end, None)
        structure.reducing_end = 1
        self.assertEqual(structure.reducing_end, 1)

    def test_clone(self):
        structure = glycoct.loads(common_glycan).next()
        ref = structure.clone()
        structure.reducing_end = 1
        self.assertTrue(structure != ref)

    def test_indexing(self):
        structure = glycoct.loads(common_glycan).next()
        ref = structure.clone()
        for i, node in enumerate(structure.index):
            self.assertEqual(node.id, ref[i].id)
        structure.deindex()
        for i, node in enumerate(structure.index):
            self.assertNotEqual(node.id, ref[i].id)
        structure.index = None
        self.assertEqual(structure[0], structure.root)


    def test_traversal(self):
        structure = glycoct.loads(common_glycan).next()
        structure[-1].add_monosaccharide(named_structures.monosaccharides['NeuGc'])
        structure.reindex(method='dfs')
        ref = structure.clone()
        self.assertEqual(structure[-1], ref[-1])
        structure.reindex(method='bfs')
        self.assertNotEqual(structure[-1], ref[-1])

    def test_traversal_by_name(self):
        structure = glycoct.loads(common_glycan).next()
        structure[-1].add_monosaccharide(named_structures.monosaccharides['NeuGc'])
        structure.reindex(method='dfs')
        ref = structure.clone()
        self.assertEqual(structure, ref)
        structure.reindex(method='depth_first_traversal')
        self.assertEqual(structure, ref)
        self.assertRaises(AttributeError, lambda: structure.reindex(method='not_real_traversal'))

    def test_leaves(self):
        structure = glycoct.loads(common_glycan).next()
        leaves = list(structure.leaves())
        for node in leaves:
            self.assertTrue(len(list(node.children())) == 0)
        leaves = list(structure.leaves(bidirectional=True))
        for node in leaves:
            self.assertTrue(len(list(node.children())) == 0 or node == structure.root)

    def test_custom_traversal_method(self):
        def rev_sort_dfs(self, visited=None, from_node=None, *args, **kwargs):
            node_stack = list([self.root])
            visited = set()
            while len(node_stack) > 0:
                node = node_stack.pop()
                if node.id in visited:
                    continue
                visited.add(node.id)
                yield (node)
                node_stack.extend(reversed(list(terminal for pos, link in node.links.items()
                                  for terminal in link if terminal.id not in visited and
                                  len(link.child.substituent_links) < 1)))
        structure = glycoct.loads(common_glycan).next()
        structure[-3].add_monosaccharide(named_structures.monosaccharides['Hex'], 4).add_substituent('methyl', 5)
        structure.reindex()
        ref = structure.clone()
        self.assertEqual(structure[-1], ref[-1])
        structure.reindex(method=rev_sort_dfs)
        self.assertNotEqual(structure[-1], ref[-1])

    def test_topological_equality(self):
        base = glycoct.loads(branchy_glycan).next()
        a = base.clone()
        b = base.clone()
        c = base.clone()
        self.assertEqual(base, b)
        self.assertEqual(a, b)
        list(a.leaves())[0].add_monosaccharide(monosaccharides["NeuGc"])
        list(b.leaves())[1].add_monosaccharide(monosaccharides["NeuGc"])
        list(c.leaves())[2].add_monosaccharide(monosaccharides["NeuGc"])
        d = base.clone()
        d_children = list(d.leaves())
        d_children[0].add_monosaccharide(monosaccharides["NeuGc"])
        d_children[1].add_monosaccharide(monosaccharides["NeuGc"])
        self.assertEqual(a, b)
        self.assertNotEqual(a, c)
        self.assertNotEqual(a, base)
        self.assertNotEqual(a, d)
        self.assertNotEqual(b, d)


class SubstituentTests(unittest.TestCase):
    def test_substituent_substituent(self):
        parent = substituent.Substituent("n-acetyl")
        pmass = parent.mass()
        child = substituent.Substituent("methyl")
        cmass = child.mass()

        parent.add_substituent(child, 2)
        self.assertAlmostEqual(parent.mass(), (pmass + cmass - Composition(H=2).mass))
        self.assertEqual(parent.total_composition(), parent.composition + child.composition)
        self.assertEqual(parent.children().next()[1], child)
        self.assertEqual(child.parents().next()[1], parent)
        parent.drop_substiuent(2)
        self.assertAlmostEqual(parent.mass(), pmass)

    def test_equality(self):
        sub_1 = substituent.Substituent("n-acetyl")
        sub_2 = substituent.Substituent("n-acetyl")
        self.assertEqual(sub_1, sub_2)
        parent = named_structures.monosaccharides['Hex']
        parent.add_substituent(sub_1, 3)
        self.assertNotEqual(sub_1, sub_2)

    def test_clone(self):
        parent = substituent.Substituent("n-acetyl")
        child = substituent.Substituent("methyl")

        parent.add_substituent(child, 2)

        dup = parent.clone()
        self.assertEqual(parent, dup)
        self.assertNotEqual(child, child.clone())


class MultiMapTests(unittest.TestCase):
    def test_iterators(self):
        from collections import Counter
        mm = multimap.MultiMap(a=1, b=3)
        mm['a'] = 3
        self.assertTrue(set(mm.keys()) == {'a', 'b'})
        self.assertTrue(set(mm.items()) == {('a', 1), ('a', 3), ('b', 3)})
        self.assertTrue(Counter(mm.values()) == Counter({1: 1, 3: 2}))

    def test_ordered_multimap(self):
        mm = multimap.MultiMap(a=1, b=3)
        mm['a'] = 3
        omm = multimap.OrderedMultiMap(a=1, b=3)
        omm['a'] = 3
        self.assertEqual(mm, omm)
        omm['c'] = 1
        self.assertNotEqual(mm, omm)
        self.assertNotEqual(omm, mm)
        self.assertFalse('c' in mm)


class CompositionTests(unittest.TestCase):

    def test_derivativize_bare(self):
        permethylated_reduced_mass = 1716.9033
        glycan = glycoct.loads(common_glycan).next()
        glycan.reducing_end = 1
        composition_transform.derivatize(glycan, 'methyl')
        self.assertAlmostEqual(glycan.mass(), permethylated_reduced_mass, 3)

    def test_composition_equality(self):
        self.assertEqual(Composition("H2O"), Composition("H2O"))
        self.assertEqual(Composition("H2O"), Composition("(H2O)"))

    def test_composition_substraction(self):
        self.assertEqual(Composition("NH2O") - Composition("N"), Composition("H2O"))

    def test_isotope_parsing(self):
        self.assertFalse(Composition("O[18]") == Composition("O"))
        self.assertAlmostEqual(Composition("O[18]").mass, 17.999, 3)
        self.assertRaises(composition.ChemicalCompositionError, lambda: Composition("O[18.5]"))

    def test_inits(self):
        self.assertEqual(Composition(O=1, H=2), Composition(formula='H2O'))
        self.assertEqual(Composition(O=1, H=2), Composition("OH2"))
        self.assertRaises(composition.ChemicalCompositionError, lambda: Composition("O2.2"))

    def test_operators(self):
        case = Composition("H2O")
        self.assertEqual(case + {'O': 1}, Composition("H2O2"))
        self.assertEqual({'O': 1} + case, Composition("H2O2"))

        self.assertEqual(case - {'O': 1}, Composition("H2"))
        self.assertEqual({'O': 1, "H": 4} - case, Composition("H2"))

        self.assertEqual(case * 3, Composition("H6O3"))
        self.assertEqual(case * -1, -case)
        self.assertRaises(composition.ChemicalCompositionError, lambda: case * 5.2)

    def test_massing(self):
        case = Composition("H2O")
        mono_mass = case.calc_mass()
        avg_mass = case.calc_mass(average=True)
        protonated = case + {"H+": 1}
        prot_mass = protonated.calc_mass()
        self.assertNotEqual(mono_mass, avg_mass)
        self.assertAlmostEqual(mono_mass, 18.0105, 3)
        self.assertAlmostEqual(avg_mass, 18.01528, 3)
        self.assertNotEqual(mono_mass, prot_mass)
        self.assertAlmostEqual(prot_mass, 19.01784, 3)
        self.assertAlmostEqual(case.calc_mass(charge=1), 19.01784, 3)
        self.assertRaises(composition.ChemicalCompositionError, lambda: protonated.calc_mass(charge=1))


class ConstantTests(unittest.TestCase):
    def test_translate(self):
        self.assertTrue(constants.Modification.d == constants.Modification['d'])
        self.assertTrue(constants.Modification.d == constants.Modification['D'])
        self.assertTrue(constants.Modification.d == constants.Modification[constants.Modification.d.value])

    def test_compare(self):
        self.assertTrue(constants.Modification.d == 'd')
        self.assertTrue(constants.Modification.d == constants.Modification.d.value)
        self.assertNotEqual(constants.Modification.d, constants.Stem[constants.Modification.d.value])
        self.assertRaises(KeyError, lambda: constants.SuperClass[1])


class LinkTests(unittest.TestCase):
    def test_link_equality(self):
        parent = named_structures.monosaccharides['Hex']
        child = named_structures.monosaccharides['Hex']
        other = named_structures.monosaccharides['Hex']
        link_1 = link.Link(parent, child, parent_position=3, child_position=3, parent_loss='H', child_loss='OH')
        link_2 = link.Link(child, other, parent_position=6, child_position=3, parent_loss='H', child_loss='OH')
        self.assertEqual(link_1, link_1)
        self.assertNotEqual(link_1, link_2)
        self.assertFalse(link_1 is None)

    def test_loss_composition(self):
        parent = named_structures.monosaccharides['Hex']
        child = named_structures.monosaccharides['Hex']

        link_1 = link.Link(parent, child, parent_position=3, child_position=3, parent_loss='H', child_loss='OH')

        self.assertEqual(link_1.parent_loss, Composition(formula="H"))
        self.assertEqual(link_1.child_loss, Composition(formula="OH"))

    def test_break_and_reconnect(self):
        parent = named_structures.monosaccharides['Hex']
        child = named_structures.monosaccharides['Hex']

        link_1 = link.Link(parent, child, parent_position=3, child_position=3, parent_loss='H', child_loss='OH')
        link_1.break_link(refund=True, reorient_fn=identity)
        self.assertTrue(len(parent.links[3]) == 0)
        self.assertTrue(len(child.links[3]) == 0)

        link_1.reconnect(refund=True, reorient_fn=identity)
        self.assertTrue(len(parent.links[3]) == 1)
        self.assertTrue(len(child.links[3]) == 1)

    def test_traversal(self):
        parent = named_structures.monosaccharides['Hex']
        child = named_structures.monosaccharides['Hex']

        link_1 = link.Link(parent, child, parent_position=3, child_position=3, parent_loss='H', child_loss='OH')
        self.assertTrue(link_1.is_parent(parent))
        self.assertTrue(link_1.is_child(child))
        self.assertEqual(link_1.to(parent), child)
        self.assertEqual(link_1.to(child), parent)
        self.assertRaises(KeyError, lambda: link_1.to(1))

if __name__ == '__main__':
    unittest.main()
