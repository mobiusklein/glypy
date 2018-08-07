import unittest

import glypy
from glypy.io import iupac
from glypy.enzyme import (
    MultiprocessingGlycome, make_n_glycan_pathway, EnzymeGraph)


class GlycomeTests(unittest.TestCase):

    def _make_glycome(self):
        glycosylases, glycosyltransferases, seeds = make_n_glycan_pathway()
        glycosyltransferases.pop("gntE")
        glycosyltransferases.pop('agal13galt')
        glycosyltransferases.pop('siat2_3')
        glycosyltransferases.pop('siat2_6')
        glycosyltransferases.pop("fuct3")
        glycome = MultiprocessingGlycome(
            glycosylases, glycosyltransferases, seeds)
        return glycome

    def test_path_between(self):
        with open("test_data/enzyme_graph.json", 'rt') as fh:
            graph = EnzymeGraph.load(fh)
        p = list(graph.seeds)[0]
        c = sorted(graph)[-1].child
        path = graph.path_between(p, c)
        assert len(path) == 13

    def test_synthesize_glycome(self):
        glycome = self._make_glycome()
        for i, gen in enumerate(glycome.run()):
            pass
        assert i == 20
        assert len(gen) == 1
        eg = EnzymeGraph(glycome.enzyme_graph)
        assert eg.node_count() == 498
        assert eg.edge_count() == 1427
        assert eg.parentless() == eg.seeds

        with open("test_data/enzyme_graph.json", 'rt') as fh:
            graph = EnzymeGraph.load(fh)
        assert eg == graph
        assert eg._dump() == graph._dump()
        seed = list(eg.seeds)[0]
        ref = list(graph.seeds)[0]
        assert seed == ref
        assert set(eg.children(seed)) == set(graph.children(seed))

    def test_galt(self):
        _, glycosyltransferases, _ = make_n_glycan_pathway()
        galt = glycosyltransferases['galt']
        target = iupac.loads(
            "a-D-Manp-(1-6)-[b-D-Glcp2NAc-(1-2)-a-D-Manp-(1-3)]b-D-Manp-(1-4)-b-D-Glcp2NAc-(1-4)-b-D-Glcp2NAc")
        for instance in galt(target):
            res = instance.mass() - glypy.monosaccharide_residues.Gal.mass()
            self.assertAlmostEqual(res, target.mass())


if __name__ == '__main__':
    unittest.main()
