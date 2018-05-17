import unittest

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


if __name__ == '__main__':
    unittest.main()
