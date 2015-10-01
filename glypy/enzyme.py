from glypy import Glycan
from glypy import monosaccharides
from glypy.structure.constants import Anomer, Configuration
from glypy.algorithms.similarity import commutative_similarity


class Glycosyltransferase(object):
    def __init__(self, parent_position, child_position, anomer, parent, child, comparator=commutative_similarity):
        self.parent_position = parent_position
        self.child_position = child_position
        self.anomer = anomer
        self.parent = parent
        self.child = child
        self.comparator = comparator

    def _traverse(self, structure):
        for node in structure:
            if (self.parent is None or self.comparator(node, self.parent)) and not node.is_occupied(
                    self.parent_position):
                yield node

    def traverse(self, structure):
        return list(self._traverse(structure))

    def apply(self, monosaccharide, parent_position=None, child_position=None):
        new_monosaccharide = self.child.clone()
        new_monosaccharide.anomer = self.anomer
        monosaccharide.add_monosaccharide(
            new_monosaccharide, position=self.parent_position,
            child_position=self.child_position)

    def __call__(self, structure, parent_position=None, child_position=None):
        for node in self.traverse(structure):
            new_structure = structure.clone()
            node = new_structure.get(node.id)
            self.apply(node, parent_position, child_position)
            yield new_structure


class Glycosidase(object):
    def __init__(self, parent_position, child_position, anomer, parent, child, comparator=commutative_similarity):
        self.parent_position = parent_position
        self.child_position = child_position
        self.anomer = anomer
        self.parent = parent
        self.child = child
        self.comparator = comparator

    def _traverse(self, structure):
        for p, link in structure.iterlinks():
            if (self.parent is None or self.comparator(link.parent, self.parent)) and\
               (self.child is None or self.comparator(link.child, self.child)) and\
               (link.parent_position == self.parent_position or self.parent_position is None) and (
                    link.child_position == self.child_position or self.child_position is None) and (
                    link.child.anomer == self.anomer or self.anomer is None):
                yield link

    def traverse(self, structure):
        return list(self._traverse(structure))

    def apply(self, link, refund=False):
        link.break_link(refund=refund)

    def __call__(self, structure, refund=False):
        for link in self.traverse(structure):
            new_structure = structure.clone()
            link = new_structure.get_link(link.id)
            self.apply(link, refund)
            yield Glycan(
                link.parent, index_method=None).reroot(index_method=None), Glycan(
                link.child, index_method=None)


m = monosaccharides.GlcNAc
m.anomer = "beta"
alpha_1_3_mannosyl_glycoprotein_2_beta_n_acetylglucosaminyltransferase = Glycosyltransferase(
    3, 1, Anomer.alpha, parent=m, child=monosaccharides.Man)

alpha_1_6_mannosyl_glycoprotein_2_beta_n_acetylglucosaminyltransferase = Glycosyltransferase(
    6, 1, Anomer.alpha, parent=m, child=monosaccharides.Man)
