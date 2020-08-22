import warnings

from six import string_types as basestring

import glypy

from glypy.io import iupac

from glypy.structure.base import MoleculeBase
from glypy.algorithms import subtree_search
from glypy.algorithms.similarity import commutative_similarity

from glypy.utils import root as proot

from .ec import EnzymeInformation, EnzymeDatabase


def rejecting(*args):
    def checker(structure):
        for subtree in args:
            if subtree_search.subtree_of(subtree, structure, exact=True):
                return False
        return True
    return checker


def reject_on_path(*args):
    def checker(structure, selected_node):
        for subtree in args:
            path = {
                v.id: k for k, v in subtree_search.walk_with(structure, subtree)
            }
            if selected_node.id in path:
                return False
        return True
    return checker


def reject_on_parent(*args):
    def checker(structure, selected_node):
        for t in args:
            for _, parent in selected_node.parents():
                if parent == t:
                    return False
        return True
    return checker


class Glycoenzyme(object):

    def __init__(self, parent_position, child_position, parent, child, terminal=True,
                 identifying_information=None,
                 comparator=None,
                 validator=None):
        if comparator is None:
            comparator = commutative_similarity
        if validator is None:
            def validator(structure):  # pylint: disable=function-redefined
                return True
        if isinstance(identifying_information, basestring):
            identifying_information = EnzymeInformation(
                identifying_information)

        self._parent_position = ()
        self._child_position = ()
        self._parents = ()
        self._child = None

        self.parent_position = parent_position
        self.child_position = child_position
        self.parents = parent
        self.child = child
        self.terminal = terminal
        self.comparator = comparator
        self.validators = self._conform_validator(validator)
        self.identifying_information = identifying_information

    @property
    def parent_position(self):
        return self._parent_position

    @parent_position.setter
    def parent_position(self, value):
        self._parent_position = self._conform_position(value)

    @property
    def child_position(self):
        return self._child_position

    @child_position.setter
    def child_position(self, value):
        self._child_position = self._conform_position(value)

    @property
    def parents(self):
        return self._parent

    @parents.setter
    def parents(self, value):
        self._parent = self._conform_molecule(value)

    @property
    def child(self):
        return self._child

    @child.setter
    def child(self, value):
        self._child = value

    def _conform_validator(self, fn):
        try:
            iter(fn)
        except TypeError:
            fn = (fn,)
        return fn

    def validate_structure(self, structure):
        for fn in self.validators:
            if not fn(structure):
                return False
        return True

    def _conform_position(self, value):
        if value is None:
            return None
        try:
            return tuple(value)
        except TypeError:
            return (value,)

    def _conform_molecule(self, value):
        if value is None:
            return value
        if isinstance(value, MoleculeBase):
            value = (value, )
        else:
            value = tuple(value)
        return value

    def _traverse(self, structure):
        raise NotImplementedError()

    def traverse(self, structure):
        if self.validate_structure(structure):
            return list(self._traverse(structure))
        else:
            return []


class Transferase(Glycoenzyme):

    def __init__(self, parent_position, child_position, parent, child, terminal=True,
                 identifying_information=None, comparator=None, validator=None,
                 parent_node_id=None, site_validator=None, exact=True):
        super(Transferase, self).__init__(
            parent_position, child_position, parent, child, terminal,
            identifying_information, comparator, validator)
        if site_validator is None:
            def site_validator(structure, selected_node):  # pylint: disable=function-redefined
                return True
        if parent_node_id is None:
            parent_node_id = proot(self.parent).id
        self.parent_node_id = parent_node_id
        self.site_validators = self._conform_validator(site_validator)
        self.exact = exact

    def _traverse(self, structure):
        if self.parents is not None:
            for parent in self.parents:
                for node in subtree_search.find_matching_subtree_roots(parent, structure, exact=True):
                    node = self._get_paired_node(node, parent)
                    for parent_position in self.parent_position:
                        if not node.is_occupied(parent_position) and (
                           (len(node.children()) == 0 and self.terminal) or (
                                not self.terminal)) and self.validate_site(structure, node):
                            yield node

    def validate_site(self, structure, node):
        for fn in self.site_validators:
            if not fn(structure, node):
                return False
        return True

    def _get_paired_node(self, node, parent):
        return {k.id: v for k, v in subtree_search.walk_with(node, parent)
                }.get(self.parent_node_id)

    def make_bond(self, node, parent_position=None, child_position=None):
        raise NotImplementedError()

    def apply(self, structure, parent_position=None, child_position=None):
        for node in self.traverse(structure):
            new_structure = structure.clone()
            node = new_structure.get(node.id)
            self.make_bond(node, parent_position, child_position)
            # reset ids, rebuild the index, and standardize traversal order
            new_structure.reindex()
            new_structure.canonicalize()
            yield new_structure

    def __call__(self, structure, parent_position=None, child_position=None):
        return self.apply(
            structure, parent_position=parent_position, child_position=child_position)


class Glycosyltransferase(Transferase):

    def make_bond(self, node, parent_position=None, child_position=None):
        new_monosaccharide = self.child.clone()
        node.add_monosaccharide(
            new_monosaccharide, position=self.parent_position[0],
            child_position=self.child_position[0])


class Substituentransferase(Glycosyltransferase):

    def make_bond(self, node, parent_position=None, child_position=None):
        new_substituent = self.child.clone()
        node.add_substituent(
            new_substituent, position=self.parent_position[0],
            child_position=self.child_position[0])


class Glycosylase(Glycoenzyme):

    @property
    def child(self):
        return self._child

    @child.setter
    def child(self, value):
        self._child = self._conform_molecule(value)

    def _test(self, link):
        return (
            (self.parents is None or any(self.comparator(link.parent, parent) for parent in self.parents)),
            (self.child is None or any(self.comparator(link.child, child) for child in self.child)),
            (self.parent_position is None or link.parent_position in self.parent_position),
            (self.child_position is None or link.child_position in self.child_position)
            (len(link.child.children()) == 0 and self.terminal) or (not self.terminal)
        )

    def _traverse(self, structure):
        for p, link in structure.iterlinks():
            if link.is_ambiguous():
                warnings.warn(
                    "Glycosylase do not support ambiguous linkages at this time.")
            else:
                if (self.parents is None or any(self.comparator(link.parent, parent) for parent in self.parents)) and\
                   (self.child is None or any(self.comparator(link.child, child) for child in self.child)) and\
                   (self.parent_position is None or link.parent_position in self.parent_position) and\
                   (self.child_position is None or link.child_position in self.child_position) and\
                   ((len(link.child.children()) == 0 and self.terminal) or (not self.terminal)):
                    yield link

    def digest(self, link, refund=False):
        link.break_link(refund=refund)

    def apply(self, structure, refund=False):
        for link in self.traverse(structure):
            new_structure = structure.clone()
            link = new_structure.get_link(link.id)
            self.digest(link, refund)
            parent, child = structure.__class__(
                link.parent, index_method=None).reroot(index_method='dfs'), structure.__class__(
                link.child, index_method='dfs')
            yield parent.canonicalize(), child.canonicalize()

    def __call__(self, structure, refund=False):
        return self.apply(structure, refund=refund)


def make_n_glycan_pathway():
    enzdb = EnzymeDatabase._from_static()

    bisecting_glcnac = iupac.loads(
        "b-D-Glcp2NAc-(1-4)-b-D-Manp-(1-4)-b-D-Glcp2NAc-(1-4)-b-D-Glcp2NAc")
    galactose = iupac.loads("b-D-Galp")

    parent = iupac.loads(
        "a-D-Manp-(1-3)-b-D-Manp-(1-4)-b-D-Glcp2NAc-(1-4)-b-D-Glcp2NAc")
    child = iupac.loads("b-D-Glcp2NAc")
    gntI = Glycosyltransferase(2, 1, parent, child, terminal=1,
                               identifying_information=enzdb[2, 4, 1, 101], parent_node_id=6)

    parent = iupac.loads(
        "a-D-Manp-(1-6)-[b-D-Glcp2NAc-(1-2)-a-D-Manp-(1-3)]b-D-Manp-(1-4)-b-D-Glcp2NAc-(1-4)-b-D-Glcp2NAc")
    child = iupac.loads("b-D-Glcp2NAc")
    gntII = Glycosyltransferase(2, 1, parent, child, terminal=True,
                                identifying_information=enzdb[2, 4, 1, 143],
                                parent_node_id=6,
                                validator=rejecting(bisecting_glcnac, galactose))

    parent = iupac.loads(
        "beta-D-Glcp2NAc-(1->2)-alpha-D-Manp-(1->3)-"
        "[alpha-D-Manp-(1->6)]"
        "-beta-D-Manp-(1->4)-beta-D-Glcp2NAc-(1->4)-beta-D-Glcp2NAc")
    child = iupac.loads("b-D-Glcp2NAc")
    gntIII = Glycosyltransferase(4, 1, parent, child, terminal=False, parent_node_id=5,
                                 identifying_information=enzdb[2, 4, 1, 144],
                                 validator=rejecting(bisecting_glcnac))

    parent = iupac.loads(
        "b-D-Glcp2NAc-(1-2)-a-D-Manp-(1-3)-b-D-Manp-(1-4)-b-D-Glcp2NAc-(1-4)-b-D-Glcp2NAc")
    child = iupac.loads("b-D-Glcp2NAc")
    gntIV = Glycosyltransferase(4, 1, parent, child, terminal=False, parent_node_id=6,
                                identifying_information=enzdb[2, 4, 1, 145],
                                validator=rejecting(bisecting_glcnac, galactose))

    parent = iupac.loads(
        "b-D-Glcp2NAc-(1-2)-a-D-Manp-(1-6)-b-D-Manp-(1-4)-b-D-Glcp2NAc-(1-4)-b-D-Glcp2NAc")
    child = iupac.loads("b-D-Glcp2NAc")
    gntV = Glycosyltransferase(6, 1, parent, child, terminal=False, parent_node_id=6,
                               identifying_information=enzdb[2, 4, 1, 155],
                               validator=rejecting(bisecting_glcnac, galactose))

    parent = iupac.loads(
        "b-D-Glcp2NAc-(1-6)-[b-D-Glcp2NAc-(1-2)]"
        "a-D-Manp-(1-6)-[a-D-Manp-(1-3)]b-D-Manp-(1-4)-b-D-Glcp2NAc-(1-4)-b-D-Glcp2NAc")
    child = iupac.loads("b-D-Glcp2NAc")
    gntVI = Glycosyltransferase(4, 1, parent, child, terminal=False, parent_node_id=6,
                                identifying_information=enzdb[2, 4, 1, 145],
                                validator=rejecting(bisecting_glcnac, galactose))

    parent = iupac.loads("b-D-Glcp2NAc")
    child = iupac.loads("b-D-Galp")

    galt = Glycosyltransferase(4, 1, parent, child, identifying_information=enzdb[
                               2, 4, 1, 38], parent_node_id=1,
                               site_validator=reject_on_path(bisecting_glcnac))

    parent = iupac.loads("b-D-Galp-(1-4)-b-D-Glcp2NAc")
    child = iupac.loads("a-D-Galp")

    agal13galt = Glycosyltransferase(3, 1, parent, child, identifying_information=None,
                                     parent_node_id=3)

    parent = iupac.loads("b-D-Galp-(1-4)-b-D-Glcp2NAc")
    child = iupac.loads("b-D-Glcp2NAc")

    gntE = Glycosyltransferase(3, 1, parent, child, identifying_information=enzdb[
                               2, 4, 1, 149], parent_node_id=3)

    parent = iupac.loads("b-D-Galp-(1-4)-b-D-Glcp2NAc")
    child = iupac.loads("a-D-Neup5Ac")

    siat2_6 = Glycosyltransferase(6, 2, parent, child, identifying_information=enzdb[
                                  2, 4, 99, 1], parent_node_id=3)
    siat2_3 = Glycosyltransferase(3, 2, parent, child, identifying_information=enzdb[
                                  2, 4, 99, 6], parent_node_id=3)
    # siat2_8 = Glycosyltransferase(
    #     8, 2, child.clone(), child, identifying_information=enzdb[2, 4, 99, 8],
    #     site_validator=reject_on_parent(child.clone()), parent_node_id=1)

    # parent = iupac.loads("b-D-Glcp2NAc")
    # child = iupac.loads("b-D-Galp2NAc")
    # b4galnact = Glycosyltransferase(4, 1, parent, child, parent_node_id=1, identifying_information=None)

    parent = iupac.loads("b-D-Galp-(1-4)-b-D-Glcp2NAc")
    child = iupac.loads("a-L-Fucp")
    fuct3 = Glycosyltransferase(3, 1, parent, child, terminal=False,
                                parent_node_id=1, identifying_information=enzdb[2, 4, 1, 152])

    parent = iupac.loads(
        "beta-D-Glcp2NAc-(1->2)-alpha-D-Manp-(1->3)-[beta-D-Glcp2NAc-(1->2)-alpha-D-Manp-(1->6)]"
        "-beta-D-Manp-(1->4)-beta-D-Glcp2NAc-(1->4)-beta-D-Glcp2NAc")
    child = iupac.loads("a-L-Fucp")
    fuct6 = Glycosyltransferase(6, 1, parent, child, terminal=False,
                                parent_node_id=1, identifying_information=enzdb[2, 4, 1, 68])

    manI = Glycosylase((2,), 1, iupac.loads("a-D-Manp"), iupac.loads("a-D-Manp"),
                       identifying_information=enzdb["3.2.1.113"])
    manII = Glycosylase((3, 6), 1, iupac.loads("a-D-Manp"), iupac.loads("a-D-Manp"),
                        identifying_information=enzdb["3.2.1.114"],
                        validator=rejecting(bisecting_glcnac))
    glucosidaseI = Glycosylase(None, None, iupac.loads("alpha-D-Manp"), iupac.loads("?-?-Glcp"),
                               identifying_information=enzdb["3.2.1.106"])
    glucosidaseII = Glycosylase((3,), 1, iupac.loads("alpha-D-Glcp"), iupac.loads("alpha-D-Glcp"),
                                identifying_information=enzdb['3.2.1.84'])

    glycosylases = {
        "manI": manI,
        "manII": manII,
        "glucosidaseI": glucosidaseI,
        "glucosidaseII": glucosidaseII,
    }

    glycosyltransferases = {
        "gntI": gntI,
        "gntII": gntII,
        "gntIII": gntIII,
        "gntIV": gntIV,
        "gntV": gntV,
        "gntVI": gntVI,
        "galt": galt,
        "gntE": gntE,
        "siat2_6": siat2_6,
        "siat2_3": siat2_3,
        "fuct3": fuct3,
        "fuct6": fuct6,
        "agal13galt": agal13galt,
        # pending literature search: https://www.ncbi.nlm.nih.gov/pubmed/7881179
        # "b4galnact": b4galnact,
    }

    starting_structure_iupac = (
        "a-D-Manp-(1-2)-a-D-Manp-(1-6)-[a-D-Manp-(1-2)-a-D-Manp-(1-3)]"
        "a-D-Manp-(1-6)-[a-D-Glcp-(1-3)-a-D-Glcp-(1-3)-a-D-Glcp-(1-3)-"
        "a-D-Manp-(1-2)-a-D-Manp-(1-2)-a-D-Manp-(1-3)]"
        "b-D-Manp-(1-4)-b-D-Glcp2NAc-(1-4)-b-D-Glcp2NAc")
    starting_structure = iupac.loads(starting_structure_iupac)

    return glycosylases, glycosyltransferases, [starting_structure]


def make_mucin_type_o_glycan_pathway():
    enzdb = EnzymeDatabase._from_static()

    parent = iupac.loads("a-D-Galp2NAc")
    child = iupac.loads("b-D-Galp")
    c1galt1 = Glycosyltransferase(
        3, 1, parent, child, identifying_information=enzdb[2, 4, 1, 122])

    parent = iupac.loads("a-D-Galp2NAc")
    child = iupac.loads("b-D-Glcp2NAc")
    b3gnt6 = Glycosyltransferase(3, 1, parent, child, terminal=False,
                                 identifying_information=enzdb[2, 4, 1, 147])

    parent = iupac.loads("a-D-Galp2NAc")
    child = iupac.loads("b-D-Glcp2NAc")
    gcnt1 = Glycosyltransferase(6, 1, parent, child, terminal=False,
                                identifying_information=enzdb[2, 4, 1, 102])

    parent = iupac.loads('b-D-Glcp2NAc-(1-6)-a-D-Galp2NAc')
    child = iupac.loads("b-D-Galp")
    b4galt5 = Glycosyltransferase(4, 1, parent, child, parent_node_id=3)

    parent = iupac.loads('b-D-Glcp2NAc-(1-3)-a-D-Galp2NAc')
    child = iupac.loads('b-D-Glcp2NAc')
    gcnt3 = Glycosyltransferase(6, 1, parent, child, parent_node_id=1, terminal=False,
                                identifying_information=enzdb[2, 4, 1, 148])

    parent = iupac.loads("a-D-Galp2NAc")
    child = iupac.loads("a-D-Galp2NAc")
    core7_galnact = Glycosyltransferase(6, 1, parent, child)

    parent = iupac.loads("a-D-Galp2NAc")
    child = iupac.loads("a-D-Galp")
    core8_galt = Glycosyltransferase(3, 1, parent, child)

    parent = iupac.loads("a-D-Galp2NAc")
    child = iupac.loads("a-D-Neup5Ac")
    st6gal1 = Glycosyltransferase(6, 2, parent, child, terminal=False,
                                  identifying_information=enzdb[2, 4, 99, 3])

    parent = iupac.loads("b-D-Galp-(1-3)-a-D-Galp2NAc")
    child = iupac.loads("a-D-Neup5Ac")
    st3gal2 = Glycosyltransferase(3, 2, parent, child, parent_node_id=3)

    parent = iupac.loads("b-D-Galp-(1-3)-a-D-Galp2NAc")
    child = iupac.loads("a-D-Neup5Ac")
    st6gal2 = Glycosyltransferase(6, 2, parent, child, parent_node_id=3)

    parent = iupac.loads('b-D-Galp-(1-4)-b-D-Glcp2NAc')
    child = iupac.loads("a-L-Fucp")
    fuct2 = Glycosyltransferase(2, 1, parent, child, terminal=False, parent_node_id=3,
                                identifying_information=enzdb['2.4.1.69'])

    parent = iupac.loads("b-D-Galp-(1-4)-b-D-Glcp2NAc")
    child = iupac.loads("a-L-Fucp")
    fuct3 = Glycosyltransferase(3, 1, parent, child, terminal=False,
                                parent_node_id=1, identifying_information=enzdb[2, 4, 1, 152])

    parent = iupac.loads("b-D-Galp-(1-4)-b-D-Glcp2NAc")
    child = iupac.loads("a-D-Neup5Ac")
    siat2_6 = Glycosyltransferase(6, 2, parent, child, identifying_information=enzdb[
                                  2, 4, 99, 1], parent_node_id=3)

    parent = iupac.loads("b-D-Glcp2NAc")
    child = iupac.loads("b-D-Galp")
    galt4 = Glycosyltransferase(4, 1, parent, child, identifying_information=enzdb[
        2, 4, 1, 86])

    parent = iupac.loads("b-D-Galp-(1-4)-b-D-Glcp2NAc")
    child = iupac.loads("b-D-Glcp2NAc")
    gntE = Glycosyltransferase(3, 1, parent, child, identifying_information=enzdb[
        2, 4, 1, 149], parent_node_id=3)

    parent = iupac.loads("b-D-Galp-(1-4)-b-D-Glcp2NAc-(1-3)-b-D-Galp")
    child = iupac.loads("b-D-Glcp2NAc")
    gntII = Glycosyltransferase(6, 1, parent, child, terminal=False,
                                identifying_information=enzdb[2, 4, 1, 150], parent_node_id=1)

    parent = iupac.loads("b-D-Galp-(1-3)-a-D-Galp2NAc")
    child = iupac.loads("b-D-Glcp2NAc")
    b3gnt3 = Glycosyltransferase(3, 1, parent, child, identifying_information=enzdb[
        "2.4.1.102"], parent_node_id=3, terminal=False)

    glycosyltransferases = {
        "c1galt1": c1galt1,
        "b3gnt6": b3gnt6,
        "gcnt1": gcnt1,
        "b4galt5": b4galt5,
        "gcnt3": gcnt3,
        "core7_galnact": core7_galnact,
        "core8_galt": core8_galt,
        "st6gal1": st6gal1,
        "st3gal2": st3gal2,
        "st6gal2": st6gal2,
        "fuct2": fuct2,
        "fuct3": fuct3,
        "siat2_6": siat2_6,
        "galt4": galt4,
        "gntE": gntE,
        "gntII": gntII,
        "b3gnt3": b3gnt3
    }

    glycosylases = {}

    seeds = [glypy.Glycan(iupac.loads("a-D-Galp2NAc"))]
    return glycosylases, glycosyltransferases, seeds


ABS_Sialidase = Glycosylase((3, 6, 8, 9), 2, parent=None, child=iupac.loads("a-D-Neup5Ac"))
NAN1_Sialidase = Glycosylase((3,), 2, parent=None, child=iupac.loads("a-D-Neup5Ac"))
BKF_Fucosidase = Glycosylase((2, 3, 4, 6), 1, parent=None, child=iupac.loads("a-L-Fucp"))
XMF_Fucosidase = Glycosylase((2,), 1, parent=None, child=iupac.loads("a-L-Fucp"))
AMF_Fucosidase = Glycosylase((3, 4), 1, parent=None, child=iupac.loads("a-L-Fucp"))
BTG_Galactosidase = Glycosylase((3, 4), 1, parent=None, child=iupac.loads("b-D-Galp"))
SPG_Galactosidase = Glycosylase((4,), 1, parent=None, child=iupac.loads("b-D-Galp"))
CBG_Galactosidase = Glycosylase((3, 4, 6), 1, parent=None, child=iupac.loads("b-D-Galp"))
JBM_Mannosidase = Glycosylase((3, 4, 6), 1, parent=None, child=iupac.loads("a-D-Manp"))

GUH_N_Acetylhexosaminidase = Glycosylase(
    (2, 3, 4, 6), 1, parent=None, child=iupac.loads("b-D-GlcpNAc"),
    # Add site_validator-like option to Glycosylase
    # Block action on bisecting GlcNAc
    # site_validator=reject_on_path(
    #     iupac.loads("b-D-Glcp2NAc-(1-4)-b-D-Manp-(1-4)-b-D-Glcp2NAc-(1-4)-b-D-Glcp2NAc"))
)

JBH_N_Acetylhexosaminidase = Glycosylase((2, 3, 4, 6), 1, parent=None, child=(
    iupac.loads("b-D-GlcpNAc"),
    iupac.loads("b-D-GalpNAc")))

glycodigest_rules = {
    "ABS": ABS_Sialidase,
    "NAN1": NAN1_Sialidase,
    "BKF": BKF_Fucosidase,
    "XMF": XMF_Fucosidase,
    "AMF": AMF_Fucosidase,
    "BTG": BTG_Galactosidase,
    "SPG": SPG_Galactosidase,
    "CBG": CBG_Galactosidase,
    "JBM": JBM_Mannosidase,
    "GUH": GUH_N_Acetylhexosaminidase,
    "JBH": JBH_N_Acetylhexosaminidase
}
