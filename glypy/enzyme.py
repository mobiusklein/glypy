import warnings
import pkg_resources
import json

from collections import defaultdict

from six import string_types as basestring

from glypy.algorithms.similarity import commutative_similarity
from glypy.algorithms import subtree_search

from glypy.utils import root as proot, StringIO

import glypy.io
from glypy.io import iupac
from glypy.io.glycoct import DistinctGlycanSet


class EnzymeCommissionNumber(object):
    def __init__(self, category, group, activity, identity):
        self.category = int(category)
        self.group = int(group)
        self.activity = int(activity)
        self.identity = int(identity)
        self._str = '.'.join(map(str, self))

    def __eq__(self, other):
        return str(self) == str(other)

    def __hash__(self):
        return hash(self._str)

    def __ne__(self, other):
        return not (self == other)

    def __getitem__(self, i):
        if i == 0:
            return self.category
        elif i == 1:
            return self.group
        elif i == 2:
            return self.activity
        elif i == 3:
            return self.identity
        elif isinstance(i, slice) or i < 0:
            return tuple(self)[i]
        else:
            raise IndexError(i)

    def __iter__(self):
        yield self.category
        yield self.group
        yield self.activity
        yield self.identity

    def __str__(self):
        return self._str

    def __repr__(self):
        return ("{self.__class__.__name__}({self.category}, {self.group},"
                " {self.activity}, {self.identity})").format(self=self)

    @classmethod
    def parse(cls, string):
        string = str(string)
        parts = string.split(".")
        if len(parts) != 4:
            raise ValueError("EC Numbers must have 4 parts (found %d): %r" % (len(parts, string)))
        try:
            parts = tuple(map(int, parts))
        except ValueError:
            raise ValueError("EC Numbers must be integers: %r" % (parts,))
        return cls(*parts)


class EnzymeInformation(object):
    def __init__(self, name, ec_number=None, alternative_names=None, **kwargs):
        if isinstance(ec_number, basestring):
            ec_number = EnzymeCommissionNumber.parse(ec_number)
        elif isinstance(ec_number, (list, tuple)):
            ec_number = EnzymeCommissionNumber(*ec_number)

        self.name = name
        self.ec_number = ec_number
        self.alternative_names = alternative_names or ()
        self.extra_information = kwargs

    def __repr__(self):
        return ("{self.__class__.__name__}({self.name!r}, {self.ec_number!r}, "
                "{self.alternative_names!r}, **{self.extra_information!r})").format(self=self)

    def __eq__(self, other):
        if isinstance(other, basestring):
            return other in (self.name, self.ec_number) or other in self.alternative_names
        else:
            try:
                return other.ec_number == self.ec_number
            except AttributeError:
                return False

    def _to_dict(self):
        store = dict()
        store["name"] = self.name
        store['ec_number'] = str(self.ec_number)
        store['alternative_names'] = tuple(self.alternative_names)
        store['extra_information'] = dict(self.extra_information)
        return store

    def __hash__(self):
        return hash(self.name)


class EnzymeDatabase(object):
    _expasy_url = "ftp://ftp.expasy.org/databases/enzyme/enzyme.dat"

    def __init__(self, fp=None, format='json'):
        self.layered_store = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
        self.direct_store = dict()

        if fp is not None:
            if format == 'expasy':
                self._parse_expasy(fp)
            elif format == 'json':
                self._parse_json(fp)
            else:
                raise ValueError("Unrecognized format: %r" % (format,))

    def add(self, enzyme_info):
        parts = enzyme_info.ec_number
        store = self.layered_store
        for k in parts[:-1]:
            store = store[k]
        store[parts[-1]] = enzyme_info
        self.direct_store[str(parts)] = enzyme_info

    def __getitem__(self, key):
        if isinstance(key[0], int):
            key = EnzymeCommissionNumber(*key)
        else:
            key = EnzymeCommissionNumber.parse(str(key))
        store = self.layered_store
        for k in key[:-1]:
            store = store[k]
        try:
            enzyme_info = store[key[-1]]
        except KeyError:
            enzyme_info = self.direct_store[str(key)]
        return enzyme_info

    def _parse_json(self, fp):
        enzymes = json.load(fp)

        for enz in enzymes:
            self.add(EnzymeInformation(**enz))

    def _dump(self, fp):
        json.dump([e._to_dict() for e in self.direct_store.values()], fp, indent=2, sort_keys=True)

    @classmethod
    def _build(cls):
        from urllib import urlopen
        fp = urlopen(cls._expasy_url)
        return cls(fp, "expasy")

    def _parse_expasy(self, fp):
        enzymes = []

        def new_store():
            return defaultdict(list)

        def postprocess(enzyme_info):
            enzyme_info['name'] = ' '.join(enzyme_info['name']).strip(".")
            enzyme_info['alternative_names'] = [
                n for n in ''.join(enzyme_info['alternative_names']).split(".") if n]
            up = enzyme_info['uniprot']
            up = [p.strip() for ent in up for p in ent.split(";") if p]
            enzyme_info['uniprot'] = up
            pro = enzyme_info['prosite']
            pro = [p.strip() for ent in pro for p in ent.split(";") if p and p != 'PROSITE']
            enzyme_info['prosite'] = pro
            enzyme_info['catalytic_activity'] = ' '.join(enzyme_info['catalytic_activity'])
            comments = enzyme_info['comments']
            comments = ' '.join(map(str.strip, map(str.rstrip, comments))).split("-!-")
            enzyme_info['comments'] = ''.join([c for c in comments if c.strip()])
            enzyme_info['ec_number'] = EnzymeCommissionNumber.parse(enzyme_info['id'])
            return enzyme_info

        current_enzyme = new_store()
        for line in fp:
            line = line.strip()
            sigil = line[:2]
            line = line[5:]
            if sigil == '//':
                if current_enzyme['id']:
                    try:
                        enzymes.append(postprocess(current_enzyme))
                    except ValueError as e:
                        if "n" in current_enzyme['id']:
                            pass
                        else:
                            print e, current_enzyme['name']
                current_enzyme = new_store()
            elif sigil == 'ID':
                current_enzyme['id'] = line
            elif sigil == 'DE':
                current_enzyme['name'].append(line.strip(" "))
            elif sigil == 'AN':
                current_enzyme['alternative_names'].append(line.strip(" "))
            elif sigil == 'CA':
                current_enzyme['catalytic_activity'].append(line.strip(" "))
            elif sigil == 'CF':
                current_enzyme['cofactors'].append(line.strip(" "))
            elif sigil == 'CC':
                current_enzyme['comments'].append(line)
            elif sigil == 'PR':
                current_enzyme['prosite'].append(line)
            elif sigil == 'DR':
                current_enzyme['uniprot'].append(line)
            else:
                print(current_enzyme.get('id'), sigil, line)

        for enz in enzymes:
            self.add(EnzymeInformation(**enz))

    @classmethod
    def _from_static(cls):
        data_buffer = pkg_resources.resource_string(glypy.io.__name__, "data/enzyme.json")
        if isinstance(data_buffer, bytes):
            data_buffer = data_buffer.decode("utf-8")
        return cls(StringIO(data_buffer), format='json')


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


class Glycoenzyme(object):
    def __init__(self, parent_position, child_position, parent, child, terminal=True,
                 identifying_information=None,
                 comparator=None,
                 validator=None):
        if comparator is None:
            comparator = commutative_similarity
        if validator is None:
            def validator(structure):
                return True
        if isinstance(identifying_information, basestring):
            identifying_information = EnzymeInformation(identifying_information)

        self._parent_position = ()
        self._child_position = ()

        self.parent_position = self._conform_position(parent_position)
        self.child_position = self._conform_position(child_position)
        self.parent = parent
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
                 parent_node_id=None, site_validator=None):
        super(Transferase, self).__init__(
            parent_position, child_position, parent, child, terminal,
            identifying_information, comparator, validator)
        if site_validator is None:
            def site_validator(structure, selected_node):
                return True
        if parent_node_id is None:
            parent_node_id = proot(self.parent).id
        self.parent_node_id = parent_node_id
        self.site_validators = self._conform_validator(site_validator)

    def _traverse(self, structure):
        for node in subtree_search.find_matching_subtree_roots(self.parent, structure, exact=True):
            node = self._get_paired_node(node)
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

    def _get_paired_node(self, node):
        return {k.id: v for k, v in subtree_search.walk_with(node, self.parent)
                }.get(self.parent_node_id)

    def apply(self, node, parent_position=None, child_position=None):
        raise NotImplementedError()

    def __call__(self, structure, parent_position=None, child_position=None):
        for node in self.traverse(structure):
            new_structure = structure.clone()
            node = new_structure.get(node.id)
            self.apply(node, parent_position, child_position)
            new_structure.reindex(hard=True)
            for new_node in new_structure:
                assert new_node.id < 1000
            yield new_structure


class Glycosyltransferase(Transferase):
    def apply(self, node, parent_position=None, child_position=None):
        new_monosaccharide = self.child.clone()
        node.add_monosaccharide(
            new_monosaccharide, position=self.parent_position[0],
            child_position=self.child_position[0])


class Substituentransferase(Glycosyltransferase):
    def apply(self, node, parent_position=None, child_position=None):
        new_substituent = self.child.clone()
        node.add_substituent(
            new_substituent, position=self.parent_position[0],
            child_position=self.child_position[0])


class Glycosylase(Glycoenzyme):

    def _test(self, link):
        return (
            (self.parent is None or self.comparator(link.parent, self.parent)),
            (self.child is None or self.comparator(link.child, self.child)),
            (self.parent_position is None or link.parent_position in self.parent_position),
            (self.child_position is None or link.child_position in self.child_position)
            (len(link.child.children()) == 0 and self.terminal) or (not self.terminal)
        )

    def _traverse(self, structure):
        for p, link in structure.iterlinks():
            if link.is_ambiguous():
                warnings.warn("Glycosidases do not support ambiguous linkages at this time.")
            else:
                if (self.parent is None or self.comparator(link.parent, self.parent)) and\
                   (self.child is None or self.comparator(link.child, self.child)) and\
                   (self.parent_position is None or link.parent_position in self.parent_position) and\
                   (self.child_position is None or link.child_position in self.child_position) and\
                   ((len(link.child.children()) == 0 and self.terminal) or (not self.terminal)):
                    yield link

    def apply(self, link, refund=False):
        link.break_link(refund=refund)

    def __call__(self, structure, refund=False):
        for link in self.traverse(structure):
            new_structure = structure.clone()
            link = new_structure.get_link(link.id)
            self.apply(link, refund)
            yield structure.__class__(
                link.parent, index_method=None).reroot(index_method='dfs'), structure.__class__(
                link.child, index_method='dfs')


def make_n_glycan_pathway():
    enzdb = EnzymeDatabase._from_static()

    bisecting_glcnac = iupac.loads("b-D-Glcp2NAc-(1-4)-b-D-Manp-(1-4)-b-D-Glcp2NAc-(1-4)-b-D-Glcp2NAc")
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
    child = iupac.loads("b-D-Glcp2NAc")

    gntE = Glycosyltransferase(3, 1, parent, child, identifying_information=enzdb[
                               2, 4, 1, 149], parent_node_id=3)

    parent = iupac.loads("b-D-Galp-(1-4)-b-D-Glcp2NAc")
    child = iupac.loads("a-D-Neup5Ac")

    siat2_6 = Glycosyltransferase(6, 2, parent, child, identifying_information=enzdb[
                                  2, 4, 99, 1], parent_node_id=3)
    siat2_3 = Glycosyltransferase(3, 2, parent, child, identifying_information=enzdb[
                                  2, 4, 99, 6], parent_node_id=3)

    parent = iupac.loads("b-D-Galp-(1-4)-b-D-Glcp2NAc")
    child = iupac.loads("a-L-Fuc")
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

    glycosylases = {
        "manI": manI,
        "manII": manII,
        "glucosidaseI": glucosidaseI
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
        "fuct6": fuct6
    }

    starting_structure_iupac = (
        "a-D-Manp-(1-2)-a-D-Manp-(1-6)-[a-D-Manp-(1-2)-a-D-Manp-(1-3)]"
        "a-D-Manp-(1-6)-[a-D-Glcp-(1-3)-a-D-Manp-(1-2)-a-D-Manp-(1-2)-a-D-Manp-(1-3)]"
        "b-D-Manp-(1-4)-b-D-Glcp2NAc-(1-4)-b-D-Glcp2NAc")
    starting_structure = iupac.loads(starting_structure_iupac)

    return glycosylases, glycosyltransferases, [starting_structure]


class Glycome(object):
    def __init__(self, glycosylases, glycosyltransferases, seeds, track_generations=False,
                 limits=None):
        if limits is None:
            limits = []
        self.glycosylases = glycosylases
        self.glycosyltransferases = glycosyltransferases
        self.seeds = seeds

        self.track_generations = track_generations
        self.enzyme_graph = defaultdict(lambda: defaultdict(set))
        self.history = []
        self.current_generation = DistinctGlycanSet(seeds)

        self.limits = limits

    def save_generation(self, generation):
        if self.track_generations:
            self.history.append(generation)

    def run(self, n=50):
        for i in range(n):
            generation = self.step()
            if not generation:
                break
            yield generation

    def within_limits(self, structure):
        for limiter in self.limits:
            if not limiter(structure):
                return False
        return True

    def step(self):
        next_generation = DistinctGlycanSet()
        for species in self.current_generation:
            parentkey = None
            for enzkey, enz in self.glycosylases.items():
                products = [root for root, leaf in enz(species, refund=1) if self.within_limits(root)]
                if products:
                    if parentkey is None:
                        parentkey = str(species)
                    for product in products:
                        childkey = str(product)
                        self.enzyme_graph[parentkey][childkey].add(enzkey)
                        next_generation.add(product)
            for enzkey, enz in self.glycosyltransferases.items():
                products = [root for root in enz(species) if self.within_limits(root)]
                if products:
                    if parentkey is None:
                        parentkey = str(species)
                    for product in products:
                        childkey = str(product)
                        self.enzyme_graph[parentkey][childkey].add(enzkey)
                        next_generation.add(product)
        # next_generation = list(next_generation.values())
        self.save_generation(self.current_generation)
        self.current_generation = next_generation
        return next_generation
