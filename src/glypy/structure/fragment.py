import itertools
import re

from glypy.composition import Composition

_fragment_shift = {
    "B": Composition(O=1, H=2),
    "Y": Composition(),
    "C": Composition(),
    "Z": Composition(H=2, O=1),
}


def _link_ids_splitter(fragment, link_ids, kind):  # pragma: no cover
    ion_types = re.findall(r"(\d+,\d+)?(\S)", kind)
    links_broken = link_ids

    pairings = zip(ion_types, links_broken)

    fragment.link_ids = {
        link_id: ion_type for ion_type,
        link_id in pairings if ion_type[0] == ""}
    fragment.crossring_cleavages = {
        node_id: ion_type for ion_type,
        node_id in pairings if ion_type[0] != ""}


latex_symbol_map = {
    "a": r"\alpha",
    "b": r"\beta",
    "c": r"\sigma",
    "d": r"\delta",
    "e": r"\epsilon",
    "f": r"\zeta",
    "g": r"\gamma",
    "h": r"\eta",
    "i": r"\iota",
    "j": r"\lambda",
    "k": r"\kappa",
    "l": r"\lambda",
    "m": r"\mu",
    "n": r"\nu",
    "o": r"\omicron",
    "p": r"\pi",
    "q": r"\sigma",
    "r": r"\rho",
    "s": r"\upsilon",
    "t": r"\tau",
    "u": r"\upsilon",
    "v": r"\omega",
}


class GlycanFragment(object):
    '''
    A simple container for a fragment ion, produced by :meth:`Glycan.fragments`

    Attributes
    ----------
    kind: |str|
        One of A, B, C, X, Y, or Z for each link broken or ring cleaved

    link_ids: |dict| of |int| -> |tuple|
        The :attr:`id` value of each link cleaved to the corresponding cleavage type

    included_nodes: |list| of |int|
        The :attr:`id` value of each |Monosaccharide| contained in the fragment

    mass: |float|
        The mass of the fragment

    name: |str|
        The fragment name under the branching nomenclature

    crossring_cleavages: |dict| of |int| -> |tuple|
        The :attr:`id` value of each link cleaved to the corresponding cleavage type, including ring coordinates

    composition: |Composition|
        The elemental composition of the fragment

    See Also
    --------
    :meth:`Glycan.fragments`
    '''
    __slots__ = [
        "mass",
        "kind",
        "included_nodes",
        "link_ids",
        "name",
        "crossring_cleavages",
        "composition",
    ]

    def __init__(self, kind, link_ids, included_nodes, mass,
                 name=None, crossring_cleavages=None, composition=None):
        self.mass = mass
        self.kind = kind
        self.link_ids = link_ids
        self.included_nodes = included_nodes
        self.crossring_cleavages = crossring_cleavages
        self.name = name
        self.composition = composition

    def copy(self):
        """Create a deep copy of the fragment.

        Returns
        -------
        GlycanFragment
        """
        return self.__class__(
            self.kind, self.link_ids.copy(), self.included_nodes.copy(), self.mass,
            self.name, self.crossring_cleavages.copy(), self.composition.copy())

    def clone(self):
        """Create a deep copy of the fragment.

        Returns
        -------
        GlycanFragment
        """
        return self.copy()

    def is_reducing(self):
        """Is this fragment from the reducing end

        Returns
        -------
        |bool|
        """
        return bool(set(self.kind) & set("XYZ"))

    def is_non_reducing(self):
        """Is this fragment from the non-reducing end

        Returns
        -------
        |bool|
        """
        return bool(set(self.kind) & set("ABC"))

    def is_internal(self):
        """Is the fragment composed of both reducing and non-reducing bond
        cleavage events.

        Returns
        -------
        bool
        """
        return bool(self.is_reducing() and self.is_non_reducing())

    @property
    def fname(self):
        """A `LaTeX` formatted version of the fragment's name

        Returns
        -------
        str
        """
        buff = []
        for c in self.name:
            if c in latex_symbol_map:
                buff.append("$_{}$".format(latex_symbol_map[c]))
            else:
                buff.append(c)
        return ''.join(buff)

    def _make_tex_name(self):
        buff = []
        for c in self.name:
            if c in latex_symbol_map:
                buff.append("_{}".format(latex_symbol_map[c]))
            else:
                buff.append(c)
        return ''.join(buff)

    @property
    def series(self):
        """An alias for :attr:`kind`

        Returns
        -------
        str
        """
        return self.kind

    def __getstate__(self):  # pragma: no cover
        d = {}
        for a in self.__slots__:
            d[a] = getattr(self, a)
        return d

    def __setstate__(self, state):  # pragma: no cover
        for a, v in state.items():
            setattr(self, a, v)

    def __hash__(self):
        return hash((self.name, int(self.mass)))

    def __eq__(self, other):  # pragma: no cover
        for field in self.__slots__:
            if getattr(self, field) != getattr(other, field, NotImplemented):
                return False
        return True

    def __ne__(self, other):  # pragma: no cover
        return not self == other

    def __repr__(self):  # pragma: no cover
        rep = "<GlycanFragment "
        for f in self.__slots__:
            rep += " {}={}".format(f, getattr(self, f))
        rep += ">"
        return rep

    @property
    def __dict__(self):
        return self.__getstate__()

    @property
    def break_count(self):
        """The number of cleavage events to create this fragment.

        Returns
        -------
        int
        """
        return len(self.link_ids) + len(self.crossring_cleavages)

    @property
    def residues_contained(self):
        """The number of :class:`~.Monosaccharide` nodes included in the fragmented sub-tree.

        Returns
        -------
        int
        """
        return len(self.included_nodes)


Fragment = GlycanFragment


class GlycanSubstructure(object):
    """Represent a subgraph of a :class:`~.Glycan` produced
    from another :class:`~.Glycan` by breaking linkages.

    Supports the :func:`~.root` and :func:`~.tree` protocols.

    Attributes
    ----------
    tree : :class:`~.Glycan`
        The substructure
    include_nodes : set
        The :attr:`~.Monosaccharide.id` of all included monosaccharides
    link_ids : list
        The :attr:`~.Link.id` of all linkages broken to form this substructure
    parent_breaks : dict
        Mapping from :attr:`~.Link.id` to :attr:`~.Monosaccharide.id` if that monosaccharide
        was the :attr:`~.Link.parent` of that linkage and that linkage was broken.
    child_breaks : dict
        Mapping from :attr:`~.Link.id` to :attr:`~.Monosaccharide.id` if that monosaccharide
        was the :attr:`~.Link.child` of that linkage and that linkage was broken.
    crossring_cleavages : dict
        Mapping from :attr:`~.Monosaccharide.id` to :class:`str` denoting the ring coordinate pair that
        was cleaved to produce any included crossring cleavage events.

    """

    def __init__(self, tree, include_nodes, link_ids,
                 parent_breaks, child_breaks, crossring_cleavages=None):
        self.tree = tree
        self.include_nodes = include_nodes
        self.link_ids = link_ids
        self.parent_breaks = parent_breaks
        self.child_breaks = child_breaks
        self.crossring_cleavages = crossring_cleavages or {}

    def contains_reducing_end(self):
        """Whether or not the sub-tree includes a :class:`ReducedEnd`

        Returns
        -------
        bool
        """
        return self.tree.root.reducing_end is not None

    def to_fragments(self, kind="BY", average=False, charge=None, mass_data=None, include_composition=True,
                     traversal_method='dfs'):
        """Transform an instance of :class:`GlycanSubstructure` into every combination of
        :class:`GlycanFragment` allowed under `kind`.

        Parameters
        ----------
        kind : Iterable, optional
            The types of fragments to emit. Defaults to "BY"
        average : bool, optional
            Calculate masses with average isotopic composition
        charge : int, optional
            Calculate `m/z` instead of neutral mass, with `z = charge`
        mass_data : dict, optional
            If mass_data is None, standard NIST mass and isotopic abundance data are used. Otherwise the
            contents of mass_data are assumed to contain elemental mass and isotopic abundance information.
            Defaults to :const:`None`.
        include_composition: bool, optional
            Whether or not to populate the `composition` attribute of the fragment. Defaults to :const:`True`

        Yields
        -------
        :class:`GlycanFragment`
        """
        parent_type = set("YZ") & set(kind)
        child_type = set("BC") & set(kind)
        crossring_type = set("AX") & set(kind)
        parent_shifts = [(lid, parent_type)
                         for lid, nid in self.parent_breaks.items()]
        child_shifts = [(lid, child_type)
                        for lid, nid in self.child_breaks.items()]

        shift_masses = {k: _fragment_shift[k].calc_mass(average=average,
                                                        charge=charge,
                                                        mass_data=mass_data)
                        for k in parent_type | child_type}
        crossring_contained = {
            kind for pos,
            kind in self.crossring_cleavages.values()}
        # Do not proceed if crossring fragments are included and
        # were not requested.
        if len(crossring_contained & crossring_type) == 0 and len(
                crossring_contained) != 0:
            # raise StopIteration()
            return

        all_shifts = parent_shifts + child_shifts
        all_link_ids = [i for i, t in all_shifts]
        frag_types = [t for i, t in all_shifts]
        base_mass = self.tree.mass(
            average=average,
            charge=charge,
            mass_data=mass_data, method=traversal_method)
        if include_composition:
            base_composition = self.tree.total_composition(method=traversal_method)
        # product of splat of empty list is a list of the empty list. So a fragment with
        # no glycosidic cleavages still enters this outer loop, letting only crossring-cleavage
        # GlycanSubstructure instances through without issue
        for shift_set in itertools.product(*frag_types):
            mass_offset = 0.0
            composition_offset = Composition()
            link_ids = {}
            # The type of fragment being produced, expressed a collection of ABCXYZs
            kind = [] + [''.join(cr_kind)
                         for cr_kind in self.crossring_cleavages.values()]
            i = 0
            shift_set = list(shift_set)
            for shift in shift_set:
                link_id = all_link_ids[i]
                shift = shift[0]
                mass_offset -= shift_masses[shift]
                composition_offset -= _fragment_shift[shift]
                link_ids[link_id] = ("", shift)
                kind.append(shift)
                i += 1

            if include_composition:
                fragment_composition = base_composition + composition_offset

            yield GlycanFragment(kind=''.join(kind), link_ids=link_ids, included_nodes=self.include_nodes,
                                 mass=base_mass + mass_offset, name=None,
                                 crossring_cleavages=self.crossring_cleavages,
                                 composition=fragment_composition)

    def __repr__(self):  # pragma: no cover
        rep = ("<GlycanSubstructure include_nodes={} link_ids={} parent_breaks={}"
               " child_breaks={} crossring_cleavages={}>\n{}").format(
                   self.include_nodes, self.link_ids, self.parent_breaks,
                   self.child_breaks, self.crossring_cleavages, self.tree)
        return rep

    def __root__(self):  # pragma: no cover
        return self.tree.root

    def __tree__(self):  # pragma: no cover
        return self.tree


Subtree = GlycanSubstructure
