import re


def link_ids_splitter(fragment, link_ids, kind):
    ion_types = re.findall(r"(\d+,\d+)?(\S)", kind)
    links_broken = link_ids

    pairings = zip(ion_types, links_broken)

    fragment.link_ids = {link_id: ion_type for ion_type, link_id in pairings if ion_type[0] == ""}
    fragment.crossring_cleavages = {node_id: ion_type for ion_type, node_id in pairings if ion_type[0] != ""}


class Fragment(object):
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
        The mass or `m/z` of the fragment

    name: |str|
        The fragment name under the branching nomenclature

    crossring_cleavages: |dict| of |int| -> |tuple|
        The :attr:`id` value of each link cleaved to the corresponding cleavage type, including ring coordinates

    score: |float|
        A score value assigned to the fragment structure by an application

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
        "score"
    ]

    def __init__(self, kind, link_ids, included_nodes, mass,
                 name=None, crossring_cleavages=None, score=0.0):
        self.mass = mass
        self.kind = kind
        self.link_ids = link_ids
        self.included_nodes = included_nodes
        self.crossring_cleavages = crossring_cleavages
        self.name = name
        self.score = score
        # The data fed from :meth:`Glycan.fragments` lumps cross ring cleavages
        # together with link ids. Separate them out.
        if crossring_cleavages is None:
            link_ids_splitter(self, self.link_ids, self.kind)

    def is_reducing(self):
        """Is this fragment from the reducing end

        Returns
        -------
        |bool|
        """
        return set(self.kind) in set("XYZ")

    def is_non_reducing(self):
        """Is this fragment from the non-reducing end

        Returns
        -------
        |bool|
        """
        return set(self.kind) in set("ABC")

    def __getstate__(self):
        d = {}
        for a in self.__slots__:
            d[a] = getattr(self, a)
        return d

    def __setstate__(self, state):
        if isinstance(state, tuple):
            kind, link_ids, included_nodes, mass, name = state
            link_ids_splitter(self, link_ids, kind)
            self.kind = kind
            self.included_nodes = included_nodes
            self.mass = mass
            self.name = name
            self.score = 0.0
        else:
            for a, v in state.items():
                setattr(self, a, v)

    def __eq__(self, other):
        for field in self.__slots__:
            if getattr(self, field) != getattr(other, field, NotImplemented):
                return False
        return True

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        rep = "<Fragment "
        for f in self.__slots__:
            rep += " {}={}".format(f, getattr(self, f))
        rep += ">"
        return rep
