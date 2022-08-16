
def reduced_end_compat(reduced_end):
    if reduced_end is None:
        return
    if not hasattr(reduced_end, "base_composition"):
        if len(reduced_end.links) == 0:
            reduced_end.base_composition = reduced_end.composition.clone()
        else:
            reduced_end.base_composition = reduced_end.composition.clone()
            for pos, link in reduced_end.links.items():
                reduced_end.base_composition += link.parent_loss
    if reduced_end.valence != 1:
        try:
            reduced_end.drop_substituent(2)
        except IndexError:
            pass
        reduced_end.valence = 1


#: Copies of the functions of the same name instrumented for determining where exactly an
#: equality search failed.
def exact_ordering_equality(self, other, substituents=True, visited=None):
    '''
    Performs equality testing between two monosaccharides where
    the exact position (and ordering by sort) of links must to match between
    the input |Monosaccharide| objects

    Returns
    -------
    |bool|
    '''
    if visited is None:
        visited = set()
    if (self.id, other.id) in visited:
        return True

    visited.add((self.id, other.id))

    if self._flat_equality(other):
        if substituents:
            for a_sub, b_sub in zip(self.substituents(), other.substituents()):
                if a_sub != b_sub:
                    return False
        for a_mod, b_mod in zip(self.modifications.items(), other.modifications.items()):
            if a_mod != b_mod:
                return False
        for a_child, b_child in zip(self.children(), other.children()):
            if a_child[0] != b_child[0]:
                return False
            if not a_child[1].exact_ordering_equality(b_child[1],
                                                      substituents=substituents,
                                                      visited=visited):
                return False
        return True
    return False


def topological_equality(self, other, substituents=True, visited=None):
    '''
    Performs equality testing between two monosaccharides where
    the exact ordering of child links does not have to match between
    the input |Monosaccharide|s, so long as an exact match of the
    subtrees is found

    Returns
    -------
    |bool|
    '''
    if visited is None:
        visited = set()
    if (self.id, other.id) in visited:
        print ("Already visited ", (self.id, other.id))
        return True
    if self._flat_equality(other) and (not substituents or self._match_substituents(other)):
        taken_b = set()
        b_children = list(other.children())
        a_children = list(self.children())
        for a_pos, a_child in a_children:
            matched = False
            for b_pos, b_child in b_children:
                if (b_pos, b_child.id) in taken_b:
                    continue
                if a_child.topological_equality(b_child,
                                                substituents=substituents,
                                                visited=visited):
                    matched = True
                    taken_b.add((b_pos, b_child.id))
                    break
            if not matched and len(a_children) > 0:
                print ("not matched and len(a_children) > 0", self.id, other.id)
                return False
        if len(taken_b) != len(b_children):
            print("len(taken_b) != len(b_children)", self.id, other.id)
            return False
        return True
    print("Not Flat-Equal", self.id, other.id)
    return False


def _match_substituents(self, other):
    '''
    Helper method for matching substituents in an order-independent
    fashion. Used by :meth:`topological_equality`
    '''
    taken_b = set()
    b_num_unknown = len(other.substituent_links[-1])
    unknown_cntr = 0
    b_substituents = list(other.substituents())
    cntr = 0
    for a_pos, a_substituent in self.substituents():
        matched = False
        cntr += 1
        for b_pos, b_substituent in b_substituents:
            if b_pos in taken_b:
                if b_pos != -1:
                    continue
                else:
                    if unknown_cntr < b_num_unknown:
                        unknown_cntr += 1
                    else:
                        continue
            if b_substituent == a_substituent:
                matched = True
                taken_b.add(b_pos)
                break
        if not matched and cntr > 0:
            return False
    if len(taken_b) + unknown_cntr != len(b_substituents):
        return False
    return True


def reparse_database(database):
    from glypy.io import glycoct
    for record in database:
        ct_str = str(record.structure)
        structure = glycoct.loads(ct_str)
        assert structure.mass() == record.structure.mass()
        record.structure = structure
        record.update()
