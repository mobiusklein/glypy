
class MoleculeBase(object):
    __slots__ = ()
    _order = 0
    node_type = object()

    def order(self, deep=False):
        '''
        Return the "graph theory" order of this molecule

        Returns
        -------
        int
        '''
        return self._order

    def has_undefined_linkages(self):
        return True

    def copy(self, *args, **kwargs):
        return self.clone(*args, **kwargs)


class SaccharideBase(MoleculeBase):
    __slots__ = ()
    node_type = object()


class SaccharideCollection(SaccharideBase):
    __slots__ = ()

    def _derivatized(self, substituent, id_base):
        pass

    def _strip_derivatization(self):
        pass


class SubstituentBase(MoleculeBase):
    __slots__ = ()
    node_type = object()
    pass


class ModificationBase(MoleculeBase):
    __slots__ = ()
    node_type = object()
    pass
