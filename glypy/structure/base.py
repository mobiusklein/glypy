
class MoleculeBase(object):
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


class SaccharideBase(MoleculeBase):
    node_type = object()


class SaccharideCollection(SaccharideBase):

    def _derivatized(self, substituent, id_base):
        pass

    def _strip_derivatization(self):
        pass


class SubstituentBase(MoleculeBase):
    node_type = object()
    pass


class ModificationBase(MoleculeBase):
    node_type = object()
    pass
