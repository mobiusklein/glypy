
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


class SaccharideBase(MoleculeBase):
    node_type = object()
    pass


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
