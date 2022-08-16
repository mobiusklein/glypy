
oget = object.__getattribute__
oset = object.__setattr__


whitelist = {
    "_source",
    "_initializer",
    "_prepare"
}


class ProxyObject(object):
    '''A lazy object proxy.

    A simple lazy object for loading resources only at the last possible minute.
    Given 0-argument callable :attr:`_initializer`, construct an object that
    will wait until an attribute or item is called for to invoke
    :attr:`_initializer` and store the resulting object as :attr:`_source`, and
    thereafter serve attributes and items from :attr:`_source`.
    '''
    def __init__(self, initializer=None):
        self._initializer = initializer
        self._source = None

    def _prepare(self):
        self._source = self._initializer()
        self.__getitem__ = self._source.__getitem__
        self.__setitem__ = self._source.__setitem__

    def __getattribute__(self, name):
        if name in whitelist:
            return oget(self, name)
        else:
            if self._source is None:
                self._prepare()
            return getattr(self._source, name)

    def __setattr__(self, name, value):
        if name in whitelist:
            oset(self, name, value)
        else:
            setattr(self._source, name, value)

    def __getitem__(self, key):
        if self._source is None:
            self._prepare()
        return self._source[key]

    def __setitem__(self, key, value):
        if self._source is None:
            self._prepare()
        self._source[key] = value

    def __repr__(self):  # pragma: no cover
        rep = r"<ProxyObject>{}{}"
        sep = ""
        val = ""
        if self._source is not None:
            sep = r'~'
            val = repr(self._source)
        return rep.format(sep, val)

    def __dir__(self):  # pragma: no cover
        if self._source is None:
            self._prepare()
        return dir(self._source)
