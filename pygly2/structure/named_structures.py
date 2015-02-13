import pkg_resources
import json
import uuid
import re

from copy import deepcopy

from pygly2.utils import StringIO, identity
from pygly2.io import glycoct


class MonosaccharideIndex(dict):
    def __init__(self, stream=None, key_transform=identity):
        if stream is None:
            stream = pkg_resources.resource_stream(__name__, "data/monosaccharides.json")
        self.update(json.load(stream))
        for k, v in self.items():
            self[key_transform(k)] = iter(glycoct.loads(v)).next().root
        self.key_transform = key_transform

    def __getitem__(self, key):
        x = dict.__getitem__(self, key)
        ret = deepcopy(x)
        ret.id = uuid.uuid4().int
        return ret

    def __getattr__(self, name):
        try:
            res = object.__getattr__(self, name)
            return res
        except AttributeError:
            return self[name]

    def __dir__(self):
        return list(self.__dict__) + list(k for k in self if (not k[0].isdigit()) and (not re.search(r"[:,-\s\(\)\[\]]", k)))

    def __repr__(self):  # pragma: no cover
        rep = StringIO()
        fmt = "{0}\n{1}\n{2}\n\n"
        for k, v in self.items():
            rep.write(fmt.format(k, v.to_glycoct(), v.mass(True)))
        return rep.getvalue()

monosaccharides = MonosaccharideIndex()


class GlycanIndex(dict):
    def __init__(self, stream=None):
        if stream is None:
            stream = pkg_resources.resource_stream(__name__, "data/glycans.json")
        self.update(json.load(stream))
        for k, v in self.items():
            self[k] = iter(glycoct.loads(v)).next()

    def __getattr__(self, name):
        try:
            res = object.__getattr__(self, name)
            return res
        except AttributeError:
            return self[name]

    def __getitem__(self, key):
        x = dict.__getitem__(self, key)
        ret = deepcopy(x)
        return ret

    def __repr__(self):  # pragma: no cover
        rep = StringIO()
        fmt = "{0}\n{1}\n{2}\n\n"
        for k, v in self.items():
            rep.write(fmt.format(k, v.to_glycoct(), v.mass()))
        return rep.getvalue()

glycans = GlycanIndex()
