import pkg_resources
import json

from copy import deepcopy

from pygly2.utils import StringIO
from pygly2.io import glycoct


class MonosaccharideIndex(dict):
    def __init__(self, stream=None):
        if stream is None:
            stream = pkg_resources.resource_stream(__name__, "data/monosaccharides.json")
        self.update(json.load(stream))
        for k, v in self.items():
            self[k] = iter(glycoct.loads(v)).next().root

    def __getitem__(self, key):
        x = dict.__getitem__(self, key)
        ret = deepcopy(x)
        assert(id(x) != id(ret))
        return ret

    def __repr__(self):
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

    def __getitem__(self, key):
        x = dict.__getitem__(self, key)
        ret = deepcopy(x)
        assert(id(x) != id(ret))
        return ret

    def __repr__(self):
        rep = StringIO()
        fmt = "{0}\n{1}\n{2}\n\n"
        for k, v in self.items():
            rep.write(fmt.format(k, v.to_glycoct(), v.mass()))
        return rep.getvalue()