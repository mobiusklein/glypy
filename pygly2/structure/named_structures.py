import pkg_resources
import json
import uuid
import re

from copy import deepcopy

from pygly2.utils import StringIO, identity
from pygly2.utils.lazy import ProxyObject
from pygly2.io import glycoct


class StructureIndex(dict):
    def __init__(self, stream, key_transform=identity, value_transform=identity):
        self.update(json.load(stream))
        for k, v in self.items():
            self[key_transform(k)] = value_transform(glycoct.loads(v).next())
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
        return list(self.__dict__) + list(k for k in self
                                          if (not k[0].isdigit()) and (not re.search(r"[:,\-\s\(\)\[\]]", k)))

    def __repr__(self):  # pragma: no cover
        rep = StringIO()
        fmt = "{0}\n{1}\n{2}\n\n"
        for k, v in self.items():
            rep.write(fmt.format(k, v.to_glycoct(), v.mass(True)))
        return rep.getvalue()


class MonosaccharideIndex(StructureIndex):
    def __init__(self, stream=None, key_transform=identity, value_transform=lambda x: x.root):
        if stream is None:
            stream = pkg_resources.resource_stream(__name__, "data/monosaccharides.json")
        super(MonosaccharideIndex, self).__init__(stream, key_transform, value_transform)

monosaccharides = (MonosaccharideIndex)()


class GlycanIndex(StructureIndex):
    def __init__(self, stream=None, key_transform=identity, value_transform=identity):
        if stream is None:
            stream = pkg_resources.resource_stream(__name__, "data/glycans.json")
        super(GlycanIndex, self).__init__(stream, key_transform, value_transform)

glycans = ProxyObject(GlycanIndex)


class MotifIndex(StructureIndex):
    def __init__(self, stream=None, key_transform=identity, value_transform=identity):
        if stream is None:
            stream = pkg_resources.resource_stream(__name__, "data/motifs.json")
        data = json.load(stream)
        for motif in data:
            name = motif['name']
            motif_class = motif['class']
            motif_category = motif['category']
            motif_structure = glycoct.loads(motif['glycoct']).next()
            motif_structure.motif_name = name
            motif_structure.motif_class = motif_class
            motif_structure.motif_category = motif_category
            self[name] = motif_structure

motifs = ProxyObject(MotifIndex)
