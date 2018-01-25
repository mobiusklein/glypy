import pkg_resources
import hjson
import re

from glypy.utils import StringIO, identity, uid
from glypy.utils.lazy import ProxyObject
from glypy.io import glycoct
from glypy.structure.glycan import NamedGlycan


class StructureIndex(dict):
    def __init__(self, stream, key_transform=identity, value_transform=identity):
        self.update(hjson.load(stream))
        for k, v in self.items():
            self[key_transform(k)] = value_transform(glycoct.loads(v))
        self.key_transform = key_transform

    def __getitem__(self, key):
        x = dict.__getitem__(self, key)
        # ret = deepcopy(x)
        ret = x.clone()
        ret.id = uid()
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

    def __repr__(self):
        return "<{}>".format(self.__class__.__name__)

    def __str__(self):  # pragma: no cover
        rep = StringIO()
        fmt = "{0}\n{1}\n{2}\n\n"
        for k, v in self.items():
            rep.write(fmt.format(k, v.to_glycoct(), v.mass(True)))
        return rep.getvalue()


class MonosaccharideIndex(StructureIndex):
    def __init__(self, stream=None, key_transform=identity, value_transform=lambda x: x.root):
        if stream is None:
            stream = pkg_resources.resource_stream(__name__, "data/monosaccharides.hjson")
        with stream:
            super(MonosaccharideIndex, self).__init__(stream, key_transform, value_transform)


monosaccharides = (MonosaccharideIndex)()


class MonosaccharideResidueIndex(MonosaccharideIndex):
    def __init__(self, stream=None, **kwargs):
        from glypy.structure.glycan_composition import MonosaccharideResidue

        key_transform = identity

        def value_transform(x):
            return MonosaccharideResidue.from_monosaccharide(x.root)

        if stream is None:
            stream = pkg_resources.resource_stream(__name__, "data/monosaccharides.hjson")
        super(MonosaccharideIndex, self).__init__(stream, key_transform, value_transform)


monosaccharide_residues = ProxyObject(MonosaccharideResidueIndex)


class GlycanIndex(StructureIndex):
    def __init__(self, stream=None, key_transform=identity, value_transform=identity):
        if stream is None:
            stream = pkg_resources.resource_stream(__name__, "data/glycans.hjson")
        super(GlycanIndex, self).__init__(stream, key_transform, value_transform)


glycans = ProxyObject(GlycanIndex)


class MotifIndex(StructureIndex):
    def __init__(self, stream=None, key_transform=identity, value_transform=identity):
        if stream is None:
            stream = pkg_resources.resource_stream(__name__, "data/motifs.hjson")
        with stream:
            data = hjson.load(stream)
        motif_classes = set()
        motif_categories = set()
        for motif in data:
            name = motif['name']
            motif_class = motif['class']
            motif_category = motif['category']
            motif_structure = NamedGlycan(name=name, root=glycoct.loads(motif['glycoct']).root, index_method=None)
            motif_structure.motif_name = name
            motif_structure.motif_class = motif_class
            motif_structure.motif_category = motif_category
            motif_structure.is_core_motif = motif["core_motif"]
            self[name] = motif_structure
            motif_classes.add(motif_class)
            motif_categories.add(motif_category)
        self._category_map = {}
        self._class_map = {}
        self.motif_classes = motif_classes
        self.motif_categories = motif_categories

    def motif_category(self, name):
        if name in self._category_map:
            return self._category_map[name]

        mapping = {}
        for k, v in self.items():
            if v.motif_category == name:
                mapping[k] = v
        if len(mapping) > 0:
            self._category_map[name] = mapping
        return mapping

    def motif_class(self, name):
        if name in self._class_map:
            return self._class_map[name]
        mapping = {}
        for k, v in self.items():
            if v.motif_class == name:
                mapping[k] = v
        if len(mapping) > 0:
            self._class_map[name] = mapping
        return mapping


motifs = ProxyObject(MotifIndex)
