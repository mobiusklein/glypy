import pkg_resources
import json


def omit_slice(seq, i):
    if i == 0:
        return seq[1:]
    return seq[:i] + seq[i+1:]

class SynonymIndex(dict):
    def __init__(self, data):
        for group in data:
            for i, name in enumerate(group):
                self[name] = omit_slice(group, i)

class MonosaccharideSynonymIndex(SynonymIndex):
    def __init__(self, stream=None):
        if stream is None:
            stream = pkg_resources.resource_stream(__name__, "data/monosaccharide_synonyms.json")
        super(MonosaccharideSynonymIndex, self).__init__(json.load(stream))

monosaccharides = MonosaccharideSynonymIndex()
