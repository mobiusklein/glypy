import pkg_resources
import json


def omit_slice(seq, i):
    if i == 0:
        return seq[1:]
    return seq[:i] + seq[i + 1:]


class SynonymIndex(dict):
    def __init__(self, data):
        for group in data:
            for i, name in enumerate(group):
                self[name] = omit_slice(group, i)


class MonosaccharideSynonymIndex(SynonymIndex):
    def __init__(self, stream=None):
        if stream is None:
            data_buffer = pkg_resources.resource_string(__name__, "data/monosaccharide_synonyms.json")
            if isinstance(data_buffer, bytes):
                data_buffer = data_buffer.decode("utf-8")
        super(MonosaccharideSynonymIndex, self).__init__(json.loads(data_buffer))


#: A mapping of monosaccharide names to their synonyms
monosaccharides = MonosaccharideSynonymIndex()
