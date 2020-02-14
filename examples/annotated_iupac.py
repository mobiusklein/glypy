'''An example of extending the IUPAC reader/writer to support annotated monosaccharides.
'''
import re
import json

import glypy

from glypy.structure.monosaccharide import AnnotatedMonosaccharide
from glypy.io import iupac


class AnnotatedMonosaccharideDeserializer(iupac.MonosaccharideDeserializer):
    '''Define a IUPAC-ish parser which can read an extra group of key=value
    pairs separated by semicolons enclosed by curly-braces. The values are
    JSON-valued, but do not support objects, for simplicity of parsing.
    '''
    _pattern = r'''(?:(?P<anomer>[abo?]|alpha|beta|\u03B1|\u03B2)-)?
                   (?P<configuration>[LD?])-
                   (?P<modification>[a-z0-9_\-,]*?)
                   (?P<base_type>(?:[A-Z][a-z]{2}?|(?:[a-z]{3}[A-Z][a-z]{2})))
                   (?P<ring_type>[xpfo?])?
                   (?P<substituent>[^-{]*?)
                   (?:\{(?P<annotations>[^}]*?)\})?
                   (?P<linkage>-\([0-9?/]+->?[0-9?/]+\)-?)?
                   $'''
    try:
        # convert to unicode for Py2
        _pattern = _pattern.decode("raw_unicode_escape")
    except AttributeError:
        pass
    pattern = re.compile(_pattern, re.VERBOSE | re.UNICODE)

    def build_residue(self, match_dict):
        residue, linkage = super(
            AnnotatedMonosaccharideDeserializer, self).build_residue(match_dict)
        residue = residue.clone(monosaccharide_type=AnnotatedMonosaccharide)
        annotations = match_dict.get('annotations')
        if annotations:
            for token in annotations.split(";"):
                if not token:
                    continue
                key, value = token.split("=")
                value = json.loads(value)
                residue.annotations[key] = value
        return residue, linkage


class AnnotatedMonosaccharideSerializer(iupac.MonosaccharideSerializer):
    '''Converts (Annotated)Monosaccharides, preserving and propagating any
    annotations.
    '''

    def monosaccharide_to_iupac(self, residue):
        residue_string = super(
            AnnotatedMonosaccharideSerializer, self).monosaccharide_to_iupac(residue)
        try:
            annotations = residue.annotations
            if annotations:
                tokens = [
                    "%s=%s" % (key, json.dumps(value)) for key, value in annotations.items()
                ]
                residue_string += ("{%s}" % ";".join(tokens))
        except AttributeError:
            pass
        return residue_string


# Create full reader and writer structures
ags = iupac.GlycanSerializer(AnnotatedMonosaccharideSerializer())
agd = iupac.GlycanDeserializer(AnnotatedMonosaccharideDeserializer())


if __name__ == "__main__":

    string = "a-D-Manp{arbitrary=true;json_valued=[1,2,3,4];annotations=\"here\"}"
    print("Parsing %r" % string)
    residue, linkage = AnnotatedMonosaccharideDeserializer()(string)
    print(residue, residue.annotations)
    print("\n")
    # Transcode a glycan made up of normal Monosaccharide instances into one of
    # AnnotatedMonosaccharide instances.
    print("Annotating N-Linked Core")
    ex = agd(ags(glypy.glycans['N-Linked Core']))
    ex[4].annotations['terminal'] = True
    ex[3].annotations['terminal'] = True
    print(ags(ex))
