from . import ring_shape
from pygly2 import monosaccharides
from pygly2.composition import composition_transform
from lxml import etree as ET
from collections import defaultdict

derivatize = composition_transform.derivatize


class FragmentRecord(object):
    def __init__(self, name, mass, permethylated_mass, cleave_1, cleave_2, kind):
        self.name = name
        self.mass = float(mass)
        self.permethylated_mass = float(permethylated_mass)
        self.cleave_1 = int(cleave_1)
        self.cleave_2 = int(cleave_2)
        self.kind = kind
        self.composition = {}

    def __repr__(self):
        return "FragmentRecord({name} {mass})".format(**self.__dict__)


def load_all():
    tree = ET.parse(__file__ + '/../Cross-ring-data.xml')
    sugars = {}
    by_name = {}
    for cls in tree.iterfind(".//class"):
        sugar_class = {}
        for case in cls.iterfind(".//sugar"):
            isomers = [tag.attrib['abbr'] for tag in case.findall(".//isomer")]
            frags = defaultdict(dict)
            for frag in case.iterfind(".//fragment"):
                data = [frag.attrib[n] for n in ["name", "mass_mono", "mass_pm_mono", "cleav1", "cleav2", "type"]]
                rec = FragmentRecord(*data)
                rec.composition = dict(map(lambda x: (x.attrib['element'], x.attrib['count']), frag.findall(".//composition")))
                frags[rec.cleave_1, rec.cleave_2][rec.kind] = rec
            sugar_class[case.attrib['subclass']] = frags
            for name in isomers:
                by_name[name] = frags
        sugars[cls.attrib['name']] = sugar_class
    return by_name


default_test_cases = [
    'Glc',
    "GlcNAc",
    "Xyl",
    'GalA',
    'Fuc',
    'IdoA',
    #"KDN"
]


def run_tests(test_cases=None):
    if test_cases is None:
        test_cases = default_test_cases
    store = load_all()
    for case in test_cases:
        test = store[case]
        target = monosaccharides[case]
        print(case)
        for k, v in test.items():
            target_d = {t.kind: t for t in ring_shape.crossring_fragments(target, k[0], k[1])}
            target_d_permethylated = {t.kind: t for t in ring_shape.crossring_fragments(
                derivatize(target.clone(), "methyl"), k[0], k[1])}
            for kind in {"A", "X"}:
                if "%0.4f" % v[kind].mass != "%0.4f" % target_d[kind].mass():
                    print(k, v[kind], target_d[kind])
                    print("------------------------")

                if "%0.4f" % v[kind].permethylated_mass != "%0.4f" % target_d_permethylated[kind].mass():
                    print('p', k, v[kind], target_d_permethylated[kind])
                    print("------------------------")
