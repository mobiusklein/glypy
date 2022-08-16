import pkg_resources
import json

try:
    from urllib import urlopen
except ImportError:
    from urllib.request import urlopen

from collections import defaultdict

from six import string_types as basestring

import glypy.io
from glypy.utils import StringIO
from glypy.utils.lazy import ProxyObject


class EnzymeCommissionNumber(object):
    """Represents an Enzyme Commission (E.C.) number specifying
    a unique identifier for a particular enzyme

    Attributes
    ----------
    activity : int
        The second-most precise number, describing the method of action
        from the particular :attr:`group`
    category : int
        The broadest number
    group : int
        The second broadest number
    identity : int
        The most precise number, describes which enzyme with the particular
        activity designated by :attr:`activity`
    """

    def __init__(self, category, group, activity, identity):
        self.category = int(category)
        self.group = int(group)
        self.activity = int(activity)
        self.identity = int(identity)
        self._str = '.'.join(map(str, self))

    def __eq__(self, other):
        return str(self) == str(other)

    def __hash__(self):
        return hash(self._str)

    def __ne__(self, other):
        return not (self == other)

    def __getitem__(self, i):
        if i == 0:
            return self.category
        elif i == 1:
            return self.group
        elif i == 2:
            return self.activity
        elif i == 3:
            return self.identity
        elif isinstance(i, slice) or i < 0:
            return tuple(self)[i]
        else:
            raise IndexError(i)

    def __iter__(self):
        yield self.category
        yield self.group
        yield self.activity
        yield self.identity

    def __str__(self):
        return self._str

    def __repr__(self):
        return ("{self.__class__.__name__}({self.category}, {self.group},"
                " {self.activity}, {self.identity})").format(self=self)

    @classmethod
    def parse(cls, string):
        string = str(string)
        parts = string.split(".")
        if len(parts) != 4:
            raise ValueError(
                "EC Numbers must have 4 parts (found %d): %r" % (len(parts, string)))
        try:
            parts = tuple(map(int, parts))
        except ValueError:
            raise ValueError("EC Numbers must be integers: %r" % (parts,))
        return cls(*parts)


class EnzymeInformation(object):
    """A collection of descriptive and identifying information about a single enzyme

    Attributes
    ----------
    alternative_names : :class:`list`
        Alternative names for the described enzyme
    ec_number : :class:`EnzymeCommisionNumber`
        The unique EC number for the described enzyme
    extra_information : :class:`dict`
        Arbitrary key-value information about described enzyme
    name : :class:`str`
        A human-readable name for the described enzyme
    """

    def __init__(self, name, ec_number=None, alternative_names=None, **kwargs):
        if isinstance(ec_number, basestring):
            ec_number = EnzymeCommissionNumber.parse(ec_number)
        elif isinstance(ec_number, (list, tuple)):
            ec_number = EnzymeCommissionNumber(*ec_number)

        self.name = name
        self.ec_number = ec_number
        self.alternative_names = alternative_names or ()
        self.extra_information = kwargs

    def __repr__(self):
        return ("{self.__class__.__name__}({self.name!r}, {self.ec_number!r}, "
                "{self.alternative_names!r}, **{self.extra_information!r})").format(self=self)

    def __eq__(self, other):
        if isinstance(other, basestring):
            return other in (self.name, self.ec_number) or other in self.alternative_names
        else:
            try:
                return other.ec_number == self.ec_number
            except AttributeError:
                return False

    def _to_dict(self):
        store = dict()
        store["name"] = self.name
        store['ec_number'] = str(self.ec_number)
        store['alternative_names'] = tuple(self.alternative_names)
        store['extra_information'] = dict(self.extra_information)
        return store

    def __hash__(self):
        return hash(self.name)


class EnzymeDatabase(object):
    _expasy_url = "ftp://ftp.expasy.org/databases/enzyme/enzyme.dat"

    def __init__(self, fp=None, format='json'):
        self.layered_store = defaultdict(
            lambda: defaultdict(lambda: defaultdict(dict)))
        self.direct_store = dict()

        if fp is not None:
            if format == 'expasy':
                self._parse_expasy(fp)
            elif format == 'json':
                self._parse_json(fp)
            else:
                raise ValueError("Unrecognized format: %r" % (format,))

    def add(self, enzyme_info):
        parts = enzyme_info.ec_number
        store = self.layered_store
        for k in parts[:-1]:
            store = store[k]
        store[parts[-1]] = enzyme_info
        self.direct_store[str(parts)] = enzyme_info

    def __getitem__(self, key):
        if isinstance(key[0], int):
            key = EnzymeCommissionNumber(*key)
        else:
            key = EnzymeCommissionNumber.parse(str(key))
        store = self.layered_store
        for k in key[:-1]:
            store = store[k]
        try:
            enzyme_info = store[key[-1]]
        except KeyError:
            enzyme_info = self.direct_store[str(key)]
        return enzyme_info

    def _parse_json(self, fp):
        enzymes = json.load(fp)

        for enz in enzymes:
            self.add(EnzymeInformation(**enz))

    def _dump(self, fp):
        json.dump([e._to_dict() for e in self.direct_store.values()],
                  fp, indent=2, sort_keys=True)

    @classmethod
    def _build(cls):
        fp = urlopen(cls._expasy_url)
        return cls(fp, "expasy")

    def _parse_expasy(self, fp):
        enzymes = []

        def new_store():
            return defaultdict(list)

        def postprocess(enzyme_info):
            enzyme_info['name'] = ' '.join(enzyme_info['name']).strip(".")
            enzyme_info['alternative_names'] = [
                n for n in ''.join(enzyme_info['alternative_names']).split(".") if n]
            up = enzyme_info['uniprot']
            up = [p.strip() for ent in up for p in ent.split(";") if p]
            enzyme_info['uniprot'] = up
            pro = enzyme_info['prosite']
            pro = [p.strip() for ent in pro for p in ent.split(";")
                   if p and p != 'PROSITE']
            enzyme_info['prosite'] = pro
            enzyme_info['catalytic_activity'] = ' '.join(
                enzyme_info['catalytic_activity'])
            comments = enzyme_info['comments']
            comments = ' '.join(
                map(str.strip, map(str.rstrip, comments))).split("-!-")
            enzyme_info['comments'] = ''.join(
                [c for c in comments if c.strip()])
            enzyme_info['ec_number'] = EnzymeCommissionNumber.parse(enzyme_info[
                                                                    'id'])
            return enzyme_info

        current_enzyme = new_store()
        for line in fp:
            line = line.strip()
            sigil = line[:2]
            line = line[5:]
            if sigil == '//':
                if current_enzyme['id']:
                    try:
                        enzymes.append(postprocess(current_enzyme))
                    except ValueError as e:
                        if "n" in current_enzyme['id']:
                            pass
                        else:
                            print(e, current_enzyme['name'])
                current_enzyme = new_store()
            elif sigil == 'ID':
                current_enzyme['id'] = line
            elif sigil == 'DE':
                current_enzyme['name'].append(line.strip(" "))
            elif sigil == 'AN':
                current_enzyme['alternative_names'].append(line.strip(" "))
            elif sigil == 'CA':
                current_enzyme['catalytic_activity'].append(line.strip(" "))
            elif sigil == 'CF':
                current_enzyme['cofactors'].append(line.strip(" "))
            elif sigil == 'CC':
                current_enzyme['comments'].append(line)
            elif sigil == 'PR':
                current_enzyme['prosite'].append(line)
            elif sigil == 'DR':
                current_enzyme['uniprot'].append(line)
            else:
                print(current_enzyme.get('id'), sigil, line)

        for enz in enzymes:
            self.add(EnzymeInformation(**enz))

    @classmethod
    def _from_static(cls):
        data_buffer = pkg_resources.resource_string(
            glypy.io.__name__, "data/enzyme.json")
        if isinstance(data_buffer, bytes):
            data_buffer = data_buffer.decode("utf-8")
        return cls(StringIO(data_buffer), format='json')


expasy_enzyme_db = ProxyObject(EnzymeDatabase._from_static)
