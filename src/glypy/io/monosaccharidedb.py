# pragma: no cover

import requests
from lxml import etree

template_url = "http://www.monosaccharidedb.org/display_monosaccharide.action?output=xml&id={0}"


def get(monosaccharide_id):
    r = requests.request(url=template_url.format(monosaccharide_id),
                         method='GET', headers={'Accept-Language': 'en-US,en;q=0.9'})
    r.raise_for_status()
    root = etree.fromstring(r.text)
    return _unpack_element(root)


def _local_name(element):
    return element.tag


schema_info = {
    'lists': {
        'atom',
        'connection',
        'position',
        'alias',
        'core_modification',
    },
    'strings': {
        'stereocode',
    }
}


def typify(val):
    if val in ("true", "false"):
        if val == 'true':
            return True
        elif val == 'false':
            return False
    try:
        return int(val)
    except Exception:
        try:
            return float(val)
        except Exception:
            return str(val)


def _unpack_element(element, **kwargs):
    try:
        name = kwargs.pop('ename')
    except KeyError:
        name = _local_name(element)
    info = {k: typify(v) for k, v in dict(element.attrib).items()}
    # process subelements
    for child in element.iterchildren():
        cname = _local_name(child)
        if cname not in schema_info['lists']:
            info[cname] = _unpack_element(child, ename=cname, **kwargs)
        else:
            info.setdefault(cname, []).append(_unpack_element(child, ename=cname, **kwargs))

    # process element text
    if element.text:
        stext = element.text.strip()
        if stext:
            val = typify(stext) if name not in schema_info['strings'] else stext
            if info:
                info[name] = val
            else:
                return val
    return info
