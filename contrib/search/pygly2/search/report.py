import re
from itertools import cycle
from collections import defaultdict
import base64
try:
    from lxml import etree as ET
except:
    try:
        from xml.etree import cElementTree as ET
    except:
        from xml.etree import ElementTree as ET

import matplotlib
from matplotlib import pyplot as plt
from matplotlib.colors import cnames

from jinja2 import Environment, PackageLoader, Undefined

from pygly2.composition import composition_transform
from pygly2.utils import StringIO
from pygly2 import plot

from pygly2.search.spectra import spectrum_model

matplotlib.rcParams['svg.fonttype'] = 'path'


def collect_fragments(record):
    matches = defaultdict(list)
    for match in record.matches:
        matches[match.match_key.split(":")[0]].append(match)
    return matches.keys()


def fetch_all_matches(db, attr='mass', order="ASC"):
    order_by = " ORDER BY {} {}".format(attr, order)
    return db.from_sql(db.execute("SELECT * FROM {table_name} WHERE \
        (tandem_scan_count + precursor_scan_count) > 0 " + order_by))


def strip_derivatize_glycoct(record):
    s = record.structure.clone()
    composition_transform.strip_derivatization(s)
    return(str(s))


def cfg_plot(record):
    if "svg_plot" in record.report_data:
        return base64.decodestring(record.report_data["svg_plot"])
    s = record.structure.clone()
    composition_transform.strip_derivatization(s)
    dtree, ax = plot.plot(s, orientation='h', squeeze=1.4, scale=.135)
    fmap = {f.name: f for f in record.fragments}
    for match in record.matches:
        match_key = match.match_key.split(":")[0]
        order = len(match_key.split("-"))
        if order == 1:
            dtree.draw_cleavage(ax=ax, fragment=fmap[match_key], color='red', label=True)
        else:
            for key in match_key.split("-"):
                dtree.draw_cleavage(fragment=fmap[key], ax=ax, color='orange', label=True)

    ax.axis('off')
    fig = ax.get_figure()
    fig.tight_layout(pad=0.2)
    img_buffer = StringIO()
    fig.savefig(img_buffer, format="svg")
    plt.close(fig)

    root, ids = ET.XMLID(img_buffer.getvalue())
    root.set("id", dtree.uuid)
    svg = ET.tostring(root)
    record.report_data["svg_plot"] = base64.encodestring(svg)
    record.update()
    return svg


def simple_plot(record):
    dtree, ax = plot.plot(record, orientation='h', squeeze=1.4, scale=.135)
    ax.axis('off')
    fig = ax.get_figure()
    fig.tight_layout(pad=0.2)
    img_buffer = StringIO()
    fig.savefig(img_buffer, format="svg")
    plt.close(fig)
    root, ids = ET.XMLID(img_buffer.getvalue())
    root.set("id", dtree.uuid)
    svg = ET.tostring(root)
    return svg


greek_symbol_map = {
    "a": "&alpha;",
    "b": "&beta;",
    "c": "&gamma;",
    "d": "&delta;",
    "e": "&epsilon;",
    "f": "&zeta;",
    "g": "&eta;",
    "h": "&iota;",
    "i": "&kappa;",
    "j": "&lambda;",
    "k": "&mu;",
    "l": "&nu;",
    "m": "&xi;",
    "n": "&varsigma;",
    "o": "&pi;",
    "p": "&rho;",
    "q": "&sigma;",
    "r": "&tau;",
    "s": "&upsilon;",
    "t": "&phi;",
    "u": "&psi;",
    "v": "&omega;"
}


def greek_fragment_names(fragment_name):
    buff = []
    for c in fragment_name:
        if c in greek_symbol_map:
            buff.append(greek_symbol_map[c])
        else:
            buff.append(c)
    return ''.join(buff)


def spectrum_plot(precursor_spectrum):
    fig = spectrum_model.plot_observed_spectra(precursor_spectrum, annotation_source=None)
    img_buffer = StringIO()
    fig.savefig(img_buffer, format='svg')
    plt.close(fig)
    root, ids = ET.XMLID(img_buffer.getvalue())
    svg = ET.tostring(root)
    return svg


def scientific_notation(num):
    if num is None or isinstance(num, Undefined):
        return "N/A"
    return "%0.3e" % num


def limit_sigfig(num):
    if num is None or isinstance(num, Undefined):
        return "N/A"
    return "%0.4f" % num


def unique(iterable):
    return set(iterable)


def css_escape(css_string):
    return re.sub(r"[\+\:,\s]", r'-', css_string)


def prepare_environment():
    loader = PackageLoader("pygly2", "search/results_template")
    env = Environment(loader=loader)
    env.filters["collect_fragments"] = collect_fragments
    env.filters["all_matches"] = fetch_all_matches
    env.filters["strip_derivatize"] = strip_derivatize_glycoct
    env.filters["scientific_notation"] = scientific_notation
    env.filters["cfg_plot"] = cfg_plot
    env.filters["simple_cfg_plot"] = simple_plot
    env.filters["min"] = min
    env.filters["max"] = max
    env.filters["unique"] = unique
    env.filters["limit_sigfig"] = limit_sigfig
    env.filters['css_escape'] = css_escape
    env.filters['greek_fragment_name'] = greek_fragment_names
    return env


def prepare_template():
    env = prepare_environment()

    template = env.get_template("results.templ")
    return template


def render(database, experimental_statistics=None, settings=None, live=False, **kwargs):
    template = prepare_template()
    return template.render(
        database=database,
        experimental_statistics=experimental_statistics,
        settings=settings, live=live, **kwargs)
