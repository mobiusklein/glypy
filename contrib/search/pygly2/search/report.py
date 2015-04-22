from collections import defaultdict

from matplotlib import pyplot as plt
from jinja2 import Environment, PackageLoader

from pygly2.composition import composition_transform
from pygly2.utils import StringIO
from pygly2 import plot


def collect_fragments(record):
    matches = defaultdict(list)
    for match in record.matches:
        matches[match.match_key.split(":")[0]].append(match)
    return matches.keys()


def strip_derivatize_glycoct(record):
    s = record.structure.clone()
    composition_transform.strip_derivatization(s)
    return(str(s))


def cfg_plot(record):
    s = record.structure.clone()
    composition_transform.strip_derivatization(s)

    dtree, ax = plot.plot(s, orientation='h', squeeze=1.4, scale=.135)
    fmap = {f.name: f for f in record.fragments}
    for match in record.matches:
        match_key = match.match_key.split(":")[0]
        dtree.draw_cleavage(ax, fmap[match_key], label=True)

    ax.axis('off')
    fig = ax.get_figure()
    fig.tight_layout(pad=0.2)
    img_buffer = StringIO()
    fig.savefig(img_buffer, format="svg")
    plt.close(fig)

    return img_buffer.getvalue()


def scientific_notation(num):
    return "%0.3e" % num


def limit_sigfig(num):
    return "%0.4f" % num


def create_environment():
    loader = PackageLoader("pygly2", "search")
    env = Environment(loader=loader)
    env.filters["collect_fragments"] = collect_fragments
    env.filters["strip_derivatize"] = strip_derivatize_glycoct
    env.filters["scientific_notation"] = scientific_notation
    env.filters["cfg_plot"] = cfg_plot
    env.filters["limit_sigfig"] = limit_sigfig

    template = env.get_template("results.templ")
    return template


def render(matches, settings=None):
    template = create_environment()
    return template.render(matches=matches, settings=settings)
