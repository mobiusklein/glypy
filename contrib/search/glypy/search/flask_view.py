from flask import Flask, request, session, g, redirect, url_for, \
     abort, render_template, flash, Markup, make_response, jsonify, Response
import argparse

from glypy.search import matching
from glypy.search import report
from glypy.search.spectra import spectrum_model


app = Flask(__name__)
report.prepare_environment(app.jinja_env)


DATABASE = None
DEBUG = True
SECRETKEY = 'TG9yZW0gaXBzdW0gZG90dW0'


def connect_results_db():
    return matching.ResultsDatabase(DATABASE)


def connect_spectrum_db():
    return spectrum_model.MSMSSqlDB(DATABASE)


@app.route("/internal/clear_cache")
def clear_cache():
    app.jinja_env.fragment_cache.clear()
    return jsonify(status=True)


@app.route("/internal/update_preferences", methods=["POST"])
def update_preferences():
    session['show_sidebar'] = request.values.get('showSidebar', True)
    session['score_threshold'] = request.values.get("scoreThreshold", 0.)
    session.modified = True
    return jsonify(status=True)


@app.route("/tandem_spectrum/<int:spectrum_id>")
def plot_spectrum(spectrum_id):
    spectra = g.spectrum_db[spectrum_id]
    try:
        if len(spectra) == 1:
            spectra = spectra[0]
    except:
        pass
    return report.spectrum_plot(spectra.tandem_data)


@app.route("/matches")
def get_id_list():
    return jsonify(matches=(map(dict, g.results_db.execute("select glycan_id from {table_name};"))))


@app.route("/mass_range/<lower>,<upper>")
def mass_range(lower, upper):
    print "Begin ", lower, upper
    results = g.results_db.from_sql(g.results_db.execute(
        "SELECT * FROM {table_name} WHERE mass BETWEEN ? AND ? AND tandem_score >= ?",
        map(float, (lower, upper, session.get('score_threshold', 0)))))
    print results
    print "Fetch Record Complete"
    metadata = g.results_db.get_metadata()
    print "Fetch Metadata Complete"
    template = render_template("results.templ", matches=results, live=True, database=g.results_db, **metadata)
    print "Render Template Complete"
    return Response(template)


@app.route("/composition/<composition>")
def by_composition(composition):
    print "Begin ", composition
    results = g.results_db.from_sql(g.results_db.execute(
        "SELECT * FROM {table_name} WHERE composition = ? AND tandem_score >= ?",
        (composition, session.get('score_threshold', 0))))
    print "Fetch Record Complete"
    metadata = g.results_db.get_metadata()
    print "Fetch Metadata Complete"
    template = render_template("results.templ", matches=results, live=True, database=g.results_db, **metadata)
    print "Render Template Complete"
    return Response(template)


@app.route("/structure/<structure_id>")
def match_entry(structure_id):
    print "Begin ", structure_id
    results = g.results_db[structure_id.split(',')]
    print results
    print "Fetch Record Complete"
    metadata = g.results_db.get_metadata()
    print "Fetch Metadata Complete"
    template = render_template("results.templ", matches=results, live=True, database=g.results_db, **metadata)
    print "Render Template Complete"
    return Response(template)


@app.before_request
def before_request():
    print session
    g.results_db = connect_results_db()
    g.spectrum_db = connect_spectrum_db()


@app.teardown_request
def teardown_request(exception):
    db = getattr(g, 'db', None)
    if db is not None:
        db.close()


def main():
    parser = argparse.ArgumentParser('view-results')
    parser.add_argument("results_database")
    args = parser.parse_args()
    global DATABASE
    DATABASE = args.results_database
    app.debug = DEBUG
    app.secret_key = SECRETKEY
    app.run()

if __name__ == "__main__":
    main()
