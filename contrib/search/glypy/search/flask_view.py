import sys
import traceback
import json
from flask import Flask, request, session, g, redirect, url_for, \
     abort, render_template, flash, Markup, make_response, jsonify
import argparse

from glypy.search import matching
from glypy.search import report
from glypy.search.spectra import spectrum_model

app = Flask(__name__)


DATABASE = None
DEBUG = True


def connect_results_db():
    return matching.ResultsDatabase(DATABASE)


def connect_spectrum_db():
    return spectrum_model.MSMSSqlDB(DATABASE)


@app.route("/tandem_spectrum/<int:spectrum_id>")
def plot_spectrum(spectrum_id):
    try:
        spectra = g.spectrum_db[spectrum_id]
        try:
            if len(spectra) == 1:
                spectra = spectra[0]
        except:
            pass
        return report.spectrum_plot(spectra.tandem_data)
    except Exception, e:
        template = '<br>'.join((str(e), traceback.format_exc(sys.exc_info()[2])))
        print e
        return template


@app.route("/matches")
def get_id_list():
    return jsonify(matches=(map(dict, g.results_db.execute("select glycan_id from {table_name};"))))


@app.route("/mass_range/<lower>,<upper>")
def mass_range(lower, upper):
    try:
        print "Begin ", lower, upper
        results = g.results_db.from_sql(g.results_db.execute(
            "SELECT * FROM {table_name} WHERE mass BETWEEN ? AND ?", map(float, (lower, upper))))
        print results
        print "Fetch Record Complete"
        metadata = g.results_db.get_metadata()
        print "Fetch Metadata Complete"
        template = report.render(matches=results, live=True, database=g.results_db, **metadata)
        print "Render Template Complete"
    except Exception, e:
        template = '<br>'.join((str(e), traceback.format_exc(sys.exc_info()[2])))
        print e
    return (template)


@app.route("/<structure_id>")
def match_entry(structure_id):
    try:
        print "Begin ", structure_id
        results = g.results_db[structure_id.split(',')]
        print results
        print "Fetch Record Complete"
        metadata = g.results_db.get_metadata()
        print "Fetch Metadata Complete"
        template = report.render(matches=results, live=True, database=g.results_db, **metadata)
        print "Render Template Complete"
    except Exception, e:
        template = '<br>'.join((str(e), traceback.format_exc(sys.exc_info()[2])))
        print e
    return (template)


@app.before_request
def before_request():
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
    app.run()

if __name__ == "__main__":
    main()
