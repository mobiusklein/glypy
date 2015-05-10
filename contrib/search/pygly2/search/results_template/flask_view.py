import sys
import traceback
import json
from flask import Flask, request, session, g, redirect, url_for, \
     abort, render_template, flash, Markup, make_response, jsonify
import argparse

from pygly2.search import matching
from pygly2.search import report

app = Flask(__name__)


DATABASE = None
DEBUG = True


def connect_db():
    return matching.ResultsDatabase(DATABASE)


@app.route("/matches")
def get_id_list():
    return jsonify(matches=(map(dict, g.db.execute("select glycan_id, match_count from {table_name};"))))


@app.route("/<int:structure_id>")
def match_entry(structure_id):
    try:
        results = [g.db[structure_id]]
        metadata = g.db.get_metadata()
        print len(results)
        template = report.render(results, **metadata)
    except Exception, e:
        template = '<br>'.join((str(e), traceback.format_exc(sys.exc_info()[2])))
        print e
    print len(results)
    return (template)


@app.before_request
def before_request():
    g.db = connect_db()


@app.teardown_request
def teardown_request(exception):
    db = getattr(g, 'db', None)
    if db is not None:
        db.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser('view-results')
    parser.add_argument("results_database")
    args = parser.parse_args()
    DATABASE = args.results_database
    app.debug = DEBUG
    app.run()
