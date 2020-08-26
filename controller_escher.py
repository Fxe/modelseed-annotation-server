import logging
from __main__ import app, escher_manager
from flask import request, jsonify
import json

logger = logging.getLogger(__name__)


@app.route("/escher/dataset", methods=["GET"])
def list_datasets():
    return jsonify(list(escher_manager.list_datasets()))


@app.route("/escher/dataset/<dataset>/model", methods=["GET"])
def list_model(dataset):
    return jsonify(list(escher_manager.list_datasets(dataset)))


@app.route("/escher/dataset/<dataset>/map", methods=["GET"])
def list_map(dataset):
    return jsonify(list(escher_manager.list_maps(dataset)))


@app.route("/escher/dataset/<dataset>/map/<map_id>", methods=["GET"])
def get_escher_map(dataset, map_id):
    if '.' in map_id:
        a, b = map_id.split('.', 1)
        em = escher_manager.get_map(dataset, a, b)
        return jsonify(em.escher_map)
    return jsonify(list(escher_manager.list_maps(dataset)))


@app.route("/escher/dataset/<dataset>/map/<map_id>/save", methods=["POST"])
def save_escher_map(dataset, map_id):
    escher_path = escher_manager.escher.get_cache_dir()
    data = request.get_json()

    if 'map' in data:
        filename = "{}/maps/{}/{}.{}.json".format(escher_path, dataset, dataset, map_id)
        with open(filename, 'w') as fh:
            fh.write(json.dumps(data['map']))

        return jsonify(data['map'][0])

    return jsonify(False)
