import logging
from modelseed_escher.model.escher_from_modelseed import EscherModelSEEDFactory
from modelseed_escher.core import EscherMap
from modelseed_escher.map.model_merge import RefitMap
from __main__ import app, escher_manager, modelseed_local
from flask import request, jsonify, abort
from escher_grid_merge import generate_integration_report, merge_with_layer
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


@app.route("/escher/dataset/<dataset>/map/<map_id>", methods=["POST"])
def flask_post_refit_escher_map(dataset, map_id):
    params = request.get_json()
    if 'cmp_config' not in params:
        abort(400)
    cmp_config_str = params['cmp_config']
    cmp_config = dict(map(lambda x: x.split(':'), cmp_config_str.split(';')))
    em = None
    if 'escher_map' in params:
        em = EscherMap(params['escher_map'])
    elif '.' in map_id:
        a, b = map_id.split('.', 1)
        em = escher_manager.get_map(dataset, a, b)

    if em:
        refit = RefitMap(modelseed_local)
        em = refit.refit(em, cmp_config)
        return jsonify(em.escher_map)

    abort(400)


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


@app.route("/escher/cluster", methods=["POST"])
def escher_cluster_map():
    cluster_data = request.get_json()
    # with open('/Users/fliu/workspace/jupyter/data/www/annotation/data/latest_cluster.json', 'w') as f:
    #    f.write(json.dumps(cluster_data))
    #print(cluster_data.keys())
    #print(type(cluster_data['escher_map']))
    with open('/Users/fliu/workspace/jupyter/data/escher/temp_maps/server_debug_cluster.json', 'w') as fh:
        fh.write(json.dumps(cluster_data))
    report = generate_integration_report(cluster_data, escher_manager)

    return jsonify(report)


@app.route("/escher/merge", methods=["POST"])
def escher_merge_map():
    cluster_data = request.get_json()
    em_merge = merge_with_layer(cluster_data, escher_manager, modelseed_local)

    return jsonify(em_merge.escher_map)


@app.route("/escher/biochem", methods=["POST"])
def flask_escher_build_biochem():
    params = request.get_json()
    f = EscherModelSEEDFactory(modelseed_local)
    f.cmp_config = params['cmp_config']
    model = f.build()
    return jsonify(model)
