import logging
import json
import cobra
from modelseed_escher.model.escher_from_modelseed import EscherModelSEEDFactory
from modelseed_escher.core import EscherMap
from modelseed_escher.map.model_merge import RefitMap
from biosapi.io.bios_model_builder import BiosModelToCobraBuilder
from __main__ import app, escher_manager, modelseed_local, bios
from flask import request, jsonify, abort
from escher_grid_merge import generate_integration_report, merge_with_layer

logger = logging.getLogger(__name__)

BIOCHEM_CACHE = {}


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

    if frozenset(params['cmp_config'].items()) in BIOCHEM_CACHE:
        print('flask_escher_build_biochem::from cache')
        return jsonify(BIOCHEM_CACHE[frozenset(params['cmp_config'].items())])

    f = EscherModelSEEDFactory(modelseed_local)
    f.cmp_config = params['cmp_config']
    model = f.build()
    if model:
        BIOCHEM_CACHE[frozenset(params['cmp_config'].items())] = model
    print('flask_escher_build_biochem::from build')
    return jsonify(model)


@app.route("/escher/biochem/<biochem_type>/<biochem_id>", methods=["POST"])
def flask_escher_build_biochem_model(biochem_type, biochem_id):
    params = request.get_json()
    if biochem_type == 'database':
        if biochem_id == 'modelseed':
            f = EscherModelSEEDFactory(modelseed_local)
            f.cmp_config = params['cmp_config']
            model = f.build()
            return jsonify(model)
        else:
            abort(400)
    elif biochem_type == 'sbml':
        b = BiosModelToCobraBuilder.from_api(biochem_id, bios)
        model = b.build()
        model = json.loads(cobra.io.to_json(model))
        return jsonify(model)
    else:
        abort(400)


@app.route("/escher/build/grid", methods=["POST"])
def build_grid_map():
    from escher_factory_api import process_build_data_input, build_escher_factory_api
    build_data = request.get_json()
    print(build_data)
    map_assembly = process_build_data_input(build_data['maps'])
    grid_x = build_data['x']
    grid_y = build_data['y']
    print(map_assembly)

    model_ids = set(map(lambda x: x['sbml_id'], map_assembly))
    escher_factory = build_escher_factory_api(model_ids, escher_manager, bios, 'ModelSeed', 3, 'ModelSeedReaction', 3)
    #escher_factory.model_path = CACHE_BASE_FOLDER
    #escher_factory.cpd_mapping = MODEL_CPD_MAPPING
    #escher_factory.rxn_mapping = MODEL_RXN_MAPPING
    #escher_factory.cmp_mapping = MODEL_CMP_MAPPING
    #escher_factory.bios_cache = bios.model_data

    master = escher_factory.build_grid(map_assembly, (grid_x, grid_y))
    # print(content)
    # build_data = request.json

    # [ "iMM904;c;c;ModelSEED.NAD(P) Biosynthesis", "iMM904;c;c;ModelSEED.Pentose and Glucuronate Interconversions", "iMM904;c;c;ModelSEED.Mannitol Utilization", "iJDZ836;CCO__45__CYTOSOL;c;ModelSEED.NAD(P) Biosynthesis", "iJDZ836;CCO__45__CYTOSOL;c;ModelSEED.Pentose and Glucuronate Interconversions", "iJDZ836;CCO__45__CYTOSOL;c;ModelSEED.Mannitol Utilization" ]

    # base = '/Users/fliu/Library/Caches/escher/1-0-0/5/maps/'
    # m = escher_manager.get_map('ModelSEED', 'ModelSEED', 'Chorismate Synthesis')
    return jsonify(master.escher_map)
