import logging
import json
from __main__ import app, escher_manager, modelseed_local, bios
from flask import request, jsonify, abort
from biosapi.report.report_model_mapping import report_species_mapping
from biosapi.core.model.merge_model_builder import MergeModelBuilder
from utils import clear_nan, fat_load_compound


@app.route("/bios/sbml", methods=["GET"])
def bios_list_models():
    res = None
    try:
        res = bios.get_models()
    except Exception as e:
        abort(400)
    return jsonify(res)


# @app.route("/model", methods=["GET"])
# def list_models():
#
#    return jsonify([])
@app.route("/model/<id>/cmp", methods=["GET"])
def list_model_compartments(id):
    resp = bios.get_model_compartments(id)

    return jsonify(resp)


@app.route("/model/<id>/spi", methods=["GET"])
def list_model_species(id):
    resp = bios.get_model_species(id)

    return jsonify(resp)


@app.route("/model/<id>/rxn", methods=["GET"])
def list_model_reactions(id):
    resp = bios.get_model_reactions(id)

    return jsonify(resp)


@app.route("/bios/merge-model", methods=["POST"])
def bios_get_merge_model():
    data = request.get_json()
    if 'model_ids' not in data:
        abort(400)

    merge_model_builder = MergeModelBuilder(bios)
    for model_id in data['model_ids']:
        merge_model_builder.with_model(model_id)
    merge_model = merge_model_builder.build()
    return jsonify(merge_model)


@app.route("/bios/sbml/<model_id>/spi/<spi_id>", methods=["GET"])
def bios_get_model_species(model_id, spi_id):
    res = None
    try:
        res = bios.get_model_specie(model_id, spi_id)
    except Exception as e:
        abort(400)
    return jsonify(res)


@app.route("/bios/sbml/<model_id>/spi", methods=["GET"])
def bios_get_all_model_species(model_id):
    res = None
    try:
        res = bios.get_model_species(model_id)
    except Exception as e:
        abort(400)
    return jsonify(res)


@app.route("/bios/sbml/<model_id>/rxn/<rxn_id>", methods=["GET"])
def bios_get_model_reaction(model_id, rxn_id):
    res = None
    try:
        res = bios.get_model_reaction(model_id, rxn_id)
    except Exception as e:
        abort(400)
    return jsonify(res)


@app.route("/bios/sbml/<model_id>/rxn", methods=["GET"])
def bios_get_all_model_reactions(model_id):
    res = None
    try:
        res = bios.get_model_reactions(model_id)
    except Exception as e:
        abort(400)
    return jsonify(res)


@app.route("/bios/biodb/<database_id>/cpd/<cpd_id>", methods=["GET"])
def bios_get_database_compound(database_id, cpd_id):
    res = None
    if database_id == 'seed.compound':
        res = fat_load_compound(cpd_id, modelseed_local)
        res = res.data
    else:
        res = bios.get_metabolite(cpd_id, database_id)

    return jsonify(res)


@app.route("/bios/biodb/<database_id>/rxn/<rxn_id>", methods=["GET"])
def bios_get_database_reaction(database_id, rxn_id):
    res = None
    if database_id == 'seed.reaction':
        res = modelseed_local.get_seed_reaction(rxn_id)
        if res is not None:
            clear_nan(res.data)
            if rxn_id in modelseed_local.reaction_aliases:
                res.data['aliases'] = {}
                for database in modelseed_local.reaction_aliases[rxn_id]:
                    res.data['aliases'][database] = list(modelseed_local.reaction_aliases[rxn_id][database])
            res = res.data
    else:
        res = bios.get_reactions([rxn_id], database_id)
        if len(res) > 0:
            res = res[0].json_data

    return jsonify(res)


@app.route("/bios/sbml-spi-mapping-report", methods=["POST"])
def bios_post_model_species_mapping_report():
    data = request.get_json()
    if 'score' not in data or 'database' not in data or 'model_ids' not in data:
        abort(400)
    score = data['score']
    report = report_species_mapping(data['database'], score, data['model_ids'], bios)
    return jsonify(report)


@app.route("/bios/sbml-spi-mapping", methods=["POST"])
def bios_post_model_species_mapping():
    data = request.get_json()
    if 'user' not in data or 'score' not in data or 'mapping-species' not in data:
        abort(400)
    result = {}
    mapping = data['mapping-species']
    user = data['user']
    score = data['score']

    database_mapping = {
        'seed.compound': 'ModelSeed'
    }

    model_ids = set()
    for database_tuple in mapping:
        for model_tuple in mapping[database_tuple]:
            spi_id, model_id = model_tuple.split('@')
            model_ids.add(model_id)
    model_spi_data = {}
    for model_id in model_ids:
        print('loading data from: ', model_id)
        model_spi_data[model_id] = {}
        species = bios.get_model_species(model_id)
        for o in species:
            if 'id' in o and o['id'] not in model_spi_data[model_id]:
                model_spi_data[model_id][o['id']] = o
            else:
                print('exclude', model_id, o['entry'])

    for database_tuple in mapping:
        result[database_tuple] = {}
        cpd_id, database_id = database_tuple.split('@')
        for model_tuple in mapping[database_tuple]:
            spi_id, model_id = model_tuple.split('@')
            if spi_id in model_spi_data[model_id]:
                spi = model_spi_data[model_id][spi_id]
                spi_uid = spi['bios_id']
                bios_database_id = database_id
                if bios_database_id in database_mapping:
                    bios_database_id = database_mapping[bios_database_id]
                resp = bios.set_annotation_model_species(model_id, spi_uid,
                                                         bios_database_id, cpd_id,
                                                         user, score)
                if not resp.status_code == 200:
                    result[database_tuple][model_tuple] = False
                else:
                    result[database_tuple][model_tuple] = True
            else:
                result[database_tuple][model_tuple] = False

    return jsonify(result)

