import cobra
import json
from __main__ import app
from flask import request, jsonify, abort


MODEL_CACHE = {}

@app.route("/optimization/model", methods=["POST"])
def flask_post_get_optmization_model():
    data = request.get_json()
    if 'user' not in data:
        abort(400)

    model_list = set()
    if data['user'] in MODEL_CACHE:
        model_list = set(MODEL_CACHE[data['user']])

    return jsonify(model_list)


@app.route("/optimization/model/<model_id>", methods=["PUT"])
def flask_put_optmization_model(model_id):
    data = request.get_json()
    if 'user' not in data or 'model' not in data:
        abort(400)

    if data['user'] not in MODEL_CACHE:
        MODEL_CACHE[data['user']] = {}

    MODEL_CACHE[data['user']][model_id] = cobra.io.from_json(json.dumps(data['model']))
    return jsonify(True)


@app.route("/optimization/model/<model_id>/media", methods=["POST"])
def flask_post_optmization_model_media(model_id):
    data = request.get_json()
    if 'user' not in data or model_id not in MODEL_CACHE[data['user']]:
        abort(400)
    pass


@app.route("/optimization/model/<model_id>/fba", methods=["POST"])
def flask_post_optmization_model_fba(model_id):
    data = request.get_json()
    if 'user' not in data or model_id not in MODEL_CACHE[data['user']]:
        abort(400)

    model = MODEL_CACHE[data['user']][model_id]
    solution = model.optimize()
    solution.fluxes.to_dict()

    return jsonify(solution)
