import logging
from __main__ import app, modelseed_local, annotation_api, annotation_api_atlas
from flask import request, jsonify, abort
from cobrakbase import KBaseAPI
from export_template_to_kbase import export_template

logger = logging.getLogger(__name__)


@app.route("/kbase/reference_info/<reference>", methods=["GET"])
def flask_get_reference_info(reference):
    token = request.args.get('token')
    if token is None:
        abort(400)
    api = KBaseAPI(token)
    result = api.get_object_info_from_ref(reference)
    return jsonify(result)


@app.route("/kbase/ws/<workspace>", methods=["GET"])
def flask_list_workspace(workspace):
    token = request.args.get('token')
    if token is None:
        abort(400)
    api = KBaseAPI(token)
    result = api.list_objects(workspace)
    return jsonify(result)


@app.route("/kbase/ws/<workspace>/<object_id>", methods=["GET"])
def flask_get_object(workspace, object_id):
    token = request.args.get('token')
    if token is None:
        abort(400)
    api = KBaseAPI(token)
    result = api.get_object(object_id, workspace)
    return jsonify(result)


@app.route("/kbase/ws/<workspace>/<object_id>", methods=["PUT"])
def flask_save_object(workspace, object_id):
    data = request.get_json()
    token = request.args.get('token')
    if token is None or 'object_type' not in data or 'object_data' not in data:
        abort(400)
    api = KBaseAPI(token)
    result = api.save_object(object_id, workspace, data['object_type'], data['object_data'])

    return jsonify(result)


@app.route("/kbase/export/template", methods=["PUT"])
def flask_export_template():
    data = request.get_json()
    print(data)
    token = data.get('token')
    input_workspace = data.get('input_workspace')
    input_id = data.get('input_id')
    output_workspace = data.get('output_workspace')
    output_id = data.get('output_id')
    annotation_namespace = data.get('annotation_namespace')
    if token is None:
        abort(400)
    api = KBaseAPI(token)

    input_template = api.get_from_ws(input_id, input_workspace)
    output_template = export_template(input_template, modelseed_local, annotation_api,
                                      annotation_api_atlas.database, annotation_namespace)

    from cobra.core.dictlist import DictList
    temp_object = {}
    for k in output_template.data.keys():
        if k not in ['data', 'info', 'provenance']:
            if type(output_template.data[k]) is DictList:
                temp_object[k] = list(output_template.data[k])
            else:
                temp_object[k] = output_template.data[k]
    
    result = api.save_object(output_id, output_workspace, 'KBaseFBA.NewModelTemplate', temp_object)

    return jsonify(result)
