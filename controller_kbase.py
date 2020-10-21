import logging
import json
from __main__ import app, modelseed_local, annotation_api, annotation_api_atlas
from flask import request, jsonify, abort
import cobra
from cobrakbase import KBaseAPI
from cobrakbase.core.kbasefba.fbamodel_builder import FBAModelBuilder
from export_template_to_kbase import export_template

logger = logging.getLogger(__name__)


@app.route("/kbase/reference_info/<reference>", methods=["GET"])
def flask_get_reference_info(reference):
    token = request.headers.get('Authorization')
    if token is None:
        abort(400)
    api = KBaseAPI(token)
    result = api.get_object_info_from_ref(reference)
    return jsonify(result)


@app.route("/kbase/ws/<workspace>", methods=["GET"])
def flask_list_workspace(workspace):
    token = request.headers.get('Authorization')
    if token is None:
        abort(400)
    api = KBaseAPI(token)
    result = api.list_objects(workspace)
    return jsonify(result)


@app.route("/kbase/ws/<workspace>/<object_id>", methods=["GET"])
def flask_get_object(workspace, object_id):
    token = request.headers.get('Authorization')
    if token is None:
        abort(400)
    api = KBaseAPI(token)

    if workspace == 'none':
        workspace = None
        object_id = object_id.replace('_', '/')

    result = api.get_object(object_id, workspace)
    return jsonify(result)


@app.route("/kbase/ws/<workspace_id>/<object_id>/<object_version>", methods=["GET"])
def flask_get_object_ref(workspace_id, object_id, object_version):
    token = request.headers.get('Authorization')
    if token is None:
        abort(400)
    api = KBaseAPI(token)

    result = api.get_object("{}/{}/{}".format(workspace_id, object_id, object_version), None)
    return jsonify(result)


@app.route("/kbase/cobra/<workspace>/<object_id>", methods=["GET"])
def flask_get_cobra_model(workspace, object_id):
    token = request.headers.get('Authorization')
    if token is None:
        abort(400)
    api = KBaseAPI(token)

    if workspace == 'none':
        workspace = None
        object_id = object_id.replace('_', '/')

    model = api.get_from_ws(object_id, workspace)
    model_json = json.loads(cobra.io.to_json(model))

    return jsonify(model_json)


@app.route("/kbase/convert/fbamodel", methods=["POST"])
def flask_convert_kbase_to_cobra_model():
    data = request.get_json()

    if 'model' not in data:
        abort(400)

    kbase_model = data['model']

    b = FBAModelBuilder.from_kbase_json(kbase_model, y, None).build()
    model = b.build()
    model_json = json.loads(cobra.io.to_json(model))

    return jsonify(model_json)


@app.route("/kbase/ws/<workspace>/<object_id>", methods=["PUT"])
def flask_save_object(workspace, object_id):
    data = request.get_json()
    token = request.headers.get('Authorization')
    if token is None or 'object_type' not in data or 'object_data' not in data:
        abort(400)
    api = KBaseAPI(token)
    result = api.save_object(object_id, workspace, data['object_type'], data['object_data'])

    return jsonify(result)


@app.route("/kbase/export/template", methods=["PUT"])
def flask_export_template():
    data = request.get_json()

    token = request.headers.get('Authorization')
    input_workspace = data.get('input_workspace')
    input_id = data.get('input_id')
    output_workspace = data.get('output_workspace')
    output_id = data.get('output_id')
    annotation_namespace = data.get('annotation_namespace')
    clear_roles = data.get('clear_roles')
    clear_complexes = data.get('clear_complexes')
    clear_reactions = data.get('clear_reactions')
    rxn_ids = data.get('rxn_ids')
    if token is None:
        abort(400)
    api = KBaseAPI(token)

    input_template = api.get_from_ws(input_id, input_workspace)
    output_template = export_template(input_template, modelseed_local, annotation_api,
                                      annotation_api_atlas.database, annotation_namespace,
                                      reaction_list=rxn_ids,
                                      clear_reactions=clear_reactions,
                                      clear_complexes=clear_complexes,
                                      clear_roles=clear_roles)

    result = api.save_object(output_id, output_workspace, 'KBaseFBA.NewModelTemplate', output_template.get_data())

    return jsonify(result)
