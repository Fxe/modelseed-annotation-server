import logging
import copy
from __main__ import app, clear_nan, annotation_api_atlas
from flask import request, jsonify, abort
from cobrakbase import KBaseAPI
from cobrakbase.core.kbasefba import NewModelTemplate

logger = logging.getLogger(__name__)


@app.route("/template/<template_id>/reaction/<rxn_id>/enable", methods=["POST"])
def set_template_reaction_enable(template_id, rxn_id):
    reaction_template_id = '{}@{}'.format(rxn_id, template_id)
    logger.warning("/template/%s/reaction/%s/enable", template_id, rxn_id)
    #data = annotation_api_atlas.collection_templates_reactions.find_one({'_id' : reaction_template_id})
    data = {}
    print(data)
    return jsonify(data)


@app.route("/template/<template_id>/reaction/<rxn_id>/disable", methods=["POST"])
def set_template_reaction_disable(template_id, rxn_id):
    reaction_template_id = '{}@{}'.format(rxn_id, template_id)
    logger.warning("/template/%s/reaction/%s/disable", template_id, rxn_id)
    #data = annotation_api_atlas.collection_templates_reactions.find_one({'_id' : reaction_template_id})
    data = {}
    print(data)
    return jsonify(data)


@app.route("/template/<template_id>/reaction/<rxn_id>/comment", methods=["GET"])
def get_template_reaction_comment(template_id, rxn_id):
    logger.debug("/template/%s/reaction/%s/comment", template_id, rxn_id)
    res = annotation_api_atlas.get_template_reaction_comment(rxn_id, template_id)
    return jsonify(res)


@app.route("/template/<template_id>/reaction/<rxn_id>/comment", methods=["POST"])
def add_template_reaction_comment(template_id, rxn_id):
    body = request.get_json()
    if not 'user_id' in body or not 'comment' in body:
        return 'terrible request this should return error code needs user_id and comment_id', 400
    logger.debug("/template/%s/reaction/%s/comment", template_id, rxn_id)
    res = annotation_api_atlas.add_template_reaction_comment(rxn_id, template_id, body['user_id'], body['comment'])
    return jsonify(res)


@app.route("/template/<template_id>/reaction/<rxn_id>/attributes", methods=["POST"])
def add_template_reaction_attribute(template_id, rxn_id):
    body = request.get_json()
    if not 'attribute' in body or not 'value' in body:
        return 'terrible request this should return error code needs attribute and value', 400
    logger.debug("/template/%s/reaction/%s/attribute", template_id, rxn_id)
    res = annotation_api_atlas.add_template_reaction_attribute(rxn_id, template_id, body['attribute'], body['value'])
    return jsonify(res)


@app.route("/template/<template_id>/reaction/<rxn_id>/attributes", methods=["GET"])
def get_template_reaction_attributes(template_id, rxn_id):
    logger.debug("/template/%s/reaction/%s/attribute", template_id, rxn_id)
    data = annotation_api_atlas.get_template_reaction_attributes(rxn_id, template_id)
    return jsonify(data)


def export_template_to_kbase(workspace, object_id):
    token = request.args.get('token')
    if token is None:
        abort(400)
    api = KBaseAPI(token)
    template_o = api.get_object('GramNegModelTemplateV2', 'NewKBaseModelTemplates')
    template = NewModelTemplate(copy.deepcopy(template_o), 'tftr', 'tcpx')
    return 0
