import logging
import copy
from __main__ import app, clear_nan, annotation_api_atlas, annotation_api
from flask import request, jsonify, abort
from cobrakbase import KBaseAPI
from cobrakbase.core.kbasefba import NewModelTemplate

logger = logging.getLogger(__name__)


@app.route("/template/<template_id>/reaction/<rxn_id>", methods=["POST"])
def set_annotation_to_template(template_id, rxn_id):
    params = request.get_json()

    res = annotation_api_atlas.add_function_to_template_rxn(
        int(params['function_id']),
        rxn_id,
        params['user_id'],
        template_id,
        params['logic'])

    if res:
        return jsonify(True)
    else:
        abort(400)


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
    if 'user_id' not in body or 'comment' not in body:
        return 'terrible request this should return error code needs user_id and comment_id', 400
    logger.debug("/template/%s/reaction/%s/comment", template_id, rxn_id)
    res = annotation_api_atlas.add_template_reaction_comment(rxn_id, template_id, body['user_id'], body['comment'])
    return jsonify(res)


@app.route("/template/<template_id>/reaction/<rxn_id>/attributes", methods=["POST"])
def add_template_reaction_attribute(template_id, rxn_id):
    body = request.get_json()
    if 'attribute' not in body or 'value' not in body:
        return 'terrible request this should return error code needs attribute and value', 400
    logger.debug("/template/%s/reaction/%s/attribute", template_id, rxn_id)
    res = annotation_api_atlas.add_template_reaction_attribute(rxn_id, template_id, body['attribute'], body['value'])
    return jsonify(res)


@app.route("/template/<template_id>/reaction/<rxn_id>/attributes", methods=["GET"])
def get_template_reaction_attributes(template_id, rxn_id):
    logger.debug("/template/%s/reaction/%s/attribute", template_id, rxn_id)
    data = annotation_api_atlas.get_template_reaction_attributes(rxn_id, template_id)
    return jsonify(data)


@app.route("/template/<template_id>/annotation/reaction/<rxn_id>/ko/<ko_id>", methods=["POST"])
def post_template_annotation_reaction_manual_ko(template_id, rxn_id, ko_id):
    body = request.get_json()
    if 'user' not in body and 'value' not in body:
        abort(400)
    o = annotation_api.get_node('KeggOrthology', ko_id)
    if o is not None:
        print(o)
        annotation_api_atlas.set_manual_ko(rxn_id, template_id, ko_id, body['value'], body['user'])
        return jsonify(True)
    abort(404)


@app.route("/template/<template_id>/annotation/reaction/<rxn_id>/function/<function_id>", methods=["POST"])
def post_template_annotation_reaction_manual_function(template_id, rxn_id, function_id):
    body = request.get_json()
    if 'user' not in body and 'value' not in body:
        abort(400)
    o = annotation_api.get_function_by_uid(function_id)
    if o is not None:
        print(o)
        annotation_api_atlas.set_manual_function(rxn_id, template_id, function_id, body['value'], body['user'])
        return jsonify(True)
    abort(404)


@app.route("/template/<template_id>/annotation/reaction/<rxn_id>/list", methods=["GET"])
def get_template_annotation_t_reaction_list(template_id, rxn_id):
    res = {}
    for doc in annotation_api_atlas.collection_templates_reactions.find():
        doc_rxn_id, doc_template_id = doc['_id'].split('@')
        if 'annotation' in doc and \
                'seed__DOT__reaction' in doc['annotation'] and \
                doc['annotation']['seed__DOT__reaction'] == rxn_id and \
                template_id == doc_template_id:
            cmp_config = '0:'
            if 'cmp' in doc:
                cmp_config = doc['cmp']
            res[doc_rxn_id] = cmp_config

    return jsonify(res)
