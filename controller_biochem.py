from __main__ import app, clear_nan, modelseed_local
from flask import jsonify

@app.route("/biochem/depict/<structure_type>/<output_format>", methods=["POST"])
def post_biochem_depict(structure_type, output_format):
    data = request.get_json()
    
    chemdust = ChemDUST(CHEMDUST_URL)
    svg_depict = chemdust.get_depict(data['structure'], structure_type, output_format)
    
    return svg_depict

@app.route("/biochem/calculator/<inchi>/svg", methods=["GET"])
def translate_inchi_svg(inchi):
    return ""

@app.route("/biochem/cpd/<id>/svg", methods=["GET"])
def get_compound_svg(id):
        
    return ""

@app.route("/biochem/rxn/<id>", methods=["GET"])
def get_seed_reaction(id):
    o = modelseed_local.get_seed_reaction(id)
    if not o == None:
        clear_nan(o.data)
        return jsonify(o.data)
    return ""

@app.route("/biochem/cpd/<id>", methods=["GET"])
def get_seed_compound(id):
    o = modelseed_local.get_seed_compound(id)
    if not o == None:
        clear_nan(o.data)
        return jsonify(o.data)
    return ""