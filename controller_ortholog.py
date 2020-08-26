from __main__ import app, annotation_orth
from flask import request, jsonify

@app.route("/annotation/ortholog/<seed_id>", methods=["GET"])
def get_ortholog_from_seed_rxn_any(seed_id):
    all_orthologs, matched_orthologs, genome_match = annotation_orth.get_orthologs_from_seed_rxn_id3(seed_id)
    result = annotation_orth.process_data2(all_orthologs, genome_match)

    return jsonify(result)

@app.route("/annotation/ortholog/<seed_id>/<s0>/<s1>/<s2>", methods=["GET"])
def get_ortholog_from_seed_rxn(seed_id, s0, s1, s2):
    all_orthologs, matched_orthologs, genome_match = annotation_orth.get_orthologs_from_seed_rxn_id3(seed_id)
    result = annotation_orth.process_data2(all_orthologs, genome_match)

    return jsonify(result)

@app.route("/annotation/ortholog/<genome_id>/<gene_id>", methods=["GET"])
def get_ortholog_from_gene_genome_id(genome_id, gene_id):
    res = annotation_orth.get_ortholog(genome_id, gene_id)
    column_data = {}
    if not res == None:
        ortholog_id = res['id']
        column_data = annotation_orth.process_data2({
            (genome_id, gene_id, ortholog_id)
        }, {})

    return jsonify(column_data)