import sys
import json
import logging
import pandas as pd
import time
import yaml

import redis
import pymongo

import cobra
import escher
import cobrakbase

import modelseed_escher

import biosapi
from bios_mock import BIOS_MOCK


from curation_api import CurationApi
from escher_factory_api import process_build_data_input, EscherFactoryApi
from escher_grid_merge import generate_integration_report, merge_with_layer
from annotation_ortholog import build_annotation_ortholog, AnnotationOrtholog
from annotation_api import AnnotationApi
from annotation_api_neo4j import AnnotationApiNeo4j
from annotation_api_redis import AnnotationApiRedisCache
from py2neo import Graph, NodeMatcher, RelationshipMatcher
from flask import Flask, request, jsonify
from flask_restful import Resource, Api
from flask_cors import CORS

#class (Resource):
#    def get(self):
#        return {}
    
logger = logging.getLogger(__name__)

app = Flask(__name__)
CORS(app)
api = Api(app)

HUGE_CACHE = {}

def clear_nan(d):
    for k in d:
        if type(d[k]) == float and pd.isna(d[k]):
            d[k] = ""

@app.route("/reload_temp", methods=["GET"])
def reload():
    bios, MODEL_CMP_MAPPING, MODEL_CPD_MAPPING, MODEL_RXN_MAPPING = load_cache_data(CACHE_BASE_FOLDER)
    return jsonify('yes!')
    
@app.route("/status", methods=["GET"])
def status():
    res = {
        'server' : True,
        'neo4j' : False,
        'mongo_atlas' : False,
        'mongo_atlas_info' : {},
        'neo4j_nodes' : -1,
    }
    try:
        node_count = len(annotation_api.neo4j_graph.nodes)
        res['neo4j'] = True
        res['neo4j_nodes'] = node_count
    except:
        pass
    try:
        cloud_status = annotation_api_atlas.server_info()
        res['mongo_atlas'] = True
        res['mongo_atlas_info'] = cloud_status['version']
    except:
        pass
    return jsonify(res)
  
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

@app.route("/escher/cluster", methods=["POST"])
def escher_cluster_map():
    cluster_data = request.get_json()
    #with open('/Users/fliu/workspace/jupyter/data/www/annotation/data/latest_cluster.json', 'w') as f:
    #    f.write(json.dumps(cluster_data))
        
    report = generate_integration_report(cluster_data, escher_manager)
        
    return jsonify(report)

@app.route("/escher/merge", methods=["POST"])
def escher_merge_map():
    cluster_data = request.get_json()
    em_merge = merge_with_layer(cluster_data, escher_manager, modelseed_local)
        
    return jsonify(em_merge.escher_map)


@app.route("/escher/build/grid", methods=["POST"])
def build_grid_map():
    build_data = request.get_json()
    print(build_data)
    map_assembly = process_build_data_input(build_data['maps'])
    grid_x = build_data['x']
    grid_y = build_data['y']
    print(map_assembly)
    
    escher_factory = EscherFactoryApi(escher_manager)
    escher_factory.model_path = CACHE_BASE_FOLDER
    escher_factory.cpd_mapping = MODEL_CPD_MAPPING
    escher_factory.rxn_mapping = MODEL_RXN_MAPPING
    escher_factory.cmp_mapping = MODEL_CMP_MAPPING
    escher_factory.bios_cache = bios.model_data
    
    master = escher_factory.build_grid(map_assembly, (grid_x, grid_y))
    #print(content)
    #build_data = request.json
        
    #[ "iMM904;c;c;ModelSEED.NAD(P) Biosynthesis", "iMM904;c;c;ModelSEED.Pentose and Glucuronate Interconversions", "iMM904;c;c;ModelSEED.Mannitol Utilization", "iJDZ836;CCO__45__CYTOSOL;c;ModelSEED.NAD(P) Biosynthesis", "iJDZ836;CCO__45__CYTOSOL;c;ModelSEED.Pentose and Glucuronate Interconversions", "iJDZ836;CCO__45__CYTOSOL;c;ModelSEED.Mannitol Utilization" ]
    
    #base = '/Users/fliu/Library/Caches/escher/1-0-0/5/maps/'
    #m = escher_manager.get_map('ModelSEED', 'ModelSEED', 'Chorismate Synthesis')
    return jsonify(master.escher_map)


@app.route("/template/<template_id>/reaction/<rxn_id>", methods=["GET"])
def get_template_reaction(template_id, rxn_id):
    reaction_template_id = '{}@{}'.format(rxn_id, template_id)
    logger.warning("/template/%s/reaction/%s", template_id, rxn_id)
    data = annotation_api_atlas.collection_templates_reactions.find_one({'_id' : reaction_template_id})
    print(data)
    return jsonify(data)

@app.route("/template/<template_id>/functions_rxn", methods=["POST"])
def post_template_function_rxns(template_id):
    body = request.get_json()
    result = {}
    for function_id in body:
        res = annotation_api_atlas.get_rxn_with_function(function_id, template_id)
        if res == None:
            res = {}
        result[function_id] = res
    return jsonify(result)

@app.route("/template/<template_id>/reaction/<rxn_id>", methods=["POST"])
def set_annotation_to_template(template_id, rxn_id):
    #print('!!', request.json)
    #print('!!', request.data)
    print(rxn_id, request.form.get('user_id'), request.form.get('function_id'), request.form.get('logic') == 'true', type(request.form.get('logic')))

    annotation_api_atlas.add_function_to_template_rxn(
        int(request.form.get('function_id')), 
        rxn_id, 
        request.form.get('user_id'), 
        template_id, 
        request.form.get('logic'))
    return ""


@app.route("/template/<template_id>/annotation/reaction/<rxn_id>/ko/<ko_id>", methods=["POST"])
def post_template_annotation_reaction_manual_ko(template_id, rxn_id, ko_id):
    body = request.get_json()
    o = annotation_api.get_node('KeggOrthology', ko_id)
    if not o == None:
        print(o)
        annotation_api_atlas.set_manual_ko(rxn_id, template_id, ko_id, True, 'fliu')
        return jsonify(True)
    return jsonify(False)

@app.route("/template/<template_id>/annotation/reaction/<rxn_id>/status", methods=["POST"])
def get_template_annotation_reaction_status(template_id, rxn_id):
    t1 = time.time()
    
    body = request.get_json()
    reaction_template_id = '{}@{}'.format(rxn_id, template_id)
    
    #genome_set = None
    genome_set_id = None
    if 'genome_set_id' in body and len(body['genome_set_id']) > 0:
        genome_set_id = body['genome_set_id']
        #genome_set = annotation_api.get_genome_set(body['genome_set_id'])
    
    
    t2 = time.time()
    
    #print('get_template_annotation_reaction_status', 'body', body)
    #print('get_template_annotation_reaction_status', 'genome_set', genome_set)
        
    rxn_annotation_manual_function = annotation_api_atlas.get_manual_function(rxn_id, template_id)
    rxn_annotation_manual_ko = annotation_api_atlas.get_manual_ko(rxn_id, template_id)
    
    
    t3 = time.time()
    #print('get_template_annotation_reaction_status', 'rxn_annotation_manual_function', rxn_annotation_manual_function)
    #print('get_template_annotation_reaction_status', 'rxn_annotation_manual_ko', rxn_annotation_manual_ko)
    
    rxn_annotation = annotation_api.get_reaction_annotation_data3(
        rxn_id, genome_set_id, 10, 
        rxn_annotation_manual_ko['ko'],
        rxn_annotation_manual_function['functions'])
    
    t4 = time.time()
    
    rxn_annotation_curation = annotation_api_atlas.collection_templates_reactions.find_one({'_id' : reaction_template_id})
    
    rxn_annotation_functions_rxn = {}
    function_ids = set()
    for name in rxn_annotation:
        function_ids.add(rxn_annotation[name]['id'])
    for function_id in function_ids:
        res = annotation_api_atlas.get_rxn_with_function(function_id, template_id)
        if res == None:
            res = {}
        rxn_annotation_functions_rxn[function_id] = res
    
    t5 = time.time()
    
    rxn = modelseed_local.get_seed_reaction(rxn_id)
    if rxn == None:
        rxn = {}
        
    t6 = time.time()
    
    print('get_template_annotation_reaction_status::genome_set', round(t2 - t1, 6))
    print('get_template_annotation_reaction_status::atlas', round(t3 - t2, 6))
    print('get_template_annotation_reaction_status::get_reaction_annotation_data', round(t4 - t3, 6))
    print('get_template_annotation_reaction_status::atlas again', round(t5 - t4, 6))
    print('get_template_annotation_reaction_status::modelseed_local', round(t6 - t5, 6))
    print('get_template_annotation_reaction_status::TOTAL', round(t6 - t1, 6))
    
    response = {
        'rxn' : {rxn_id : clear_nan(rxn.data)},
        'cpd' : {},
        'manual_function' : rxn_annotation_manual_function,
        'manual_ko' : rxn_annotation_manual_ko,
        'annotation' : rxn_annotation,
        'curation' : rxn_annotation_curation,
        'function_rxns' : rxn_annotation_functions_rxn
    }
    
    response = json.loads(json.dumps(response))
    #print(response)
    
    return jsonify(response)

@app.route("/template/<template_id>/reaction/<rxn_id>/gene", methods=["POST"])
def post_template_reaction_gene_annotation(template_id, rxn_id):
    data = request.get_json()
    if 'genome_id' in data and 'gene_id' in data \
     and 'user_id' in data and 'logic' in data and 'desc' in data:
        annotation_api_atlas.set_annotation_to_gene(
            data['genome_id'], 
            data['gene_id'], 
            rxn_id, 
            data['user_id'], 
            template_id, 
            data['logic'],
            data['desc']
        )
    else:
        print('bad params:', data.keys())
    return ""

@app.route("/template/<template_id>/reaction/<rxn_id>/gene", methods=["GET"])
def get_template_reaction_gene_annotation(template_id, rxn_id):
    logger.warning("/template/%s/reaction/%s/gene", template_id, rxn_id)
    data = annotation_api_atlas.get_reaction_gene_annotation(rxn_id, template_id)
    return jsonify(data)

@app.route("/template/<template_id>/model_reaction/<mrxn_id>/map", methods=["POST"])
def post_template_model_reaction_mapping(template_id, mrxn_id):
    data = request.get_json()
    annotation_api_atlas.add_function_to_template_rxn(
        int(request.form.get('function_id')), 
        rxn_id, 
        request.form.get('user_id'), 
        template_id, 
        request.form.get('logic'))
    return ""

@app.route("/template/<template_id>/model_compound/<cpd_id>/map", methods=["POST"])
def post_template_model_compound_mapping(template_id, rxn_id):
    data = request.get_json()
    annotation_api_atlas.add_function_to_template_rxn(
        int(request.form.get('function_id')), 
        rxn_id, 
        request.form.get('user_id'), 
        template_id, 
        request.form.get('logic'))
    return ""

@app.route("/annotation/genome_set", methods=["GET"])
def list_genome_set():
    res = annotation_api.list_genome_sets()
    
    return jsonify(list(res))

@app.route("/annotation/genome_set/<id>", methods=["GET"])
def get_genome_set(id):
    res = annotation_api.get_genome_set(id)
    resp = {}
    if not res == None:
        resp['id'] = id
        resp['genomes'] = list(res)
    
    return jsonify(resp)

@app.route("/annotation/ko/<id>", methods=["GET"])
def get_annotation_ko(id):
    functions = annotation_api.get_functional_roles(id)
    
    resp = {}
    for k in functions:
        resp[k] = len(functions[k])
    return jsonify(resp)

@app.route("/annotation/ortholog/<seed_id>", methods=["GET"])
def get_ortholog_from_seed_rxn(seed_id):
    all_orthologs, matched_orthologs, genome_match = annotation_orth.get_orthologs_from_seed_rxn_id2(seed_id)
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

@app.route("/annotation/rxn/<id>", methods=["GET"])
def get_rxn_annotation(id):
    if id in HUGE_CACHE:
        return jsonify(HUGE_CACHE[id])
    
    resp = annotation_api.get_reaction_annotation_data(id)
    
    HUGE_CACHE[id] = resp
    
    return jsonify(resp)

@app.route("/model", methods=["GET"])
def list_models():
    
    return jsonify([])

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

@app.route("/genome/kbase/<id>", methods=["PUT"])
def put_genome_from_kbase(id):
    data = request.get_json()
    return jsonify(data)

@app.route("/query/genome", methods=["POST"])
def post_query_genome():
    data = request.get_json()
    form = result = request.form
    start = int(request.form.get('start'))
    length = int(request.form.get('length'))
    print(data)
    print(form)
    sval = request.form.get('search[value]')
    total = -1
    total_filter = 50
    taxa_filter = None
    if not sval == None and len(sval.strip()) > 0:
        taxa_filter = sval
    result = annotation_api.page_genomes(0, length, taxa_filter)
    
    rows = []
    if not result == None:
        for r in result:
            total = r['total_genomes']
            row = [
                r['n'].id, 
                r['n']['key'], 
                r['n']['scientific_name']
            ]
            rows.append(row)
            #print(r['n'].id, r['n']['key'], r['n']['scientific_name'], r['total_genomes'])
            
    result = {
        'draw' : 0,
        'recordsTotal' : total,
        'recordsFiltered' : total_filter,
        'data' : rows
    }
            
    return jsonify(result)


@app.route("/query/genome/<genome_id>/genes", methods=["POST"])
def post_query_genome_genes(genome_id):
    sval = request.form.get('search[value]')
    start = int(request.form.get('start'))
    length = int(request.form.get('length'))
    page = int(start / length)
    #print(data)
    #print(form)
    print(start, length, page)
    genes = set()
    function_filter = None
    if not sval == None and len(sval.strip()) > 0:
        function_filter = sval
    
    result = annotation_api.page_genome_genes(genome_id, start, length, function_filter)
    
    rows = []
    total = 0
    total_filter = 50
    if not result == None:
        for r in result:
            total = r['total_genes']
            total_filter = r['total_genes']
            genes.add(r['n']['key'])
            #print(r['n'].id, r['n']['key'], r['function'], r['function_source'])
            gene_id, genome_id = r['n']['key'].split('@')
            row = [
                r['n'].id, 
                genome_id, 
                gene_id, 
                [r['function'], r['function_source']], 
                {}
            ]
            rows.append(row)
            
    result = {
        'draw' : request.form.get('draw'),
        'recordsTotal' : total,
        'recordsFiltered' : total_filter,
        'data' : rows
    }
            
    return jsonify(result)

@app.route("/query/function", methods=["POST"])
def post_query_function():
    data = request.get_json()
    form = result = request.form
    print('draw:', request.form.get('draw'))
    print('search[value]:', request.form.get('search[value]'))
    print('start:', request.form.get('start'))
    print('length:', request.form.get('length'))
    
    start = int(request.form.get('start'))
    length = int(request.form.get('length'))
    sval = request.form.get('search[value]')
    sparam = ""
    if not sval == None and len(sval.strip()) > 0:
        sparam = "WHERE n.key CONTAINS '{}'".format(sval)
    print(sparam)
    page = int(start / length)
    print(start, length, page)
    a, b = annotation_api.page_nodes2('Function', 
                                      page, 
                                      length, 
                                      sparam)
    
    def get_gene_counts(function_str):
        query = 'MATCH (n:Function {key:{function_str}})-[r:has_gene]->(o:KBaseGene) RETURN count(r) as total'
        c = annotation_api.matcher.graph.run(query, function_str = function_str)
        o = c.next()
        return o.data()['total']
    
    def get_function_source(n):
        function_source = set()
        n = annotation_api.matcher.get(n.id)
        for rel in n.graph.match((n, None), r_type = 'has_source'):
            function_source.add(rel.end_node['key'])
        return function_source
    
    def get_kos(n, limit):
        kos = set()
        n = annotation_api.matcher.get(n.id)
        for rel in n.graph.match((None, n), r_type="has_annotation", limit=limit):
            node_gene = rel.start_node
            for rel_ortholog in node_gene.graph.match((node_gene, None), r_type="has_ortholog"):
                kos.add(rel_ortholog.end_node['key'])
        return kos
    
    search_limit = 2000
    rows = []
    if not a == None:
        for r in a:
            n = r['n']
            #print(n.id, n['key'])
            function_source = get_function_source(n)
            gene_count = get_gene_counts(n['key'])
            kos = get_kos(n, search_limit)
            if gene_count > search_limit:
                kos.add('*')
            row = [n.id, n['key'], gene_count, list(kos), list(function_source)]
            rows.append(row)
    #print(form)
    result = {
        'draw' : request.form.get('draw'),
        'recordsTotal' : 182576,
        'recordsFiltered' : b[0].get('count'),
        'data' : rows
    }
    return jsonify(result)

@app.route("/query/ko", methods=["POST"])
def post_query_ko():
    data = request.get_json()
    print(data)
    return jsonify(data)

#@app.route("/genome/kbase/<id>", methods=["PUT"])
#def post_query_function(id):
#    data = request.get_json()
#    return jsonify(data)

def get_functional_roles(ko_id):
    functions = {}
    doc = mdb_kbase_ko_to_genes.find_one({'_id' : ko_id})
    gene_mapping = doc['genes']
    for tup in gene_mapping:
        print(tup)
        kbase_genome = mdb_kbase_genomes.find_one({'_id' : tup[0]})
        genome = cobrakbase.core.KBaseGenome(kbase_genome)
        features = list(filter(lambda f : f['id'] == tup[1], genome.data['features']))
        if len(features) == 1:
            if not 'function' in features[0]:
                features[0]['function'] = 'null'
            if not features[0]['function'] in functions:
                functions[features[0]['function']] = set()
            functions[features[0]['function']].add((tup[0], tup[1]))
        else:
            print('error', tup)
    return functions

def get_functional_roles2(ko_id, annotation_api):
    functions = {}
    doc = mdb_kbase_ko_to_genes.find_one({'_id' : ko_id})
    if doc == None:
        return {}, {}, {}
    gene_mapping = doc['genes']
    
    found = 0
    for tup in gene_mapping:
        gene_id = "{}@{}".format(tup[1], tup[0])
        gene_functions = annotation_api.collection_gene_functions.find_one({'_id' : gene_id})
        
        if not gene_functions == None:
            found += 1
            for annotation_id in gene_functions['function']:
                function_id = gene_functions['function'][annotation_id]
                if not function_id in functions:
                    functions[function_id] = set()
                functions[function_id].add(gene_id)

    functions_data = {}
    
    for function_id in functions:
        function_doc = annotation_api.collection_functions.find_one({'_id' : function_id})
        if not function_doc == None:
            functions_data[function_id] = function_doc
    
    return functions, functions_data, {'total' : len(gene_mapping), 'has_function' : found}

def get_reaction_annotation_data(rxn_id):
    kos = annotation_api.get_ko_by_seed_id(rxn_id)
    function_data = {}
    for ko in kos:
        #print(ko)
        functions, functions_data, metadata = get_functional_roles2(ko, annotation_api)
        print(ko, metadata)
        for f in functions:
            function = functions_data[f]['value']
            if not function in function_data:
                function_data[function] = {
                    'id' : f,
                    'hits' : []
                }
            function_data[function]['hits'].append({
                'score' : len(functions[f]),
                'source' : ['KEGG', ko]
            })
    for template_doc in annotation_api.collection_templates.find():
        if rxn_id in template_doc['reactions']:
            template_id = template_doc['_id']
            template_reaction_doc = template_doc['reactions'][rxn_id]
            for and_rule in template_reaction_doc['functions']:
                for function_id in and_rule:
                    function_doc = annotation_api.collection_functions.find_one({'_id' : function_id})
                    function = function_doc['value']
                    if not function in function_data:
                        function_data[function] = {
                            'id' : function_id,
                            'hits' : []
                        }
                    function_data[function]['hits'].append({
                        'score' : 0,
                        'source' : ['template', template_id]
                    })
                    #print(template_id, function_doc['value'])
        
    for function in function_data:
        function_id = function_data[function]['id']
        #print(function_id)
        function_doc = annotation_api.collection_functions.find_one({'_id' : function_id})
        function_metadata = annotation_api.collection_functions_data.find_one({'_id' : function_id})
        
        subsystem_tags = {}
        if 'subsystem' in function_doc:
            for subsystem_tag in function_doc['subsystem']:
                subsystem_tags[subsystem_tag] = [function_doc['subsystem'][subsystem_tag]]
                if 'subsystem_class' in function_doc and subsystem_tag in function_doc['subsystem_class']:
                    subsystem_tags[subsystem_tag].append(function_doc['subsystem_class'][subsystem_tag])
        #print(subsystem_tags)
        
        
        source_tags = {}
        if not function_metadata == None:
            sources = function_metadata['sources']
            for s in sources:
                source_tags[s] = [len(sources[s]['genomes']), len(sources[s]['genes']), sources[s]['genes'][:10]]
            #print(source_tags)
        function_data[function]['sources'] = source_tags
        function_data[function]['subsystems'] = subsystem_tags
    return function_data

escher_manager = None
annotation_api = None
annotation_orth = None
MODEL_CPD_MAPPING = None
MODEL_RXN_MAPPING = None
MODEL_CMP_MAPPING = None
CACHE_BASE_FOLDER = None
bios = None

def load_cache_data(folder):
    bios = BIOS_MOCK(folder + 'bios_cache_fungi.json')
    MODEL_CPD_MAPPING = None
    MODEL_RXN_MAPPING = None
    MODEL_CMP_MAPPING = None
    with open(folder + 'cpd_mapping_cache4.json', 'r') as f:
        MODEL_CPD_MAPPING = json.loads(f.read())
    with open(folder + 'rxn_mapping_cache4.json', 'r') as f:
        MODEL_RXN_MAPPING = json.loads(f.read())
    with open(folder + 'cmp_mapping_cache.json', 'r') as f:
        MODEL_CMP_MAPPING = json.loads(f.read())
    return bios, MODEL_CMP_MAPPING, MODEL_CPD_MAPPING, MODEL_RXN_MAPPING

if __name__ == '__main__':
    with open('config.yaml', 'r') as config_h:    
        config = yaml.load(config_h, Loader=yaml.FullLoader)
        
        CACHE_BASE_FOLDER = config['cache']
        MODELSEED_FOLDER = config['modelseed']['path']
        CHEMDUST_URL = config['chemapi']
        
        #bios = biosapi.BIOS()
        bios, MODEL_CMP_MAPPING, MODEL_CPD_MAPPING, MODEL_RXN_MAPPING = load_cache_data(CACHE_BASE_FOLDER)
        """
        bios = BIOS_MOCK(CACHE_BASE_FOLDER + 'bios_cache_fungi.json')
        with open(CACHE_BASE_FOLDER + 'cpd_mapping_cache4.json', 'r') as f:
            MODEL_CPD_MAPPING = json.loads(f.read())
        with open(CACHE_BASE_FOLDER + 'rxn_mapping_cache4.json', 'r') as f:
            MODEL_RXN_MAPPING = json.loads(f.read())
        with open(CACHE_BASE_FOLDER + 'cmp_mapping_cache.json', 'r') as f:
            MODEL_CMP_MAPPING = json.loads(f.read())
        """
        
    escher_manager = modelseed_escher.EscherManager(escher)
    #aclient = pymongo.MongoClient('mongodb://192.168.1.21:27017/')
    aclient = pymongo.MongoClient(config['mongo_client'])
    #annotation_api = AnnotationApi(mclient)
    annotation_api_atlas = CurationApi(aclient)
    #mclient = pymongo.MongoClient('mongodb://192.168.1.21:27017/')
    #annotation_api_atlas = CurationApi(mclient)
    
    #####      Load ModelSEED     #####
    modelseed_local = cobrakbase.modelseed.from_local(MODELSEED_FOLDER)
    
    host, port, user, pwd = (config['neo4j']['host'], 
                             config['neo4j']['port'], 
                             config['neo4j']['user'], 
                             config['neo4j']['password'])
    if len(sys.argv) > 1:
        host = sys.argv[1]
    if len(sys.argv) > 2:
        pwd = sys.argv[2]
        
    #annotation_api = AnnotationApiNeo4j(user=user, pwd=pwd, port=port, host=host)
    cache = redis.Redis(host=config['redis']['host'], 
                        port=config['redis']['port'], 
                        db=config['redis']['db'])

    annotation_api = AnnotationApiRedisCache(cache, user=user, pwd=pwd, port=port, host=host)

    annotation_api.neo4j_graph = Graph("http://neo4j:" + pwd + "@" + host + ":7474")
    annotation_api.matcher = NodeMatcher(annotation_api.neo4j_graph)
    annotation_api.r_matcher = RelationshipMatcher(annotation_api.neo4j_graph)
    annotation_api.init_constraints()
    
    kbase = cobrakbase.KBaseAPI(config['kbase']['token'])
    #annotation_orth = build_annotation_ortholog(kbase, CACHE_BASE_FOLDER, bios)
    #annotation_orth.model_rxn_mapping = MODEL_RXN_MAPPING
    #annotation_orth.model_cpd_mapping = MODEL_CPD_MAPPING
    
    import controller_biochem
    import controller_annotation
    
    app.run(port=8058, host='0.0.0.0', debug=False)
    
