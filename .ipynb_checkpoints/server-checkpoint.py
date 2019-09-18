import pymongo
import cobrakbase
from annotation_api import AnnotationApi
from flask import Flask, request, jsonify
from flask_restful import Resource, Api

#class (Resource):
#    def get(self):
#        return {}
    
app = Flask(__name__)
api = Api(app)
#api.add_resource(Employees, '/employees')

@app.route("/biochem/calculator/<inchi>/svg", methods=["GET"])
def translate_inchi_svg(inchi):
    return ""

@app.route("/biochem/cpd/<id>/svg", methods=["GET"])
def get_compound_svg(id):
    return ""

@app.route("/template/<template_id>/reaction/<rxn_id>", methods=["GET"])
def get_template_reaction(template_id, rxn_id):
    reaction_template_id = '{}@{}'.format(rxn_id, template_id)
    data = annotation_api_atlas.collection_templates_reactions.find_one({'_id' : reaction_template_id})
    return jsonify(data)

@app.route("/template/<template_id>/reaction/<rxn_id>", methods=["POST"])
def set_annotation_to_template(template_id, rxn_id):
    #print('!!', request.json)
    #print('!!', request.data)
    print(rxn_id, request.form.get('function_id'), request.form.get('logic') == 'true', type(request.form.get('logic')))

    annotation_api_atlas.add_function_to_template_rxn(
        request.form.get('function_id'), 
        rxn_id, 
        request.form.get('user_id'), 
        template_id, 
        request.form.get('logic') == 'true')
    return ""

@app.route("/annotation/ko/<id>", methods=["GET"])
def get_annotation_ko(id):
    functions = get_functional_roles(id)
    
    resp = {}
    for k in functions:
        resp[k] = len(functions[k])
    return jsonify(resp)

@app.route("/annotation/rxn/<id>", methods=["GET"])
def get_rxn_annotation(id):
    resp = get_reaction_annotation_data(id)
    return jsonify(resp)

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

if __name__ == '__main__':
    mclient = pymongo.MongoClient('mongodb://127.0.0.1:27017/')
    #mclient = pymongo.MongoClient('mongodb://192.168.1.21:27017/')
    #mongodb+srv://<username>:<password>@bios-dk66o.gcp.mongodb.net/test?retryWrites=true&w=majority
    #aclient = pymongo.MongoClient("mongodb+srv://server:dx75S3HBXX6h2U3D@bios-dk66o.gcp.mongodb.net/test?retryWrites=true&w=majority")
    annotation_api = AnnotationApi(mclient)
    annotation_api_atlas = AnnotationApi(mclient)
    database = mclient['KBase']
    mdb_kbase_ko_to_genes = database['ko_gene_mapping']
    mdb_kbase_genomes = database['genomes']
    mdb_kbase_taxa = database['taxa']
    
    app.run(port=8058, host='0.0.0.0')