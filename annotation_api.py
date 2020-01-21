import logging
import copy
import time
import hashlib
import functools
import operator
import pymongo

logger = logging.getLogger(__name__)

def convert_to_mongo_document(id, genome):
    genome_doc = {}
    #genome_doc['workspace_id'] = 'jplfaria:narrative_1524466549180'
    genome_doc['_id'] = id
    for k in genome.data:
        if not type(genome.data[k]) == list:
            genome_doc[k] = genome.data[k]
    return genome_doc

def update_genome_genes(id, genome):
    features = {}
    feature_seq_dna = {}
    feature_seq_protein = {}
    for f in genome.data['features']:
        feature_id = "{}@{}".format(f['id'], id)
        features[feature_id] = copy.deepcopy(f)
        features[feature_id]['_id'] = feature_id
        dna = None
        protein = None
        
        if 'function' in f:
            del features[feature_id]['function']
        if 'dna_sequence' in f:
            dna = f['dna_sequence']
            del features[feature_id]['dna_sequence']
        if 'protein_translation' in f:
            protein = f['protein_translation']
            del features[feature_id]['protein_translation']
        
        feature_seq_dna[feature_id] = {
            '_id' : feature_id,
            'sequence' : dna
        }
        feature_seq_protein[feature_id] = {
            '_id' : feature_id,
            'sequence' : protein
        }
    return features, feature_seq_dna, feature_seq_protein

def load_seed_reaction(ms, kegg_rxn_to_ko, annotation_api):
    for rxn_id in ms.reactions:
        kegg_ids = set()
        kegg_kos = set()
        rxn = ms.get_seed_reaction(rxn_id)
        if rxn_id in ms.reaction_aliases and 'KEGG' in ms.reaction_aliases[rxn_id]:
            kegg_ids = set(ms.reaction_aliases[rxn_id]['KEGG'])
        
        to_reduce = list(map(lambda x : list(kegg_rxn_to_ko[x]) if x in kegg_rxn_to_ko else [], kegg_ids))
        if len(to_reduce) == 0 or len(to_reduce) == 1 and len(to_reduce[0]) == 0:
            kegg_kos = []
        else:
            kegg_kos = set(functools.reduce(operator.add, to_reduce))
        #if len(kegg_kos) > 0:
        #print(rxn_id, kegg_ids, kegg_kos)

        seed_rxn_doc = {
            '_id' : rxn_id,
            'name' : rxn.data['name'] if 'name' in rxn.data else "",
            'kegg_ids' : list(kegg_ids),
            'kegg_kos' : list(kegg_kos)
        }
        
        annotation_api.collection_reactions.update_one({'_id' : rxn_id}, {'$set' : seed_rxn_doc}, upsert=True)

class AnnotationApi:
    
    def __init__(self, client, database = 'annotation'):
        self.mongodb = client[database]
        
        self.collection_reactions = self.mongodb['reactions']
        
        self.collection_functions = self.mongodb['functions_group']

        self.collection_functions = self.mongodb['functions']
        self.collection_functions_data = self.mongodb['functions_data']
        
        self.collection_templates = self.mongodb['templates']
        self.collection_templates_reactions = self.mongodb['templates_reactions']
        
        self.collection_genes = self.mongodb['genes']
        self.collection_genes_seq_dna = self.mongodb['genes_seq_dna']
        self.collection_genes_seq_protein = self.mongodb['genes_seq_protein']
        self.collection_genomes = self.mongodb['genomes']
        self.collection_genomes_genes = self.mongodb['genomes_genes']
        self.collection_gene_functions = self.mongodb['gene_functions']
        
    def get_ko_by_seed_id(self, seed_id):
        rxn_doc = self.collection_reactions.find_one({'_id' : seed_id})
        print(rxn_doc)
        if not rxn_doc == None:
            if 'kegg_kos' in rxn_doc:
                return set(rxn_doc['kegg_kos'])
            else:
                return set()
        return None
    
    def add_genome(self, genome_id, genome):
        genome_doc = self.collection_genomes.find_one({'_id' : genome_id})
        
        if genome_doc == None:
            genome_doc = convert_to_mongo_document(genome_id, genome.data)
            features, features_dna, features_protein = update_genome_genes(genome_id, genome)
            genome_doc['count_features'] = len(features)
            
            self.collection_genomes.insert_one(genome_doc)
            self.collection_genomes_genes.insert_one({
                '_id' : genome_id,
                'features' : list(map(lambda o : o['id'], genome.data['features']))
            })
            self.collection_genes.insert_many(list(features.values()))
            self.collection_genes_seq_dna.insert_many(list(features_dna.values()))
            self.collection_genes_seq_protein.insert_many(list(features_protein.values()))
            
            print(genome_id)
        else:
            print('found')
        return genome_doc
    
    def hash_value(self, value):
        m = hashlib.sha256()
        m.update(value.encode())
        return m.hexdigest()
    
    def make_hash_doc(self, doc):
        hash_doc = copy.deepcopy(doc)
        hash_doc['value'] = hash_doc['_id']
        
        m = hashlib.sha256()
        m.update(doc['_id'].encode())
        hash_doc['_id'] = m.hexdigest()
        
        return hash_doc
    
    def get_function(self, function_name):
        function_doc = {'_id' : function_name}
        function_doc = self.make_hash_doc(function_doc)
        collection_doc = self.collection_functions.find_one(function_doc)
        return collection_doc
    
    def add_function(self, function_name):
        function_doc = {'_id' : function_name}
        function_doc = self.make_hash_doc(function_doc)
        collection_doc = self.collection_functions.find_one(function_doc)
        if collection_doc == None:
            self.collection_functions.insert_one(function_doc)
        return function_doc
    
    def add_functions(self, function_names):
        docs = map(lambda f : {'_id' : f}, function_names)
        hash_docs = map(lambda doc : self.make_hash_doc(doc), docs)
        for function_doc in hash_docs:
            collection_doc = self.collection_functions.find_one(function_doc)
            if collection_doc == None:
                self.collection_functions.insert_one(function_doc)
        return hash_docs
    
    def add_annotation_from_genome(self, annotation_id, genome_id, genome):
        feature_to_function = {}
        for f in genome.data['features']:
            feature_id = "{}@{}".format(f['id'], genome_id)
            if 'function' in f:
                feature_to_function[feature_id] = f['function']
                self.add_gene_annotation(annotation_id, feature_id, f['function'])
        return feature_to_function
    
    def add_gene_annotation(self, annotation_id, feature_id, function):
        function_doc = self.add_function(function)
        gene_function_doc = self.collection_gene_functions.update_one(
            {'_id' : feature_id},
            {'$set' : {'function.' + annotation_id : function_doc['_id']}}, upsert= True
        )
        
        return gene_function_doc
    
    def add_template(self, template_id, template):
        roles = {}
        function_uids = {}
        compcompound_compartment = {}

        for compcompound in template['compcompounds']:
            compcompound_compartment[compcompound['id']] = compcompound['templatecompartment_ref'].split('/')[-1]

        for role in template['roles']:
            roles[role['id']] = role['name']
        for role_id in roles:
            function = roles[role_id]
            #print(function)
            function_doc = self.get_function(function)
            if function_doc['value'] == function:
                function_uids[function] = function_doc['_id']
            else:
                print('errooo!')
        missing_function = set()
        for role_id in roles:
            function = roles[role_id]
            if not function in function_uids:
                missing_function.add(function)

        print(len(missing_function))

        cpx_to_function_uids = {}
        for cpx in template['complexes']:
            uids = set()
            for complexrole in cpx['complexroles']:
                role_id = complexrole['templaterole_ref'].split('/')[-1]
                if role_id in roles:
                    uids.add(function_uids[roles[role_id]])
                else:
                    print('errooo!')

            cpx_to_function_uids[cpx['id']] = uids

        reactions = {}
        for rxn in template['reactions']:
            rxn_id = rxn['reaction_ref'].split('/')[-1]
            rxn_compartments = set()
            for templateReactionReagent in rxn['templateReactionReagents']:
                templatecompcompound_ref = templateReactionReagent['templatecompcompound_ref'].split('/')[-1]
                rxn_compartments.add(compcompound_compartment[templatecompcompound_ref])
            or_rule = set()
            for complex_ref in rxn['templatecomplex_refs']:
                complex_id = complex_ref.split('/')[-1]
                and_rule = set()
                if complex_id in cpx_to_function_uids:
                    and_rule = frozenset(cpx_to_function_uids[complex_id])
                or_rule.add(and_rule)

            reactions[rxn_id] = {
                'functions' : list(map(lambda x : list(x), or_rule)),
                'compartment' : list(rxn_compartments),
                'base_cost' : rxn['base_cost'],
                'direction' : rxn['direction'],
                'forward_penalty' : rxn['forward_penalty'],
                'reverse_penalty' : rxn['reverse_penalty'],
                'type' : rxn['type'],
                'GapfillDirection' : rxn['GapfillDirection'],

            }

        template_doc = {
            '_id' : template_id,
            'reactions' : reactions
        }

        self.collection_templates.insert_one(template_doc)

        return template_doc
    
    def index_gene_functions(self, gene_genome_id):
        gene_functions_doc = self.collection_gene_functions.find_one({'_id' : gene_genome_id})
        if not gene_functions_doc == None:
            gene_id, genome_id = gene_functions_doc['_id'].split('@')
            for source in gene_functions_doc['function']:
                function_id = gene_functions_doc['function'][source]
                functions_data_doc = self.collection_functions_data.find_one({'_id' : function_id})
                if functions_data_doc == None:
                    self.collection_functions_data.insert_one({
                        '_id' : function_id,
                        'sources' : {},
                        'sub_functions': [],
                    })
                self.collection_functions_data.update_one(
                    {'_id' : function_id},
                    {'$addToSet': {
                            'sources.' + source + '.genomes': genome_id,
                            'sources.' + source + '.genes': gene_functions_doc['_id'],
                        }
                    }
                )
        else:
            print('!')

            
    def get_template_reaction_function_score(self, reaction_id, template_id, function_id):
        reaction_template_id = '{}@{}'.format(reaction_id, template_id)
        o = self.collection_templates_reactions.find_one({'_id' : reaction_template_id})
        #print(o['functions'])
        if not o == None and 'functions' in o and str(function_id) in o['functions']:
            return o['functions'][str(function_id)]
        return None
            
    def add_function_to_template_rxn(self, function_id, reaction_id, user_id, template_id, logic):
        reaction_template_id = '{}@{}'.format(reaction_id, template_id)

        timestamp = int(time.time())

        #update template
        doc = self.collection_templates_reactions.find_one({'_id':reaction_template_id})
        if doc == None:
            self.collection_templates_reactions.insert_one({
                '_id' : reaction_template_id,
                'functions' : {},
                'log' : []
            })
        self.collection_templates_reactions.update_one(
            {"_id" :reaction_template_id},
            {'$set' : {"functions." + str(function_id) : logic}})

        #log action
        self.collection_templates_reactions.update_one(
            {'_id' : reaction_template_id}, 
            {'$push' : {'log' : {'timestamp' : timestamp, 'user_id' : user_id, 'action' : logic, 'target' : function_id}}}, upsert=True)
        
