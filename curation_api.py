import logging
import time

logger = logging.getLogger(__name__)

def fix_mongo_object_key(doc):
    if doc == None:
        return doc
    for k in doc:
        delete = set()
        if type(doc[k]) == dict:
            for key in doc[k]:
                if '#' in key:
                    delete.add(key)
            for key in delete:
                doc[k][key.replace('#', '.')] = doc[k][key]
            for key in delete:
                del doc[k][key]
    return doc

class CurationApi:

    def __init__(self, client, database_id='annotation'):
        self.database = client[database_id]
        self.collection_templates_reactions = self.database['templates_reactions']
        self.collection_templates_reactions_ko = self.database['templates_reactions_manual_ko']
        self.collection_templates_reactions_function = self.database['templates_reactions_manual_function']
        self.collection_reaction_gene_annotation = self.database['reaction_gene_annotation']
        self.collection_model_reaction_mapping = self.database['model_reaction_mapping']
        self.collection_model_compound_mapping = self.database['model_compound_mapping']
        
    def server_info(self):
        return self.database.client.server_info()
    
    def initialize_empty_template_reaction_record(self, reaction_template_id):
        self.collection_templates_reactions.insert_one({
            '_id' : reaction_template_id,
            'functions' : {},
            'log' : [],
            'comments' : [],
            'attributes' : {}
        })
    
    def add_function_to_template_rxn(self, function_id, reaction_id, user_id, template_id, logic):
        reaction_template_id = '{}@{}'.format(reaction_id, template_id)

        timestamp = int(time.time())

        #update template
        doc = self.collection_templates_reactions.find_one({'_id':reaction_template_id})
        if doc == None:
            self.initialize_empty_template_reaction_record(reaction_template_id)
            
        self.collection_templates_reactions.update_one(
            {"_id" :reaction_template_id},
            {'$set' : {"functions." + str(function_id) : logic}})

        #log action
        self.collection_templates_reactions.update_one(
            {'_id' : reaction_template_id}, 
            {'$push' : {'log' : {'timestamp' : timestamp, 'user_id' : user_id, 'action' : logic, 'target' : function_id}}}, upsert=True)
        
        
    def set_annotation_to_gene(self, genome_id, gene_id, reaction_id, user_id, template_id, logic, comment = ""):
        reaction_template_id = '{}@{}'.format(reaction_id, template_id)
        gene_genome_id = '{}@{}'.format(gene_id, genome_id)

        timestamp = int(time.time())

        doc = self.collection_reaction_gene_annotation.find_one({'_id':reaction_template_id})
        if doc == None:
            self.collection_reaction_gene_annotation.insert_one({
                '_id' : reaction_template_id,
                'genes' : {},
                'log' : []
            })

        self.collection_reaction_gene_annotation.update_one(
            {"_id" :reaction_template_id},
            {'$set' : {"genes." + str(gene_genome_id.replace('.', '#')) : logic}}
        )

        #log action
        self.collection_reaction_gene_annotation.update_one(
            {'_id' : reaction_template_id}, 
            {'$push' : {'log' : {
                'timestamp' : timestamp, 
                'user_id' : user_id, 
                'action' : logic, 
                'comment' : comment,
                'target' : gene_genome_id}
                       }}, 
            upsert=True
        )
    
    def set_reference_to_model_reaction(self, database, ref_id, model_reaction_id, user_id, template_id, logic, comment = ""):
        reaction_template_id = '{}@{}'.format(model_reaction_id, template_id)
        #database = database.replace('.', '#')
        ref = '{}@{}'.format(ref_id, database)

        timestamp = int(time.time())

        doc = self.collection_model_reaction_mapping.find_one({'_id':reaction_template_id})
        if doc == None:
            self.collection_model_reaction_mapping.insert_one({
                '_id' : reaction_template_id,
                'mapping' : {},
                'log' : []
            })

        self.collection_model_reaction_mapping.update_one(
            {"_id" :reaction_template_id},
            {'$set' : {"mapping." + str(ref.replace('.', '#')) : logic}}
        )

        #log action
        self.collection_model_reaction_mapping.update_one(
            {'_id' : reaction_template_id}, 
            {'$push' : {'log' : {
                'timestamp' : timestamp, 
                'user_id' : user_id, 
                'action' : logic, 
                'comment' : comment,
                'target' : ref}
                       }}, 
            upsert=True
        )

    def set_reference_to_model_compound(self, database, ref_id, model_compound_id, user_id, template_id, logic, comment = ""):
        compound_template_id = '{}@{}'.format(model_compound_id, template_id)
        #database = database.replace('.', '#')
        ref = '{}@{}'.format(ref_id, database)

        timestamp = int(time.time())

        doc = self.collection_model_compound_mapping.find_one({'_id':compound_template_id})
        if doc == None:
            self.collection_model_compound_mapping.insert_one({
                '_id' : compound_template_id,
                'mapping' : {},
                'log' : []
            })

        self.collection_model_compound_mapping.update_one(
            {"_id" :compound_template_id},
            {'$set' : {"mapping." + str(ref.replace('.', '#')) : logic}}
        )

        #log action
        self.collection_model_compound_mapping.update_one(
            {'_id' : compound_template_id}, 
            {'$push' : {'log' : {
                'timestamp' : timestamp, 
                'user_id' : user_id, 
                'action' : logic, 
                'comment' : comment,
                'target' : ref}
                       }}, 
            upsert=True
        )
        
    def get_reaction_gene_annotation(self, reaction_id, template_id):
        reaction_template_id = '{}@{}'.format(reaction_id, template_id)
        data = self.collection_reaction_gene_annotation.find_one({'_id' : reaction_template_id})
        fix_mongo_object_key(data)
        return data
    
    def get_model_compound_mapping(self, model_compound_id, template_id):
        compound_template_id = '{}@{}'.format(model_compound_id, template_id)
        data = self.collection_model_compound_mapping.find_one({'_id' : compound_template_id})
        fix_mongo_object_key(data)
        return data
    
    def get_model_reaction_mapping(self, model_reaction_id, template_id):
        reaction_template_id = '{}@{}'.format(model_reaction_id, template_id)
        data = self.collection_model_reaction_mapping.find_one({'_id' : reaction_template_id})
        fix_mongo_object_key(data)
        return data
    
    
    
    def get_rxn_with_function(self, function_id, template_id):
        result = {}
        doc = self.database['template_' + template_id].find_one({'_id' : str(function_id)})
        if not doc == None and 'mapping' in doc:
            result = doc['mapping']
        #for doc in self.collection_templates_reactions.find():
        #    doc_rxn_id, doc_template_id = doc['_id'].split('@')
        #    if doc_template_id == template_id:
        #        if 'functions' in doc and str(function_id) in doc['functions']:
        #            result[doc_rxn_id] = doc['functions'][str(function_id)]
        return result
    
    def get_manual_ko(self, reaction_id, template_id):
        reaction_template_id = '{}@{}'.format(reaction_id, template_id)
        data = self.collection_templates_reactions_ko.find_one({'_id' : reaction_template_id})
        if data == None:
            return {
                '_id' : reaction_template_id,
                'ko' : {},
                'log' : []
            }
        return data
    
    def set_manual_ko(self, reaction_id, template_id, ko_id, logic, user_id):
        reaction_template_id = '{}@{}'.format(reaction_id, template_id)
        timestamp = int(time.time())

        #update template
        doc = self.collection_templates_reactions_ko.find_one({'_id':reaction_template_id})
        if doc == None:
            self.collection_templates_reactions_ko.insert_one({
                '_id' : reaction_template_id,
                'ko' : {},
                'log' : []
            })
        self.collection_templates_reactions_ko.update_one(
            {"_id" :reaction_template_id},
            {'$set' : {"ko." + str(ko_id) : logic}})

        #log action
        self.collection_templates_reactions_ko.update_one(
            {'_id' : reaction_template_id}, 
            {'$push' : {'log' : {'timestamp' : timestamp, 'user_id' : user_id, 'action' : logic, 'target' : ko_id}}}, upsert=True)
        
    def get_manual_function(self, reaction_id, template_id):
        reaction_template_id = '{}@{}'.format(reaction_id, template_id)
        data = self.collection_templates_reactions_function.find_one({'_id' : reaction_template_id})
        if data == None:
            return {
                '_id' : reaction_template_id,
                'functions' : {},
                'log' : []
            }
        return data
        
    def set_manual_function(self, reaction_id, template_id, function_id, logic, user_id):
        reaction_template_id = '{}@{}'.format(reaction_id, template_id)
        timestamp = int(time.time())

        #update template
        doc = self.collection_templates_reactions_function.find_one({'_id':reaction_template_id})
        if doc == None:
            self.collection_templates_reactions_function.insert_one({
                '_id' : reaction_template_id,
                'functions' : {},
                'log' : []
            })
        self.collection_templates_reactions_function.update_one(
            {"_id" :reaction_template_id},
            {'$set' : {"functions." + str(function_id) : logic}})

        #log action
        self.collection_templates_reactions_function.update_one(
            {'_id' : reaction_template_id}, 
            {'$push' : {'log' : {'timestamp' : timestamp, 'user_id' : user_id, 'action' : logic, 'target' : function_id}}}, upsert=True)
        
    def add_template_reaction_comment(self, reaction_id, template_id, user_id, comment):
        reaction_template_id = '{}@{}'.format(reaction_id, template_id)
        timestamp = int(time.time())
        
        doc = self.collection_templates_reactions.find_one({'_id' : reaction_template_id})
        if doc == None:
            self.initialize_empty_template_reaction_record(reaction_template_id)
                
        self.collection_templates_reactions.update_one(
            {'_id' : reaction_template_id},
            {'$push' : {'comments' : {'timestamp' : timestamp, 'user_id' : user_id, 'comment' : comment}}}, upsert=True)
        return True
        
    def get_template_reaction_comment(self, reaction_id, template_id):
        
        reaction_template_id = '{}@{}'.format(reaction_id, template_id)
        
        doc = self.collection_templates_reactions.find_one({'_id' : reaction_template_id})
        comments = {}
        if not doc == None and 'comments' in doc:
            comments = doc['comments']
        
        return comments
    
    def add_template_reaction_attribute(self, reaction_id, template_id, attribute, value):
        reaction_template_id = '{}@{}'.format(reaction_id, template_id)
        timestamp = int(time.time())
        
        doc = self.collection_templates_reactions.find_one({'_id' : reaction_template_id})
        if doc == None:
            self.initialize_empty_template_reaction_record(reaction_template_id)

        self.collection_templates_reactions.update_one(
            {'_id' : reaction_template_id},
            {'$set' : {'attributes.' +  attribute: value}}, upsert=True)
        
        return True

        
    def get_template_reaction_attributes(self, reaction_id, template_id):
        
        reaction_template_id = '{}@{}'.format(reaction_id, template_id)
        
        doc = self.collection_templates_reactions.find_one({'_id' : reaction_template_id})
        attributes = {}
        if not doc == None:
            attributes = doc['attributes']
        
        return attributes