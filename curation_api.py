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

    def __init__(self, client, database_id = 'annotation'):
        self.database = client[database_id]
        self.collection_templates_reactions = self.database['templates_reactions']
        self.collection_reaction_gene_annotation = self.database['reaction_gene_annotation']
        self.collection_model_reaction_mapping = self.database['model_reaction_mapping']
        self.collection_model_compound_mapping = self.database['model_compound_mapping']
        
    def server_info(self):
        return self.database.client.server_info()
    
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