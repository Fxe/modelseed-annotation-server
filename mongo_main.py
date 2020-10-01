from cobrakbase.core.utils import get_cmp_token


def backup(from_database, to_database):
    for collection_name in from_database.list_collections():
        print('backup:', collection_name['name'])
        for doc in from_database[collection_name['name']].find():
            to_database[collection_name['name']].insert_one(doc)


def fixxxxx(database):
    for doc in database['templates_reactions'].find():
        doc_rxn_id, doc_template_id = doc['_id'].split('@')
        if 'functions' in doc:
            for function_id in doc['functions']:
                #print(doc_template_id, doc_rxn_id, function_id, doc['functions'][function_id])
                function_rxn_doc = database['template_' + doc_template_id].find_one({'_id':function_id})
                if function_rxn_doc == None:
                    database['template_' + doc_template_id].insert_one({'_id':function_id, 'mapping' : {}})
                database['template_' + doc_template_id].update_one(
                {"_id" :function_id},
                {'$set' : {"mapping." + doc_rxn_id : doc['functions'][function_id]}}
                )


def fix_database_reaction_annotation(local_database):
    for doc in local_database['templates_reactions'].find():
        rxn_id, template_id = doc['_id'].split('@')
        local_database['templates_reactions'].update_one(
            {"_id" :doc['_id']},
            {'$set' : {"annotation.seed__DOT__reaction": rxn_id}}
        )
        

def replace_templates_reactions_cmp_token(database):
    delete = set()
    for doc in database['templates_reactions'].find():
        if 'cmp' in doc:
            delete.add(doc['_id'])
    print(len(delete))
    replace = {}
    for old_id in delete:
        doc = database['templates_reactions'].find_one({'_id' : old_id})
        rxn_id, template_id = doc['_id'].split('@')
        seed_id = doc['annotation']['seed__DOT__reaction']
        if 'cmp' in doc:
            cmp_config = doc['cmp']
            new_id = "{}_{}@{}".format(seed_id, 
                                       get_cmp_token(set(cmp_config.values())), 
                                       template_id)
            doc['_id'] = new_id
            replace[old_id] = new_id
            database['templates_reactions'].insert_one(doc)
    for old_id in delete:
        database['templates_reactions'].find_one_and_delete({'_id' : old_id})
    return replace


def replace_templates_reactions_ko_cmp_token(database, replace):
    manual_ko_set = set()
    for doc in database['templates_reactions_manual_ko'].find():
        manual_ko_set.add(doc['_id'])
    print(len(manual_ko_set))
    for old_id in replace:
        if old_id in manual_ko_set:
            doc = database['templates_reactions_manual_ko'].find_one({'_id' : old_id})
            if doc is not None:
                doc['_id'] = replace[old_id]
                database['templates_reactions_manual_ko'].insert_one(doc)
                database['templates_reactions_manual_ko'].find_one_and_delete({'_id' : old_id})
