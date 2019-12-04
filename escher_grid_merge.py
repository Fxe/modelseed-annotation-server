import copy
from modelseed_escher.convert_utils import add_compartment
from modelseed_escher.core import EscherMap
from modelseed_escher.map import EscherCluster

def move_to_compartment(cmp_id, em):
    em.add_uid_to_reaction_metabolites()
    node_uid_cmp = {}
    node_uid_id = {}
    for node_uid in em.escher_graph['nodes']:
        
        n = em.escher_graph['nodes'][node_uid]
        if n['node_type'] == 'metabolite':
            if not 'compartment' in n:
                node_uid_cmp[node_uid] = cmp_id
                n['bigg_id'] += '_' + cmp_id
                node_uid_id[node_uid] = n['bigg_id']
    add_compartment(em, node_uid_cmp)
    for rxn_uid in em.escher_graph['reactions']:
        rnode = em.escher_graph['reactions'][rxn_uid]
        for o in rnode['metabolites']:
            o['bigg_id'] = node_uid_id[o['node_uid']]
            
        rnode['bigg_id'] += '_' + cmp_id
    return node_uid_cmp

def get_next_id(em):
    next_id = 0
    for node_id in em.escher_graph['nodes']:
        if int(node_id) >= next_id:
            next_id = int(node_id) + 1

    return next_id

def add_nodes(em1, em2, next_id, x_offset = 0, y_offset = 0):
    for uid in em1.nodes:
        n = em1.nodes[uid]
        if n['node_type'] == 'metabolite':
            n_copy = copy.deepcopy(n)
            n_copy['x'] += x_offset
            n_copy['y'] += y_offset
            n_copy['label_x'] += x_offset
            n_copy['label_y'] += y_offset
            em2.escher_graph['nodes'][next_id] = n_copy
            next_id += 1
    return next_id

def add_layers(em, cluster_grid, escher_manager):
    #em = EscherMap(cluster_params['escher_map'])
    next_id = get_next_id(em)
    for map_id in cluster_grid:
        if '.' in map_id:
            a, b = map_id.split('.', 1)
            for grid_cell in cluster_grid[map_id]:
                seed_em = escher_manager.get_map('ModelSEED', a, b)
                canvas = seed_em.escher_graph['canvas']
                x, y, w, h, cmp = grid_cell
                x_offset = (x * w) + -1 * canvas['x']
                y_offset = (y * h) + -1 * canvas['y']
                seed_em_cmp = seed_em.clone()
                move_to_compartment(cmp, seed_em_cmp)
                next_id = add_nodes(seed_em, em, next_id, 
                                    x_offset, y_offset)
                next_id = add_nodes(seed_em_cmp, em, next_id, 
                                    x_offset, y_offset)
                print(next_id, map_id, x, y, w, h, cmp)
                        
    return em

def report(cluster_data, em):
    cluster_map = em.escher_map
    any_merge = True
    result = {
        'models' : set(),
        'records' : []
    }
    for cluster in cluster_data:
        record = {
            'model_data' : {},
            'database' : {}
        }
        match_cmp = '?'
        cmps = set()
        for id in cluster:
            node_id = cluster_map[1]['nodes'][id]['bigg_id']
            if node_id.startswith('cpd') and '_' in node_id:
                cpd_id, cmp = node_id.split('_')
                cmps.add(cmp)

        if len(cmps) == 1:
            match_cmp = cmps.pop()
        for id in cluster:
            node_id = cluster_map[1]['nodes'][id]['bigg_id']
            cpd_id = node_id
            database = None
            if '@' in node_id:
                cpd_id, database = node_id.split('@')

            #print(cpd_id, database)
            #if is model
            if not database == None:
                result['models'].add(database)
                cmp = match_cmp
                if not cmp in record['model_data']:
                    record['model_data'][cmp] = {}
                if not database in record['model_data'][cmp]:
                    record['model_data'][cmp][database] = set()
                record['model_data'][cmp][database].add(cpd_id)
            else:
                database = 'seed.compound'
                if not database in record['database']:
                    record['database'][database] = set()
                if node_id.startswith('cpd') and '_' in node_id:
                    cpd_id, cmp = node_id.split('_')
                record['database'][database].add(cpd_id)
            #elseif database
        if len(record['model_data']) > 0:
            result['records'].append(record)


    result['models'] = list(result['models'])
    for r in result['records']:
        for cmp_id in r['model_data']:
            for model_id in r['model_data'][cmp_id]:
                r['model_data'][cmp_id][model_id] = list(r['model_data'][cmp_id][model_id])
        for database_id in r['database']:
            r['database'][database_id] = list(r['database'][database_id])
            
    return result

def generate_integration_report(cluster_params, escher_manager):
    em = EscherMap(cluster_params['escher_map'])
    em = add_layers(em, cluster_params['grid'], escher_manager)
    ec = EscherCluster(25)
    cluster_data = ec.cluster(em)
    result = report(cluster_data, em)
    return result