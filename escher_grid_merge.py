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



def add_nodes2(em1, em2, next_id, x_offset = 0, y_offset = 0):
    for uid in em1.nodes:
        n = em1.nodes[uid]
        if n['node_type'] == 'midmarker':
            n_copy = copy.deepcopy(n)
            n_copy['x'] += x_offset
            n_copy['y'] += y_offset
            #n_copy['label_x'] += x_offset
            #n_copy['label_y'] += y_offset
            em2.escher_graph['nodes'][next_id] = n_copy
            next_id += 1
    return next_id

def get_reaction_midmarker_uid(uid, em):
    for rxn_uid in em.escher_graph['reactions']:
        rxn_node = em.escher_graph['reactions'][rxn_uid]
        #print(rxn_node)
        for segment_uid in rxn_node['segments']:
            segment = rxn_node['segments'][segment_uid]
            if uid == segment['from_node_id']:
                return rxn_node
            if uid == segment['to_node_id']:
                return rxn_node
    return None

def tag_midmarker_reaction(em, cmp = None):
    for uid in em.nodes:
        n = em.nodes[uid]
        if n['node_type'] == 'midmarker':
            rxn_node = get_reaction_midmarker_uid(uid, em)
            if not rxn_node == None:
                #print(uid, rxn_node['bigg_id'])
                n['rxn_id'] = rxn_node['bigg_id']
                if not cmp == None:
                    n['compartment'] = cmp

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
                tag_midmarker_reaction(seed_em, cmp)
                x_offset = (x * w) + -1 * canvas['x']
                y_offset = (y * h) + -1 * canvas['y']
                seed_em_cmp = seed_em.clone()
                move_to_compartment(cmp, seed_em_cmp)
                next_id = add_nodes(seed_em, em, next_id, 
                                    x_offset, y_offset)
                next_id = add_nodes2(seed_em, em, next_id, 
                                    x_offset, y_offset)
                next_id = add_nodes(seed_em_cmp, em, next_id, 
                                    x_offset, y_offset)
                next_id = add_nodes2(seed_em_cmp, em, next_id, 
                                    x_offset, y_offset)
                print(next_id, map_id, x, y, w, h, cmp)
                        
    return em

def sort_by_database(report, db):
    m = {}
    records = []
    for r in report['records']:
        if 'seed.compound' in r['database']:
            key = tuple(sorted(r['database']['seed.compound']))
            if not key in m:
                m[key] = {}
            for cmp in r['model_data']:
                if not cmp in m[key]:
                    m[key][cmp] = {}
                for model_id in r['model_data'][cmp]:
                    if not model_id in m[key][cmp]:
                        m[key][cmp][model_id] = set()
                    m[key][cmp][model_id] |= set(r['model_data'][cmp][model_id])
        else:
            records.append(r)
    for key in m:
        record = {
            'model_data' : {},
            'database' : {db : list(key)}
        }
        for cmp in m[key]:
            record['model_data'][cmp] = {}
            for model_id in m[key][cmp]:
                record['model_data'][cmp][model_id] = list(m[key][cmp][model_id])
        records.append(record)
        #print(key)
    report['records'] = records
    return report

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
    #cluster_ids = ec.cluster_ids(cluster_data, em)
    #cluster_uids = ec.ids_to_uid(cluster_ids, em)
    result = report(cluster_data, em)
    result = sort_by_database(result, 'seed.compound')
    return result

def get_seed_and_cmp(cc, em_grid):
    cmp = None
    seed_id = None
    for uid in cc:
        n = em_grid.nodes[uid]
        if n['bigg_id'].startswith('cpd'):
            seed_id = n['bigg_id']
            if n['bigg_id'][-2:][0] == '_':
                #print(seed_id[-2:][0], seed_id[:-2], seed_id[-1:])
                seed_id = n['bigg_id'][:-2]
                cmp = n['bigg_id'][-1:]
        #print(n['bigg_id'])
        #print(n['bigg_id'], seed_id, cmp)
    return seed_id, cmp

def get_midmarker_uid(rxn_node, em):
    for segment_uid in rxn_node['segments']:
        segment = rxn_node['segments'][segment_uid]
        if em.nodes[segment['from_node_id']]['node_type'] == 'midmarker':
            return segment['from_node_id']
        if em.nodes[segment['to_node_id']]['node_type'] == 'midmarker':
            return segment['to_node_id']
    return None

def merge_with_layer(cluster_params, escher_manager, ms = None):
    em = EscherMap(cluster_params['escher_map'])
    em_merge = em.clone()
    
    em = add_layers(em, cluster_params['grid'], escher_manager)
    ec = EscherCluster(25)
    cpd_cluster = ec.cluster(em)
    rxn_clusters = ec.compute_all_metabolite_clusters(em.escher_graph, 'midmarker')
    uid_to_cluster = {}
    rxn_uid_to_cluster = {}
    for uid_set in cpd_cluster:
        for uid in uid_set:
            uid_to_cluster[uid] = uid_set
    for uid_set in rxn_clusters:
        for uid in uid_set:
            rxn_uid_to_cluster[uid] = uid_set
            
    cpd_remap = {}
    rxn_remap = {}
    for uid in em_merge.nodes:
        n = em_merge.nodes[uid]
        if n['node_type'] == 'metabolite':
            #print(n['bigg_id'])
            if uid in uid_to_cluster:
                seed_id, cmp = get_seed_and_cmp(uid_to_cluster[uid], em)
                if not ms == None and not seed_id == None:
                    seed_cpd = ms.get_seed_compound(seed_id)
                    n['name'] = seed_cpd.name
                #print(n['bigg_id'], seed_id, cmp)
                if not cmp == None:
                    seed_id += '_' + cmp
                cpd_remap[n['bigg_id']] = seed_id

    for rxn_uid in em_merge.escher_graph['reactions']:
        seed_id = None
        cmp = None
        rxn_node = em_merge.escher_graph['reactions'][rxn_uid]
        node_uid = get_midmarker_uid(rxn_node, em)
        if not node_uid == None and node_uid in rxn_uid_to_cluster:
            for maybe_seed_uid in rxn_uid_to_cluster[node_uid]:
                n = em.nodes[maybe_seed_uid]
                if 'rxn_id' in n:
                    seed_id = n['rxn_id']
                    rxn_node['name'] = seed_id
                if 'compartment' in n:
                    cmp = n['compartment']
        if not seed_id == None:
            rxn_remap[rxn_node['bigg_id']] = seed_id
            if not cmp == None:
                rxn_remap[rxn_node['bigg_id']] += '_' + cmp      
                
    em_merge.swap_ids(cpd_remap, rxn_remap)
    return em_merge
    
            
            
    