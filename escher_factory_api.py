import logging
import json
import modelseed_escher

logger = logging.getLogger(__name__)

def process_build_data_input(build_data):
    map_assembly = []
    for assembly_content in build_data:
        sbml_id, cmp_sbml, cmp_target, map_id = assembly_content.split(';')
        if (len(sbml_id) > 0 and len(cmp_sbml) > 0 and len(cmp_target) > 0 and len(map_id) > 0):
            #print(sbml_id, cmp_sbml, cmp_target, map_id)
            _, map_id = map_id.split('.', 1)
            map_assembly.append({
                'map_id' : map_id,
                'cmp_target' : cmp_target.split(':'),
                'cmp_sbml' : cmp_sbml.split(':'),
                'sbml_id' : sbml_id,
            })
        else:
            logger.warning("discard: %s", assembly_content)
    return map_assembly

def remap_map_compounds(em, bigg_to_seed, id_function = lambda x : (x, 0)):
    map_compound_remap = {}
    unmaped = set()
    node_uid_cmp = {}
    for map_uid in em.escher_graph['nodes']:
        node = em.escher_graph['nodes'][map_uid]
        if node['node_type'] == 'metabolite':
            node_id = node['bigg_id']
            #print(node_id, node)
            bigg_id, cmp = id_function(node_id)
            #print(node_id, bigg_id, cmp)
            #bigg_id = node_id[:-2]
            #cmp = node_id[-1:]
            
            #print(node_id, bigg_id, cmp)
            if bigg_id in bigg_to_seed:
                node_uid_cmp[map_uid] = cmp
                map_compound_remap[node_id] = set(bigg_to_seed[bigg_id])
                #print(map_uid, bigg_id, cmp, bigg_to_seed[bigg_id])
            else:
                unmaped.add(node_id)
    return map_compound_remap, unmaped, node_uid_cmp

def get_cmps(rxn, uid_to_spi):
    cmps = set()
    for t in rxn['bios_stoichiometry']['r']:
        if not t[1] in uid_to_spi:
            return None
        s = uid_to_spi[t[1]]
        cmps.add(s['compartment'])
        #print(t, s['compartment'], s['id'])
    for t in rxn['bios_stoichiometry']['l']:
        if not t[1] in uid_to_spi:
            return None
        s = uid_to_spi[t[1]]
        cmps.add(s['compartment'])
        #print(t, s['compartment'], s['id'])
    return cmps

def load_model(model_id, rxns, spis, model_to_compatment_to_rxn, model_to_compatment_to_spi):
    if not model_id in model_to_compatment_to_rxn:
        model_to_compatment_to_rxn[model_id] = {}
    if not model_id in model_to_compatment_to_spi:
        model_to_compatment_to_spi[model_id] = {}
    
    uid_to_spi = {}
    for o in spis:
        uid_to_spi[o['bios_id']] = o
        if 'compartment' in o:
            if not o['compartment'] in model_to_compatment_to_spi[model_id]:
                model_to_compatment_to_spi[model_id][o['compartment']] = {}
            model_to_compatment_to_spi[model_id][o['compartment']][o['id']] = o
    for o in rxns:
        cmps = get_cmps(o, uid_to_spi)
        if not cmps == None:
            cmp_key = tuple(sorted(cmps))
            if not cmp_key in model_to_compatment_to_rxn[model_id]:
                model_to_compatment_to_rxn[model_id][cmp_key] = {}
            model_to_compatment_to_rxn[model_id][cmp_key][o['id']] = o

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

def add_compartment(em, node_uid_cmp):
    for node_uid in node_uid_cmp:
        if node_uid in em.escher_graph['nodes']:
            node = em.escher_graph['nodes'][node_uid]
            node['compartment'] = node_uid_cmp[node_uid]

class EscherFactoryApi:
    
    def __init__(self, escher_manager):
        self.escher_manager = escher_manager
        self.cpd_mapping = {}
    
    def get_cpd_mapping(self, model_id, model_cmp):
        res = {}
        if model_id in self.cpd_mapping and model_cmp in self.cpd_mapping[model_id]:
            mapping = self.cpd_mapping[model_id][model_cmp]
            for sid in mapping:
                for seed_id in mapping[sid]:
                    if not seed_id in res:
                        res[seed_id] = [sid]
                    elif not sid in res[seed_id]:
                        logger.warning('%s: [%s] ! [%s]', seed_id, sid, res[seed_id])
        return res
    
    def build_me_a_map(self, escher_model_id, escher_map_id, sbml_id, target_cmp, cpd_mapping):
        escher_map = self.escher_manager.get_map('ModelSEED', escher_model_id, escher_map_id)
        em = modelseed_escher.core.EscherMap(json.loads(escher_map))
        em.add_uid_to_reaction_metabolites()

        test = em.clone()
        move_to_compartment(target_cmp, test)
        test2 = test.clone()
        #TODO: load mapping from somewhere
        #cpd_mapping = {
        #    'cpd00001' : ['M_h2o_c'],
        #    'cpd00009' : ['M_pi_c']
        #}
        map_compound_remap, unmaped_cpd, node_uid_cmp = remap_map_compounds(test2, 
                                                                            cpd_mapping, lambda x : x.split('_'))
        cpd_remap = {}
        rxn_remap = {}
        for k in map_compound_remap:
            cpd_remap[k] = list(map_compound_remap[k])[0] + '@' + sbml_id
        #for k in map_reaction_remap:
        #    rxn_remap[k] = list(map_reaction_remap[k])[0]
        test.swap_ids(cpd_remap, rxn_remap)
        return test
    
    
    def build_grid(self, map_assembly, grid_setup):
        em_list = []
        for ma in map_assembly:
            escher_model_id = 'ModelSEED'
            escher_map_id = ma['map_id']
            sbml_id = ma['sbml_id']
            target_cmp = ma['cmp_target'][0]
            cmp_sbml = ma['cmp_sbml'][0]
            
            cpd_mapping = self.get_cpd_mapping(sbml_id, cmp_sbml)
            print(escher_map_id)
            em_ = self.build_me_a_map(escher_model_id, escher_map_id, 
                                      sbml_id, target_cmp, cpd_mapping)
            em_list.append(em_)
        builder = modelseed_escher.EscherGrid()
        master = builder.build(em_list, grid_setup)
        
        return master