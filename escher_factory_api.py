import logging
import json
import cobra
import modelseed_escher
import biosapi
import os.path
from modelseed_escher.map.model_merge import validate_map_with_model

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

def match2(to_match, map_stoichiometry, map_compound_remap):
    #to_match = set(map(lambda x : x[0], rxn.cstoichiometry.keys()))
    missing = set()
    mapping = {}
    for cpd_id in map_stoichiometry:
        if cpd_id in map_compound_remap:
            for other_id in map_compound_remap[cpd_id]:
                if other_id in to_match:
                    mapping[cpd_id] = other_id
    missing =  to_match - set(mapping.values())
    return mapping, missing, to_match

def get_stoichiometry(rnode):
    stoichiometry = {}
    for m in rnode['metabolites']:
        stoichiometry[m['bigg_id']] = m['coefficient']
        
    return stoichiometry

def remap_map_reactions(em, bigg_to_seed_rxn, map_compound_remap, exclude, 
                        get_rxn, to_match_func = lambda x : x, id_func = lambda x : x):
    #print('!!!!')
    unmaped_rxn = set()
    map_reaction_remap = {}
    for map_uid in em.escher_graph['reactions']:
        rnode = em.escher_graph['reactions'][map_uid]
        node_id = id_func(rnode['bigg_id'])
        
        if node_id in bigg_to_seed_rxn:
            for db_id in bigg_to_seed_rxn[node_id]:
                map_stoichiometry = get_stoichiometry(rnode)
                #print(db_id)
                #NEED KEGG/METACYC/BIGG provenance
                #print('get_rxn', db_id)
                rxn_cstoich = get_rxn(db_id)
                logger.debug('%s:%s rxn_cstoich: %s', node_id, db_id, rxn_cstoich)
                #print(rxn_cstoich)
                if rxn_cstoich != None:
                    to_match = set(map(lambda x : to_match_func(x[0]), rxn_cstoich.keys()))
                    
                    #print(to_match)
                    #print(map_stoichiometry)
                    mapping, missing, to_match = match2(to_match, map_stoichiometry, map_compound_remap)
                    
                    #print(missing)
                    missing -= exclude
                    logger.debug('%s:%s Match: %s', node_id, db_id, to_match)
                    logger.debug('%s:%s Map S: %s', node_id, db_id, map_stoichiometry)
                    
                    if len(missing) == 0:
                        #print(map_uid, node_id, db_id, map_stoichiometry, rxn.cstoichiometry)
                        #print(mapping, missing, to_match)
                        if not node_id in map_reaction_remap:
                            map_reaction_remap[rnode['bigg_id']] = set()
                        map_reaction_remap[rnode['bigg_id']].add(db_id)
                        for s_id in rnode['segments']:
                            s = rnode['segments'][s_id]
                            #print(s_id, s)
                            break
                    else:
                        logger.debug('%s:%s Miss : %s', node_id, db_id, missing)
                        unmaped_rxn.add(rnode['bigg_id'])
        else:
            logger.debug('%s node_id in bigg_to_seed_rxn = false', node_id)
            unmaped_rxn.add(rnode['bigg_id'])
       
    #print('remap_map_reactions', map_reaction_remap, unmaped_rxn)
    return map_reaction_remap, unmaped_rxn

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
            else:
                node_uid_id[node_uid] = n['bigg_id']
    add_compartment(em, node_uid_cmp)
    for rxn_uid in em.escher_graph['reactions']:
        rnode = em.escher_graph['reactions'][rxn_uid]
        for o in rnode['metabolites']:
            o['bigg_id'] = node_uid_id[o['node_uid']]
        if not 'compartment' in rnode:
            rnode['bigg_id'] += '_' + cmp_id
    return node_uid_cmp

def add_compartment(em, node_uid_cmp):
    for node_uid in node_uid_cmp:
        if node_uid in em.escher_graph['nodes']:
            node = em.escher_graph['nodes'][node_uid]
            node['compartment'] = node_uid_cmp[node_uid]
            
def remap_map_compounds2(em, bigg_to_seed, cmp_mapping, 
                         id_function = lambda x : (x, 0), default_cmp = 'c'):
    """
    hi?
    :param em
    :type em: EscherMap
    :return EscherMap
    """
    map_compound_remap = {}
    unmaped = set()
    node_uid_cmp = {}
    for map_uid in em.escher_graph['nodes']:
        node = em.escher_graph['nodes'][map_uid]
        if node['node_type'] == 'metabolite':
            node_id = node['bigg_id']
            #print(node_id, node)
            bigg_id, cmp = id_function(node_id)
            model_cmp = cmp_mapping[cmp]
            logger.debug("[%s]: %s (%s -> %s)", node_id, bigg_id, cmp, model_cmp)
            #print(node_id, bigg_id, cmp, model_cmp)
            #bigg_id = node_id[:-2]
            #cmp = node_id[-1:]
            
            #print(node_id, bigg_id, cmp)
            model_id = None
            if model_cmp in bigg_to_seed:
                if bigg_id in bigg_to_seed[model_cmp]:
                    node_uid_cmp[map_uid] = cmp
                    model_id = bigg_to_seed[model_cmp][bigg_id]
                    logger.debug("[%s]: %s (%s -> %s)", node_id, model_cmp, bigg_id, bigg_to_seed[model_cmp][bigg_id])
                    
            if model_id == None:
                unmaped.add(node_id)
            else:
                map_compound_remap[node_id] = set(model_id)

    return map_compound_remap, unmaped, node_uid_cmp

def remap_map_reactions2(em, bigg_to_seed_rxn, map_compound_remap, exclude, 
                        get_rxn, to_match_func = lambda x : x, id_func = lambda x : x):
    #print('!!!!')
    unmaped_rxn = set()
    map_reaction_remap = {}
    
    for map_uid in em.escher_graph['reactions']:
        rnode = em.escher_graph['reactions'][map_uid]
        map_node_id = rnode['bigg_id']
        node_id = id_func(rnode['bigg_id'])
        logger.debug("[%s] %s %s", map_node_id, node_id, node_id in bigg_to_seed_rxn)
        if node_id in bigg_to_seed_rxn:
            for db_id in bigg_to_seed_rxn[node_id]:
                map_stoichiometry = get_stoichiometry(rnode)
                #print(db_id)
                #NEED KEGG/METACYC/BIGG provenance
                #print('get_rxn', db_id)
                rxn_cstoich = get_rxn(db_id)
                logger.debug('[%s]:%s rxn_cstoich: %s', map_node_id, db_id, rxn_cstoich)
                #print(rxn_cstoich)
                if rxn_cstoich != None:
                    to_match = set(map(lambda x : to_match_func(x[0]), rxn_cstoich.keys()))
                    
                    #print(to_match)
                    #print(map_stoichiometry)
                    mapping, missing, to_match = match2(to_match, map_stoichiometry, map_compound_remap)
                    
                    #print(missing)
                    missing -= exclude
                    logger.debug('%s:%s Match: %s', node_id, db_id, to_match)
                    logger.debug('%s:%s Map S: %s', node_id, db_id, map_stoichiometry)
                    
                    if len(missing) == 0:
                        #print(map_uid, node_id, db_id, map_stoichiometry, rxn.cstoichiometry)
                        #print(mapping, missing, to_match)
                        if not node_id in map_reaction_remap:
                            map_reaction_remap[rnode['bigg_id']] = set()
                        map_reaction_remap[rnode['bigg_id']].add(db_id)
                        for s_id in rnode['segments']:
                            s = rnode['segments'][s_id]
                            #print(s_id, s)
                            break
                    else:
                        logger.debug('%s:%s Miss : %s', node_id, db_id, missing)
                        unmaped_rxn.add(rnode['bigg_id'])
                else:
                    logger.warning('[%s] unable to get stoichiometry of: %s', map_node_id, db_id)
        else:
            logger.debug('%s node_id in bigg_to_seed_rxn = false', node_id)
            unmaped_rxn.add(rnode['bigg_id'])
       
    #print('remap_map_reactions', map_reaction_remap, unmaped_rxn)
    return map_reaction_remap, unmaped_rxn

            
class EscherFactoryApi:
    
    def __init__(self, escher_manager):
        self.escher_manager = escher_manager
        self.cpd_mapping = {}
        self.rxn_mapping = {}
        self.bios_cache = {}
        self.cpd_match_exclude = set()
        
    def get_model_reaction(self, sbml_id, rxn_id):
        if sbml_id in self.bios_cache:
            all_rxns = self.bios_cache[sbml_id]['rxn']
            m = list(filter(lambda x : x['id'] == rxn_id, all_rxns))
            if len(m) == 1:
                rxn = biosapi.core.model.BiosModelReaction(m[0])
                cstoich = rxn.cstoichiometry
                return cstoich
        return None
    
    def get_cmp_mapping(self, model_id):
        res = {}
        if model_id in self.cmp_mapping:
            res = self.cmp_mapping[model_id]
        return res
    
    def get_cpd_mapping(self, model_id):
        res = {}
        cmp_mapping = self.get_cmp_mapping(model_id)
        for model_cmp in cmp_mapping.values():
            res[model_cmp] = self.get_cpd_mapping_by_compartment(model_id, model_cmp)
        return res
    
    def get_cpd_mapping_by_compartment(self, model_id, model_cmp):
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
    
    def get_rxn_mapping(self, model_id):
        res = {}
        if model_id in self.rxn_mapping:
            mapping = self.rxn_mapping[model_id]
            for sid in mapping:
                for seed_id in mapping[sid]:
                    if not seed_id in res:
                        res[seed_id] = set()
                    res[seed_id].add(sid)
        return res
    
    def lambda_get_model_reaction(self, sbml_id):
        return lambda x : self.get_model_reaction(sbml_id, x)
    
    def translate_to_model2(self, escher_model_id, escher_map_id, sbml_id, target_cmp, cpd_mapping, rxn_mapping, cmp_mapping = {}):
        em = self.escher_manager.get_map('ModelSEED', escher_model_id, escher_map_id)
        em.add_uid_to_reaction_metabolites()
        move_to_compartment(target_cmp, em)

        map_compound_remap, unmaped_cpd, node_uid_cmp = remap_map_compounds2(
            em, cpd_mapping, cmp_mapping, lambda x : x.split('_'))
        map_reaction_remap, unmaped_rxn = remap_map_reactions2(
            em, rxn_mapping, map_compound_remap, 
            self.cpd_match_exclude, 
            self.lambda_get_model_reaction(sbml_id),
            lambda x : x, lambda x : x[:-2])

        logger.warning('unmaped_cpd: %d', len(unmaped_cpd))
        logger.warning('unmaped_rxn: %d', len(unmaped_rxn))

        #print(sbml_id, map_reaction_remap)
        
        cpd_remap = {}
        rxn_remap = {}
        for k in map_compound_remap:
            cpd_remap[k] = list(map_compound_remap[k])[0] + '@' + sbml_id
        for k in map_reaction_remap:
            rxn_remap[k] = list(map_reaction_remap[k])[0] + '@' + sbml_id
        #print(sbml_id, rxn_remap)
        em.swap_ids(cpd_remap, rxn_remap)
        em.delete_reactions(unmaped_rxn)
        em.delete_metabolites(unmaped_cpd)
        return em
    
    def translate_to_model(self, escher_model_id, escher_map_id, sbml_id, target_cmp, cpd_mapping, rxn_mapping):
        em = self.escher_manager.get_map('ModelSEED', escher_model_id, escher_map_id)
        em.add_uid_to_reaction_metabolites()

        move_to_compartment(target_cmp, em)

        map_compound_remap, unmaped_cpd, node_uid_cmp = remap_map_compounds(em, 
                                                                            cpd_mapping, lambda x : x.split('_'))

        logger.warning('unmaped_cpd: %d', len(unmaped_cpd))

        map_reaction_remap, unmaped_rxn = remap_map_reactions(em, rxn_mapping, map_compound_remap, self.cpd_match_exclude,
                                                              self.lambda_get_model_reaction(sbml_id),
                                                              lambda x : x, lambda x : x[:-2])

        logger.warning('unmaped_rxn: %d', len(unmaped_rxn))

        print(sbml_id, map_reaction_remap)
        
        cpd_remap = {}
        rxn_remap = {}
        for k in map_compound_remap:
            cpd_remap[k] = list(map_compound_remap[k])[0] + '@' + sbml_id
        for k in map_reaction_remap:
            rxn_remap[k] = list(map_reaction_remap[k])[0] + '@' + sbml_id
        #print(sbml_id, rxn_remap)
        em.swap_ids(cpd_remap, rxn_remap)
        em.delete_reactions(unmaped_rxn)
        em.delete_metabolites(unmaped_cpd)
        return em
    
    def build_me_a_map(self, escher_model_id, escher_map_id, sbml_id, target_cmp, cpd_mapping):
        em = self.escher_manager.get_map('ModelSEED', escher_model_id, escher_map_id)
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
    
    
    def get_model(self, model_id):
        model = None
        model_path = self.model_path + '/TempModels/' + model_id + '.json'
        if os.path.exists(model_path):
            with open(model_path, 'r') as fh:
                model = cobra.io.from_json(fh.read())
        return model
    
    def build_grid(self, map_assembly, grid_setup):
        em_list = []
        for ma in map_assembly:
            escher_model_id = 'ModelSEED'
            escher_map_id = ma['map_id']
            sbml_id = ma['sbml_id']
            target_cmp = ma['cmp_target'][0]
            cmp_sbml = ma['cmp_sbml'][0]
            
            cmp_mapping = self.get_cmp_mapping(sbml_id)
            cpd_mapping = self.get_cpd_mapping(sbml_id)
            rxn_mapping = self.get_rxn_mapping(sbml_id)
            
            print(escher_map_id)
            #em_ = self.build_me_a_map(escher_model_id, escher_map_id, 
            #                          sbml_id, target_cmp, cpd_mapping)
            
            logger.warning('translate_to_model: %s [%s > %s]', sbml_id, cmp_sbml, target_cmp)
            em_ = self.translate_to_model2(escher_model_id, 
                                                    escher_map_id, 
                                                    sbml_id,
                                                    target_cmp,
                                                    cpd_mapping,
                                                    rxn_mapping,
                                          cmp_mapping)
            model = self.get_model(sbml_id)
            if not model == None:
                validate_map_with_model(em_, model)
            else:
                logger.warning('unable to get model: %s', sbml_id)
            
            
            em_list.append(em_)
        builder = modelseed_escher.EscherGrid()
        master = builder.build(em_list, grid_setup)
        
        return master