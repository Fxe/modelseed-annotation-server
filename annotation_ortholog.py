import logging
import json
import copy
import cobrakbase
import os
from cobrakbase.core import KBaseFBAModel
from cobra.core.gene import Gene, ast2str, eval_gpr, parse_gpr

logger = logging.getLogger(__name__)

ACTIVE_GENOMES = [
    'GCF_000209165.1',
 'GCF_000182925.2',
 'GCF_000002525.2',
 'GCF_000184455.2',
 'GCF_000149615.1',
 'GCF_000226395.1',
 'GCF_000002515.2',
 'GCF_000027005.1',
 'GCF_000006335.2',
 'GCF_000146045.2',
 'GCF_000091025.4',
 'GCF_000002545.3',
 'Mucor_circinelloides_CBS277.49_v2.0',
 'GCF_000002855.3'
                 ]

def read_json(filename):
    data = None
    with open(filename, 'r') as f:
        data = json.loads(f.read())
    return data

def build_annotation_ortholog(kbase, path_to_cache, bios):
    #ref = '49164/21/1' #FungalTemplate13GenomeOrthologs-OrthoMCL janakakbase:narrative_1570052138482
    #env_ws = 'janakakbase:narrative_1570052138482'
    #ortho = kbase.get_object('FungalOrthos', env_ws)
    #info = kbase.get_object_info_from_ref(ref)
    #ortho = kbase.get_object(info.id, info.workspace_id)
    
    env_ws = 'filipeliu:narrative_1575329472745'
    ortho = kbase.get_object('FungalTemplateAndRepresentativeSetOrthosV0', env_ws)
    
    ref_to_genome = {}
    genome_id_to_ref = {}

    annotation_path = path_to_cache + '/cache/genomes/'

    for genome_ref in ortho['genome_refs']:
        info = kbase.get_object_info_from_ref(genome_ref)
        genome_data = None
        if os.path.exists(annotation_path + info.id + '.json'):
            print('load local', genome_ref)
            genome_data = read_json(annotation_path + info.id + '.json')
        else:
            print('load kbase', genome_ref)
            genome_data = kbase.get_object(info.id, info.workspace_id)
        genome = cobrakbase.core.KBaseGenome(genome_data)
        print(info.id, info.workspace_id)
        ref_to_genome[genome_ref] = genome

    for ref in ref_to_genome:
        genome = ref_to_genome[ref]
        genome_id_to_ref[genome.id] = ref
        print(ref, genome.id, genome.data['scientific_name'])
        
    ws_fungi = 'jplfaria:narrative_1510597445008'
    model_ids = [
        'iNL895_KBase2',
        'iCT646_KBase2',

        'iMM904_KBase3',
        'iTO977_KBase2',
        'iSS884_KBase3',
        'iLC915_KBase3',

        'iWV1213_KBase2',
        'iAL1006_KBase3',

        'iRL766_KBase2',
        'iMA871_KBase3',

        'iJDZ836_KBase3',
        'iWV1314_KBase3',
        'iOD907_KBase2',
        'iJL1454_KBase2',
        'iNX804_KBase2',
        'yeast_6.06_KBase2',
        'yeast_7.6_KBase2',
    ]
    kbase_models = {}
    for model_id in model_ids:
        logger.warning('load [%s:%s]', model_id, ws_fungi)
        fbamodel = kbase.get_from_ws(model_id, ws_fungi)
        kbase_models[fbamodel.id] = fbamodel
        
    annotation_orth = AnnotationOrtholog(ortho, bios)
    annotation_orth.ref_to_genome = ref_to_genome
    annotation_orth.genome_id_to_ref = genome_id_to_ref
    
    for model_id in kbase_models:
        fbamodel = kbase_models[model_id]
        genome_ref = fbamodel.data['genome_ref']
        info = kbase.get_object_info_from_ref(genome_ref)

        annotation_orth.models[fbamodel.id] = fbamodel
        annotation_orth.model_to_genome[fbamodel.id] = info.id
    annotation_orth.index_model_gene_reactions()
    
    return annotation_orth

class AnnotationOrtholog:
    
    def __init__(self, ortho, bios):
        self.model_rxn_gpr = {}
        self.model_rxn_mapping = {}
        self.model_cpd_mapping = {}
        self.bios = bios
        self.models = {}
        self.model_to_genome = {}
        self.ortho = ortho
        self.ref_to_genome = {}
        self.genome_id_to_ref = {}
        self.genome_ref_to_gene_to_ortholog_index = {}
        self.model_gene_reaction = {}
        
        for ortholog_index in range(len(ortho['orthologs'])):
            ortholog = ortho['orthologs'][ortholog_index]
            ortholog_id = ortholog['id']
            for o in ortholog['orthologs']:
                genome_ref = o[2]
                gene_id = o[0]
                if not genome_ref in self.genome_ref_to_gene_to_ortholog_index:
                    self.genome_ref_to_gene_to_ortholog_index[genome_ref] = {}
                self.genome_ref_to_gene_to_ortholog_index[genome_ref][gene_id] = ortholog_index
        pass
    
    def get_genes(self, gpr_exp):
        gpr_exp = gpr_exp.replace(' AND ', ' and ')
        gpr_exp = gpr_exp.replace(' OR ', ' or ')
        #gpr_exp = gpr_exp.replace(' and', ' and ')
        #gpr_exp = gpr_exp.replace('; ', ' or ')
        genes = set()
        try:
            exp, genes = parse_gpr(gpr_exp)
            return genes
        except:
            return None

    @staticmethod
    def get_original_id(rxn):
        sid = None
        if 'string_attributes' in dir(rxn) and 'original_id' in rxn.string_attributes:
            return rxn.string_attributes['original_id']
        return sid

    def index_model_gene_reactions(self):
        for model_id in self.models:
            self.model_gene_reaction[model_id] = {}
            for rxn in self.models[model_id].reactions:
                if 'imported_gpr' in dir(rxn):
                    gpr_exp = rxn.imported_gpr
                    if model_id == 'iMA871':
                        gpr_exp = gpr_exp.replace('; ', ' or ')
                    if model_id == 'iCT646' and self.get_original_id(rxn) == 'R_FAS100':
                        gpr_exp = 'CTRG_02936 OR (CTRG_05241 AND CTRG_02501)'
                    genes = self.get_genes(gpr_exp)
                    if genes is None:
                        logger.warning('[%s] %s: %s', model_id, self.get_original_id(rxn), rxn.imported_gpr)
                    elif len(genes) > 0:
                        for g in genes:
                            if g not in self.model_gene_reaction[model_id]:
                                self.model_gene_reaction[model_id][g] = set()
                            self.model_gene_reaction[model_id][g].add(self.get_original_id(rxn))

    def get_model_reaction_from_seed_id(self, model_id, seed_rxn_id, score,
                                        compartment=None, standard_compartment=True):
        res = {}
        model_reactions_by_database_id = self.bios.get_model_reactions_by_database_id(model_id,
                                                                                      seed_rxn_id,
                                                                                      'ModelSeedReaction')

        for p in model_reactions_by_database_id:
            author_str = p[1]
            author_score = dict(map(lambda x: x.split(':'), author_str.split(';')))
            for a in author_score:
                if int(author_score[a]) >= score:
                    res[p[0].id] = p[0]
        return res

    def get_ortholog(self, genome_id, gene_id):
        if not genome_id in self.genome_id_to_ref:
            return None
        
        genome_ref = self.genome_id_to_ref[genome_id]
        if genome_ref in self.genome_ref_to_gene_to_ortholog_index and gene_id in self.genome_ref_to_gene_to_ortholog_index[genome_ref]:
            ortholog_index = self.genome_ref_to_gene_to_ortholog_index[genome_ref][gene_id]
            return self.ortho['orthologs'][ortholog_index]
        
        return None
    
    def get_orthologs_from_seed_rxn_id(self, seed_rxn_id):
        model_rxn_grp = {}
        genome_match = {}
        for model_id in self.models:

            model_rxn_grp[model_id] = {}
            print(model_id)
            for r in self.models[model_id].reactions:
                if 'ModelSeedReaction' in r.data['dblinks'] and seed_rxn_id in r.data['dblinks']['ModelSeedReaction']:
                    model_rxn_grp[model_id][r.id] = r

        all_orthologs = set()
        matched_orthologs = set()
        for model_id in model_rxn_grp:
            genome_id = self.model_to_genome[model_id]
            if not genome_id in genome_match:
                genome_match[genome_id] = {}
            genome_match[genome_id][model_id] = []
            genome = self.ref_to_genome[self.genome_id_to_ref[genome_id]]
            for rxn_id in model_rxn_grp[model_id]:
                gpr_exp = model_rxn_grp[model_id][rxn_id].data['imported_gpr']
                genes = self.get_genes(gpr_exp)
                print(genome_id, gpr_exp, genes)

                exp_match = {
                    'rxn_id': rxn_id,
                    'gpr': model_rxn_grp[model_id][rxn_id].data['imported_gpr'],
                    'genes': {}
                }

                for gene_id in genes:
                    features = list(filter(lambda x : x['id'] == gene_id, genome.features))
                    if len(features) == 1:
                        ortholog = self.get_ortholog(genome_id, gene_id)
                        matched_orthologs.add((genome_id, gene_id, ortholog['id']))
                        for o in ortholog['orthologs']:
                            o_genome_id = self.ref_to_genome[o[2]].id
                            if o_genome_id in self.model_to_genome.values():
                                all_orthologs.add((o_genome_id, o[0], ortholog['id']))
                        print(features[0]['id'], features[0]['function'], ortholog['id'])
                        exp_match['genes'][gene_id] = {
                            'ortholog_id' : ortholog['id'],
                            'function' : features[0]['function']
                        }

                print(exp_match)
                genome_match[genome_id][model_id].append(exp_match)

        #print(all_orthologs - matched_orthologs)
        
        return all_orthologs, matched_orthologs, genome_match
    
    def bios_get_model_reaction(self, model_id, rxn_id):
        return self.bios.get_model_reaction(model_id, rxn_id)

    def bios_get_model_species(self, model_id, spi_id):
        return self.bios.get_model_specie(model_id, spi_id)

    def get_rxn_compartment(self, rxn, model_id):
        cmp = set()
        for o in rxn.cstoichiometry:
            spi = self.bios_get_model_species(model_id, o[0])
            if spi == None:
                return None
            cmp.add(spi['compartment'])
        return list(sorted(cmp))
    
    def get_kbase_reaction(self, model_id, sid):
        l = list(filter(lambda x : self.get_original_id(x) == sid, self.models[model_id].reactions))
        if len(l) == 1:
            return l[0]
        return None
    
    def get_model_rxn_grp(self, seed_rxn_id, compartment=None, standard_compartment=True):
        model_rxn_grp = {}
        for model_id in self.models:
            model_rxn_grp[model_id] = {}
            if model_id in self.model_rxn_mapping:
                for rxn_id in self.model_rxn_mapping[model_id]:
                    if len(self.model_rxn_mapping[model_id][rxn_id]) > 0:
                        rxn = self.bios_get_model_reaction(model_id, rxn_id)
                        cmp = self.get_rxn_compartment(rxn, model_id)
                        seed_id = self.get_reaction_mapping(model_id, rxn.id)
                        if seed_rxn_id in seed_id:
                            if compartment == None:
                                kbase_rxn = self.get_kbase_reaction(model_id, rxn.id)
                                if not kbase_rxn == None:
                                    model_rxn_grp[model_id][rxn.id] = kbase_rxn
                                else:
                                    print('!!!!')
                            else:
                                model_rxn_grp[model_id][rxn.id] = rxn
                    #print(rxn_id, get_rxn_compartment(rxn, model_id))
        return model_rxn_grp
    
    def get_reaction_mapping(self, model_id, rxn_id):
        if model_id in self.model_rxn_mapping and rxn_id in self.model_rxn_mapping[model_id]:
            return self.model_rxn_mapping[model_id][rxn_id]
        return None

    def get_model_rxn_grp2(self, seed_rxn_id, compartment=None, standard_compartment=True, score=5):
        model_rxn_grp = {}
        for model_id in self.models:
            print(model_id)
            # model_rxn_grp[model_id] = get_model_rxn_grp3(aa, model_id, seed_rxn_id, compartment, standard_compartment)
            model_rxn_grp[model_id] = self.get_model_reaction_from_seed_id(model_id, seed_rxn_id, score, compartment, standard_compartment)
        return model_rxn_grp

    def get_gpr(self, model_id, mrxn):
        if model_id in self.model_rxn_gpr and mrxn.id in self.model_rxn_gpr[model_id]:
            return self.model_rxn_gpr[model_id][mrxn.id]
        return ''

    def get_orthologs_from_seed_rxn_id3(self, seed_rxn_id, compartment=None, standard_compartment=True):
        model_rxn_grp = self.get_model_rxn_grp2(seed_rxn_id, compartment, standard_compartment, 5)
        #print(model_rxn_grp)
        genome_match = {}
        all_orthologs = set()
        matched_orthologs = set()
        for model_id in model_rxn_grp:
            genome_id = self.model_to_genome[model_id]
            if genome_id in self.genome_id_to_ref:
                if genome_id not in genome_match:
                    genome_match[genome_id] = {}
                genome_match[genome_id][model_id] = []
                genome = self.ref_to_genome[self.genome_id_to_ref[genome_id]]
                for rxn_id in model_rxn_grp[model_id]:
                    model_rxn = model_rxn_grp[model_id][rxn_id]
                    gpr_exp = self.get_gpr(model_id, model_rxn)

                    #print(model_id, rxn_id, gpr_exp)

                    genes = self.get_genes(gpr_exp)

                    logger.debug("%s %s %s", genome_id, gpr_exp, genes)

                    exp_match = {
                        'rxn_id': rxn_id,
                        'gpr': gpr_exp,
                        'genes': {}
                    }

                    for gene_id in genes:
                        features = list(filter(lambda x: x['id'] == gene_id, genome.features))
                        if len(features) == 1:
                            ortholog = self.get_ortholog(genome_id, gene_id)
                            matched_orthologs.add((genome_id, gene_id, ortholog['id']))
                            for o in ortholog['orthologs']:
                                o_genome_id = self.ref_to_genome[o[2]].id
                                if o_genome_id in self.model_to_genome.values():
                                    all_orthologs.add((o_genome_id, o[0], ortholog['id']))
                            gene_function = '?' if 'function' not in features[0] else features[0]['function']
                            #print(features[0]['id'], gene_function, ortholog['id'])
                            exp_match['genes'][gene_id] = {
                                'ortholog_id' : ortholog['id'],
                                'function' : gene_function
                            }
                        else:
                            print('!', gene_id, len(features))
                    logger.debug("%s", exp_match)
                    genome_match[genome_id][model_id].append(exp_match)

        return all_orthologs, matched_orthologs, genome_match

    def get_orthologs_from_seed_rxn_id2(self, seed_rxn_id, compartment = None, standard_compartment = True):
        model_rxn_grp = self.get_model_rxn_grp(seed_rxn_id, compartment, standard_compartment)
        genome_match = {}

        all_orthologs = set()
        matched_orthologs = set()
        for model_id in model_rxn_grp:
            genome_id = self.model_to_genome[model_id]
            if genome_id in self.genome_id_to_ref:
                if not genome_id in genome_match:
                    genome_match[genome_id] = {}
                genome_match[genome_id][model_id] = []
                genome = self.ref_to_genome[self.genome_id_to_ref[genome_id]]
                for rxn_id in model_rxn_grp[model_id]:
                    gpr_exp = model_rxn_grp[model_id][rxn_id].data['imported_gpr']
                    genes = self.get_genes(gpr_exp)

                    logger.debug("%s %s %s", genome_id, gpr_exp, genes)
                    #print(genome_id, gpr_exp, genes)
                    exp_match = {
                        'rxn_id' : rxn_id,
                        'gpr' : model_rxn_grp[model_id][rxn_id].data['imported_gpr'],
                        'genes' : {}
                    }

                    for gene_id in genes:
                        features = list(filter(lambda x : x['id'] == gene_id, genome.features))
                        if len(features) == 1:
                            ortholog = self.get_ortholog(genome_id, gene_id)
                            matched_orthologs.add((genome_id, gene_id, ortholog['id']))
                            for o in ortholog['orthologs']:
                                o_genome_id = self.ref_to_genome[o[2]].id
                                if o_genome_id in self.model_to_genome.values():
                                    all_orthologs.add((o_genome_id, o[0], ortholog['id']))
                            gene_function = '?' if not 'function' in features[0] else features[0]['function']
                            print(features[0]['id'], gene_function, ortholog['id'])
                            exp_match['genes'][gene_id] = {
                                'ortholog_id' : ortholog['id'],
                                'function' : gene_function
                            }
                        else:
                            print('!', gene_id, len(features))
                    logger.debug("%s", exp_match)
                    genome_match[genome_id][model_id].append(exp_match)
                    #print(rxn_id, model_rxn_grp[model_id][rxn_id], gpr_exp)
        return all_orthologs, matched_orthologs, genome_match
    
    def process_data(self, all_orthologs, genome_match):

        model_data = {
            'orthologs' : {},
            'genomes' : {},
        }
        for o in all_orthologs:
            genome_id = o[0]
            gene_id = o[1]
            ortholog_id = o[2]

            if not ortholog_id in model_data['orthologs']:
                model_data['orthologs'][ortholog_id] = {}
            if not genome_id in model_data['orthologs'][ortholog_id]:
                model_data['orthologs'][ortholog_id][genome_id] = {}
            if not gene_id in model_data['orthologs'][ortholog_id][genome_id]:
                model_data['orthologs'][ortholog_id][genome_id][gene_id] = '?'

            genome = self.ref_to_genome[self.genome_id_to_ref[genome_id]]
            features = list(filter(lambda x : x['id'] == gene_id, genome.features))
            if len(features) == 1:
                model_data['orthologs'][ortholog_id][genome_id][gene_id] = features[0]['function']
            #model_data['orthologs'][ortholog_id].append(list(o))
        for genome_id in genome_match:


            genome = self.ref_to_genome[self.genome_id_to_ref[genome_id]]
            #print(ref, genome.id, genome.data['scientific_name'])

            model_data['genomes'][genome_id] = {
                'name' : genome.data['scientific_name'],
                'source' : {}
            }
            for source_id in genome_match[genome_id]:
                model_data['genomes'][genome_id]['source'][source_id] = genome_match[genome_id][source_id]
        return model_data
    
    
    def process_data2(self, all_orthologs, genome_match):

        model_data = {
            'gene_reaction' : {},
            'orthologs' : {},
            'genomes' : {},
        }
        for o in all_orthologs:
            genome_id = o[0]
            gene_id = o[1]
            ortholog_id = o[2]

            if not ortholog_id in model_data['orthologs']:
                model_data['orthologs'][ortholog_id] = {}
                ortholog_data = list(filter(lambda x : x['id'] == ortholog_id, self.ortho['orthologs']))[0]
                for o in ortholog_data['orthologs']:
                    genome_ref = o[2]
                    genome = self.ref_to_genome[genome_ref]
                    gene_id = o[0]
                    if not genome.id in model_data['orthologs'][ortholog_id]:
                        model_data['orthologs'][ortholog_id][genome.id] = {}
                    if not gene_id in model_data['orthologs'][ortholog_id][genome.id]:
                        model_data['orthologs'][ortholog_id][genome.id][gene_id] = '?'

                    #genome = aaaaaaaaa.ref_to_genome[aaaaaaaaa.genome_id_to_ref[genome_id]]
                    features = list(filter(lambda x : x['id'] == gene_id, genome.features))
                    if len(features) == 1:
                        gene_function = '?' if not 'function' in features[0] else features[0]['function']
                        model_data['orthologs'][ortholog_id][genome.id][gene_id] = gene_function
                    else:
                        print(len(features))
            #model_data['orthologs'][ortholog_id].append(list(o))
        for genome_id in genome_match:


            genome = self.ref_to_genome[self.genome_id_to_ref[genome_id]]
            #print(ref, genome.id, genome.data['scientific_name'])

            model_data['genomes'][genome_id] = {
                'name' : genome.data['scientific_name'],
                'source' : {}
            }
            for source_id in genome_match[genome_id]:
                model_data['genomes'][genome_id]['source'][source_id] = genome_match[genome_id][source_id]

        for o in model_data['orthologs']:
            for genome_id in model_data['orthologs'][o]:
                if not genome_id in model_data['gene_reaction']:
                    model_data['gene_reaction'][genome_id] = {}
                for gene_id in model_data['orthologs'][o][genome_id]:
                    if not gene_id in model_data['gene_reaction'][genome_id]:
                        model_data['gene_reaction'][genome_id][gene_id] = []

        for model_id in self.model_to_genome:
            genome_id = self.model_to_genome[model_id]
            if genome_id in model_data['gene_reaction']:
                for gene_id in model_data['gene_reaction'][genome_id]:
                    if model_id in self.model_gene_reaction and gene_id in self.model_gene_reaction[model_id]:
                        model_data['gene_reaction'][genome_id][gene_id] = list(self.model_gene_reaction[model_id][gene_id])

        return model_data