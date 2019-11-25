import json
import copy
import cobrakbase
from cobra.core.gene import Gene, ast2str, eval_gpr, parse_gpr

def read_json(filename):
    data = None
    with open(filename, 'r') as f:
        data = json.loads(f.read())
    return data

def build_annotation_ortholog(kbase):
    ref = '49164/21/1' #FungalTemplate13GenomeOrthologs-OrthoMCL janakakbase:narrative_1570052138482
    info = kbase.get_object_info_from_ref(ref)
    ortho = kbase.get_object(info.id, info.workspace_id)
    
    ref_to_genome = {}
    genome_id_to_ref = {}

    annotation_path = '/Users/fliu/workspace/jupyter/data/www/fungi_viz/data/annotation_cache/'

    for genome_ref in ortho['genome_refs']:
        info = kbase.get_object_info_from_ref(genome_ref)
        genome_data = read_json(annotation_path + info.id + '.json')
        #genome_data = kbase.get_object(info.id, info.workspace_id)
        genome = cobrakbase.core.KBaseGenome(genome_data)
        print(info.id, info.workspace_id)
        ref_to_genome[genome_ref] = genome

    for ref in ref_to_genome:
        genome = ref_to_genome[ref]
        genome_id_to_ref[genome.id] = ref
        print(ref, genome.id, genome.data['scientific_name'])
        
    model_ids = [
        'iMM904',
        'iJDZ836',
        'iAL1006',
        #'yeast_7.6_KBase2.json',
        #'yeast_6.06_KBase2.json',
        #'iRL766_KBase2.json',
        #'iSS884_KBase2.json',
        #'iNL895_KBase2.json',
        #'iCT646_KBase2.json',
        #'iNX804_KBase2.json',
        #'iMA871_KBase2.json',
        #'iOD907_KBase2.json',
        #'iWV1314_KBase2.json',
        #'iTO977_KBase2.json',
        #'iWV1213_KBase2.json',
        #'iJL1454_KBase2.json',
        #'iLC915_KBase2.json',
        #'iJDZ836_KBase2.json',
        #'iAL1006_KBase2.json'
    ]
    
    gpr_prefix = '_KBase3.json'
    base_prefix = '_Base.json'
    
    models = {}
    model_to_genome = {}
    
    for model_id in model_ids:
        fbamodel_base = cobrakbase.core.model.KBaseFBAModel(read_json('/Users/fliu/workspace/jupyter/data/www/fungi_viz/data/kbase/' + model_id + base_prefix))
        fbamodel_gpr = cobrakbase.core.model.KBaseFBAModel(read_json('/Users/fliu/workspace/jupyter/data/www/fungi_viz/data/kbase/' + model_id + gpr_prefix))
        genome_ref = fbamodel_gpr.data['genome_ref']
        info = kbase.get_object_info_from_ref(genome_ref)
        model_to_genome[model_id] = info.id
        #print(info.id, info.workspace_id)
        oid_to_rxn = {}
        for rxn in fbamodel_gpr.reactions:
            oid = rxn.get_original_id()
            if not oid in oid_to_rxn:
                oid_to_rxn[oid] = rxn
            else:
                logger.warning('! %s %s', model_id, oid)

        deleted = set()
        for rxn in fbamodel_base.reactions:
            if rxn.id in oid_to_rxn:
                rxn.data['modelReactionProteins'] = copy.deepcopy(oid_to_rxn[rxn.id].data['modelReactionProteins'])
                #print(rxn, rxn.data['dblinks'])
                #print(oid_to_rxn[rxn.id])
            else:
                deleted.add(rxn.id)

        fbamodel_base.delete_reactions(deleted)
            #fbamodel = cobrakbase.core.KBaseFBAModel(read_json('./fungi_viz/data/kbase/iMM904_KBase2.json'))
        #cpd_db_links, rxn_db_links = load_mapping_helper(fbamodel_base.id)
        manual_comp = {}
        #ff = FungiFixer({}, {}, manual_comp)
            #ff.fix_mapping = fix_mapping
            #ff.cpd_db_links = cpd_db_links
            #ff.rxn_db_links = rxn_db_links
            #ff.merge_compounds = merges
            #ff.delete_compounds
            #ff.delete_comp = delete_comp
            #ff.manual_fix = manual_fix
        #ff.fix(fbamodel_base)
        models[model_id] = fbamodel_base
        fbamodel_base.update_indexes()
        for cmp in fbamodel_base.data['modelcompartments']:
            if 'z' in cmp['id']:
                print(model_id, cmp['id'], cmp['label'])
            #    pass
        
    annotation_orth = AnnotationOrtholog(ortho)
    annotation_orth.ref_to_genome = ref_to_genome
    annotation_orth.genome_id_to_ref = genome_id_to_ref
    annotation_orth.models = models
    annotation_orth.model_to_genome = model_to_genome
    
    
    return annotation_orth

class AnnotationOrtholog:
    
    def __init__(self, ortho):
        self.models = {}
        self.model_to_genome = {}
        self.ortho = ortho
        self.ref_to_genome = {}
        self.genome_id_to_ref = {}
        self.genome_ref_to_gene_to_ortholog_index = {}
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
    
    def get_ortholog(self, genome_id, gene_id):
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
                gpr_exp = gpr_exp.replace('AND', 'and')
                gpr_exp = gpr_exp.replace('OR', 'or')
                exp, genes = parse_gpr(gpr_exp)
                print(genome_id, gpr_exp, genes)

                exp_match = {
                    'gpr' : gpr_exp,
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
                        print(features[0]['id'], features[0]['function'], ortholog['id'])
                        exp_match['genes'][gene_id] = {
                            'ortholog_id' : ortholog['id'],
                            'function' : features[0]['function']
                        }

                print(exp_match)
                genome_match[genome_id][model_id].append(exp_match)

        #print(all_orthologs - matched_orthologs)
        
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