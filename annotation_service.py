class AnnotationService:
    
    def __init__(self, annotation_api, curation_api):
        self.annotation_api = annotation_api
        self.curation_api = curation_api
        
    def get_tooltip_annotation_data(self, rxn_id, template_id, genome_set_id):
        #get manual assigned Functions
        rxn_annotation_manual_function = self.curation_api.get_manual_function(rxn_id, template_id)
        #get manual assigned KEGG KO's
        rxn_annotation_manual_ko = self.curation_api.get_manual_ko(rxn_id, template_id)
        
        #load functions
        rxn_annotation = self.annotation_api.get_reaction_annotation_data3(
        rxn_id, genome_set_id, 10, 
        rxn_annotation_manual_ko['ko'],
        rxn_annotation_manual_function['functions'])
        
        
        #Look for other reactions
        reaction_template_id = '{}@{}'.format(rxn_id, template_id)
        rxn_annotation_curation = self.curation_api.collection_templates_reactions.find_one({'_id' : reaction_template_id})
    
        rxn_annotation_functions_rxn = {}
        function_ids = set()
        for name in rxn_annotation:
            function_ids.add(rxn_annotation[name]['id'])
        for function_id in function_ids:
            res = self.curation_api.get_rxn_with_function(function_id, template_id)
            if res == None:
                res = {}
            rxn_annotation_functions_rxn[function_id] = res
            
        #rxn = modelseed_local.get_seed_reaction(rxn_id)
        #if rxn == None:
        #    rxn = {}
            
        response = {
            #'rxn' : {rxn_id : clear_nan(rxn.data)},
            'cpd' : {},
            'manual_function' : rxn_annotation_manual_function,
            'manual_ko' : rxn_annotation_manual_ko,
            'annotation' : rxn_annotation,
            'curation' : rxn_annotation_curation,
            'function_rxns' : rxn_annotation_functions_rxn
        }
        
        return response