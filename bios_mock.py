import logging
import json

logger = logging.getLogger(__name__)

class BIOS_MOCK:
    
    def __init__(self, data_path):
        logger.info('loading cache data: %s', data_path)
        with open(data_path, 'r') as f:
            self.all_data = json.loads(f.read())
        
    
    def get_model_reactions(self, model_id):
        if not model_id in self.all_data:
            return []
        return self.all_data[model_id]['rxn']
    
    def get_model_species(self, model_id):
        if not model_id in self.all_data:
            return []
        return self.all_data[model_id]['spi']
    
    def get_model_compartments(self, model_id):
        if not model_id in self.all_data:
            return []
        return self.all_data[model_id]['cmp']