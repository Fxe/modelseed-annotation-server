import logging
import json
from biosapi.core.model import BiosModelReaction

logger = logging.getLogger(__name__)


class BIOS_MOCK:
    
    def __init__(self, data_path):
        logger.info('loading cache data: %s', data_path)
        self.model_data = {}
        self.model_spi = {}
        self.model_rxn = {}
        self.model_cmp = {}
        self.model_rxn_mapping = {}
        self.model_cpd_mapping = {}
        with open(data_path, 'r') as f:
            self.model_data = json.loads(f.read())
            
        self.model_spi, self.model_rxn = self.index_data(self.model_data)
        
    def index_data(self, model_data):
        model_spi = {}
        model_rxn = {}
        for model_id in model_data:
            model_spi[model_id] = {}
            model_rxn[model_id] = {}
            for o in model_data[model_id]['spi']:
                if 'id' in o:
                    if o['id'] in model_spi[model_id]:
                        print('1')
                    model_spi[model_id][o['id']] = o
                else:
                    logger.warning('[%s] missing id', model_id)
            for o in model_data[model_id]['rxn']:
                if 'id' in o:
                    if o['id'] in model_rxn[model_id]:
                        print('1')
                    model_rxn[model_id][o['id']] = o
                else:
                    logger.warning('[%s] missing id', model_id)
        return model_spi, model_rxn
        
    def get_models(self):
        return None
    
    def get_model_reactions(self, model_id):
        if not model_id in self.model_data:
            return []
        return self.model_data[model_id]['rxn']
    
    def get_model_species(self, model_id):
        if not model_id in self.model_data:
            return []
        return self.model_data[model_id]['spi']
    
    def get_model_compartments(self, model_id):
        if not model_id in self.model_data:
            return []
        return self.model_data[model_id]['cmp']
    
    def get_model_specie(self, model_id, spi_id):
        if model_id in self.model_spi and spi_id in self.model_spi[model_id]:
            return self.model_spi[model_id][spi_id]
        return None
    
    def get_model_reaction(self, model_id, rxn_id):
        if model_id in self.model_rxn and rxn_id in self.model_rxn[model_id]:
            return BiosModelReaction(self.model_rxn[model_id][rxn_id])
        return None
    
    def get_model_reactions_by_database_id(self, model_id, rxn_id, database):
        res = []
        if model_id in self.model_rxn_mapping and model_id in self.model_rxn:
            for mrxn_id in self.model_rxn_mapping[model_id]:
                if rxn_id in self.model_rxn_mapping[model_id][mrxn_id]:
                    res.append((BiosModelReaction(self.model_rxn[model_id][mrxn_id]), 'mock:5'))
        return res