from bios_mock import BIOS_MOCK
import json
import pandas as pd


def load_cache_data(folder):
    bios = BIOS_MOCK(folder + 'bios_cache_fungi.json')
    MODEL_CPD_MAPPING = None
    model_rxn_mapping = {}
    model_cpd_mapping = {}
    model_rxn_gpr = None
    with open(folder + 'cpd_mapping_cache4.json', 'r') as f:
        model_cpd_mapping = json.loads(f.read())
    with open(folder + 'rxn_mapping_cache4.json', 'r') as f:
        model_rxn_mapping = json.loads(f.read())
    with open(folder + 'cmp_mapping_cache.json', 'r') as f:
        MODEL_CMP_MAPPING = json.loads(f.read())
    with open(folder + 'model_rxn_gpr_cache.json', 'r') as f:
        model_rxn_gpr = json.loads(f.read())
    df = pd.read_csv(folder + '/etc_mapping.tsv', sep='\t')
    for row_id, d in df.iterrows():
        seed_id = d['ModelSEED']
        for k in d.keys():
            if k in model_rxn_mapping:
                if not pd.isna(d[k]):
                    for mrxn_id in d[k].split(';'):
                        if mrxn_id in model_rxn_mapping[k] and seed_id not in model_rxn_mapping[k][mrxn_id]:
                            # print('miss', k, mrxn_id, seed_id, rxn_mapping_cache[k][mrxn_id])
                            model_rxn_mapping[k][mrxn_id] = [seed_id]
                        elif mrxn_id in model_rxn_mapping[k] and seed_id in model_rxn_mapping[k][mrxn_id]:
                            # print('ok', k, mrxn_id, seed_id, rxn_mapping_cache[k][mrxn_id])
                            pass
                        else:
                            # print('set', k, mrxn_id, seed_id)
                            model_rxn_mapping[k][mrxn_id] = [seed_id]

    bios.model_rxn_mapping = model_rxn_mapping
    bios.model_cpd_mapping = model_cpd_mapping
    return bios, MODEL_CMP_MAPPING, MODEL_CPD_MAPPING, model_rxn_mapping, model_rxn_gpr

def clear_nan(d):
    for k in d:
        if type(d[k]) == float and pd.isna(d[k]):
            d[k] = ""
