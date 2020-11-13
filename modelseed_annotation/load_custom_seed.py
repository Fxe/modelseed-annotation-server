import pandas as pd
import copy
import math


def to_stoichiometry(cs, ms):
    l = []
    for p in cs:
        cpd = ms.get_seed_compound(p[0])
        ss = '{}:{}:{}:{}:"{}"'.format(cs[p], p[0], p[1], '0', cpd.name)
        l.append(ss)
    return ";".join(l)


def to_equation(cs):
    lhs = []
    rhs = []
    for p in cs:
        coeff = cs[p]
        if coeff < 0:
            lhs.append("({}) {}[{}]".format(math.fabs(coeff), p[0], p[1]))
        else:
            rhs.append("({}) {}[{}]".format(math.fabs(coeff), p[0], p[1]))
    return "{} <=> {}".format(' + '.join(lhs), ' + '.join(rhs))


def to_definition(cs, ms):
    lhs = []
    rhs = []
    for p in cs:
        coeff = cs[p]
        cpd = ms.get_seed_compound(p[0])
        if coeff < 0:
            lhs.append("({}) {}[{}]".format(math.fabs(coeff), cpd.name, p[1]))
        else:
            rhs.append("({}) {}[{}]".format(math.fabs(coeff), cpd.name, p[1]))
    return "{} <=> {}".format(' + '.join(lhs), ' + '.join(rhs))


base = {
    'id': 'rxn00001',
    'abbreviation': 'R00004',
    'name': 'diphosphate phosphohydrolase',
    'code': '(1) cpd00001[0] + (1) cpd00012[0] <=> (2) cpd00009[0]',
    'stoichiometry': '-1:cpd00001:0:0:"H2O";-1:cpd00012:0:0:"PPi";2:cpd00009:0:0:"Phosphate";1:cpd00067:0:0:"H+"',
    'equation': '(1) cpd00001[0] + (1) cpd00012[0] <=> (2) cpd00009[0] + (1) cpd00067[0]',
    'definition': '(1) H2O[0] + (1) PPi[0] <=> (2) Phosphate[0] + (1) H+[0]',
    'is_transport': 0,
    'reversibility': '>',
    'direction': '=',
    'abstract_reaction': None,
    'pathways': 'MetaCyc: Degradation (Degradation/Utilization/Assimilation); Glyphosate-Degradation (glyphosate degradation); Noncarbon-Nutrients (Inorganic Nutrient Metabolism); PWY-7805 ((aminomethyl)phosphonate degradation); PWY-7807 (glyphosate degradation III); Phosphorus-Compounds (Phosphorus Compound Metabolism)',
    'aliases': 'AraCyc: INORGPYROPHOSPHAT-RXN|BiGG: IPP1; PPA; PPA_1; PPAm|BrachyCyc: INORGPYROPHOSPHAT-RXN|KEGG: R00004|MetaCyc: INORGPYROPHOSPHAT-RXN|Name: Diphosphate phosphohydrolase; Inorganic diphosphatase; Inorganic pyrophosphatase; Pyrophosphate phosphohydrolase; diphosphate phosphohydrolase; inorganic diphosphatase; inorganic diphosphatase (one proton translocation); inorganicdiphosphatase; pyrophosphate phosphohydrolase',
    'ec_numbers': '3.6.1.1',
    'deltag': -3.46,
    'deltagerr': 0.05,
    'compound_ids': 'cpd00001;cpd00009;cpd00012;cpd00067',
    'status': 'OK',
    'is_obsolete': 0,
    'linked_reaction': 'rxn27946;rxn27947;rxn27948;rxn32487;rxn38157;rxn38158',
    'notes': 'GCC|HB|EQC|EQU',
    'source': 'Primary Database'
}

def load_custom_seed(path, ms):
    rxns = {}
    df = pd.read_csv(path, sep='\t')
    for row_id, d in df.iterrows():
        rxn_id = d['rxn_id']
        name = d['Name']
        cpd_id = d['cpd_id']
        coeff = d['coeff']
        cmp_token = d['cmp']
        if rxn_id not in rxns:
            rxns[rxn_id] = copy.deepcopy(base)
            rxns[rxn_id]['metabolites'] = {}
        rxns[rxn_id]['id'] = rxn_id
        rxns[rxn_id]['name'] = name
        rxns[rxn_id]['abbreviation'] = rxn_id
        cc = (cpd_id, cmp_token)
        if cc not in rxns[rxn_id]['metabolites']:
            rxns[rxn_id]['metabolites'][cc] = coeff

    for rxn_id in rxns:
        cs = rxns[rxn_id]['metabolites']
        rxns[rxn_id]['code'] = to_equation(cs)
        rxns[rxn_id]['equation'] = to_equation(cs)
        rxns[rxn_id]['stoichiometry'] = to_stoichiometry(cs, ms)
        rxns[rxn_id]['definition'] = to_definition(cs, ms)
        rxns[rxn_id]['compound_ids'] = ';'.join(set(map(lambda x: x[0], cs)))
        rxns[rxn_id]['linked_reaction'] = ''
        rxns[rxn_id]['pathways'] = 'MetaCyc: Photosynthesis light reactions'
        rxns[rxn_id]['aliases'] = ''
        del rxns[rxn_id]['metabolites']
        ms.reactions[rxn_id] = copy.deepcopy(rxns[rxn_id])
