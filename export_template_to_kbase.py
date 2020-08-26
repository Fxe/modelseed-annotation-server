import logging
import copy
from cobrakbase.core.kbasefba import NewModelTemplate, TemplateManipulator, TemplateCuration
from cobrakbase.core.kbasefba.newmodeltemplate_validator import NewModelTemplateValidator
from cobrakbase.core.kbasegenomesgenome import normalize_role

logger = logging.getLogger(__name__)


def export_template(template_o, modelseed, annotation_api, mongo_database, annotation_namespace='fungi'):
    
    from cobra.core.dictlist import DictList
    temp_object = {}
    for k in template_o.data.keys():
        if k not in ['data', 'info', 'provenance']:
            if type(template_o.data[k]) is DictList:
                temp_object[k] = list(template_o.data[k])
            else:
                temp_object[k] = template_o.data[k]
    template = NewModelTemplate(copy.deepcopy(temp_object), template_o.info, None, 'tftr', 'tcpx')
    #template = NewModelTemplate(copy.deepcopy(template_o), 'tftr', 'tcpx')
    validator = NewModelTemplateValidator(template)
    validator.validate_compounds()
    validator.validate()
    print('undeclared compounds', len(validator.undec_compounds))
    print('undeclared roles', len(validator.undec_roles))
    print('undeclared complexes', len(validator.undec_complexes))

    tm = TemplateManipulator(template, None)
    template_reactions_filter = tm.clean_template('ModelSEED')
    len(template_reactions_filter)
    updated, removed = tm.clear_orphan_roles()
    print('updated', len(updated))
    print('removed', len(removed))

    validator = NewModelTemplateValidator(template)
    validator.validate_compounds()
    validator.validate()
    print('undeclared compounds', len(validator.undec_compounds))
    print('undeclared roles', len(validator.undec_roles))
    print('undeclared complexes', len(validator.undec_complexes))

    tc = TemplateCuration(template, mongo_database, annotation_api)
    tm = TemplateManipulator(template, modelseed)
    a = tc.get_reaction_annotation()
    print(a.keys())
    print(len(a[annotation_namespace]))
    #dict_keys(['template_v3', 'fungi'])
    #6257
    search_name_to_role_id = tm.get_search_name_to_role_id()
    for k in search_name_to_role_id:
        if len(search_name_to_role_id[k]) > 1:
            print(k)
    accept, remove = tc.get_curation_data(annotation_namespace)
    test_accept = dict(filter(lambda x : len(x[1]) > 0, accept.items()))
    test_remove = dict(filter(lambda x : len(x[1]) > 0, remove.items()))

    test_accept_sn_to_roles = tc.get_roles_to_add(test_accept, search_name_to_role_id)
    for role_sn in test_accept_sn_to_roles:
        role_name = list(test_accept_sn_to_roles[role_sn])[0]
        template.add_role(role_name)
    search_name_to_role_id = tm.get_search_name_to_role_id()

    def get_compartment_token(cmp_config):
        v = cmp_config.values()
        if len(v) == 1:
            return list(v)[0]
        if len(v) == 2:
            if 'e' in v and 'c' in v:
                return 'c'
            elif 'c' in v:
                return list(filter(lambda x: not x == 'c', v))[0]
        return None

    for doc in tc.curation_api['templates_reactions'].find():
        rxn_id, template_id = doc['_id'].split('@')
        if template_id == annotation_namespace:
            if 'cmp' in doc:
                cmp_id = get_compartment_token(doc['cmp'])
                # print(rxn_id, cmp_id)
                trxn = template.get_reaction(rxn_id + '_' + cmp_id)
                if trxn == None:
                    tm.add_reaction(rxn_id, doc['cmp'])

    cmp = 'c'
    for rxn_id in set(test_remove):
        template_rxn = template.get_reaction(rxn_id + '_' + cmp)
        if template_rxn is None:
            logger.warning('%s', rxn_id)
        else:
            role_change = tc.get_role_change(rxn_id, {}, test_remove)
            # role_change = get_role_change2(tc, rxn_id, {}, test_remove)
            # print(trxn.id)
            nfunction = tc.update_roles(template_rxn, role_change, search_name_to_role_id, True)

    def refresh(template):
        template.role_set_to_cpx = {}
        template.search_name_to_role_id = {}
        for role in template.data['roles']:
            template.search_name_to_role_id[normalize_role(role['name'])] = role['id']
        for cpx in template.data['complexes']:
            roles = set()
            for complexrole in cpx['complexroles']:
                role_id = complexrole['templaterole_ref'].split('/')[-1]
                roles.add(role_id)
            # print(cpx, roles)
            template.role_set_to_cpx[';'.join(sorted(roles))] = cpx['id']

    refresh(template)

    cmp = 'c'
    for rxn_id in set(test_accept):
        template_rxn = template.get_reaction(rxn_id + '_' + cmp)
        if template_rxn is None:
            logger.warning('%s', rxn_id)
        else:
            try:
                role_change = tc.get_role_change(rxn_id, test_accept, {})
                # role_change = get_role_change2(tc, rxn_id, test_accept, {})
                nfunction = tc.update_roles(template_rxn, role_change, search_name_to_role_id, True)
            except Exception as e:
                print(template_rxn.id, e)

    reaction_annotation = tc.get_reaction_annotation()
    role_uids = set()
    system_accept_role_uids = set()
    for rxn_id in reaction_annotation[annotation_namespace]:
        for role_id in reaction_annotation[annotation_namespace][rxn_id]['user']:
            score = reaction_annotation[annotation_namespace][rxn_id]['current'][str(role_id)]
            if score == 'opt_score1':
                user_log = reaction_annotation[annotation_namespace][rxn_id]['user'][role_id]
                if len(user_log) > 1 or 'system' not in user_log:
                    role_uids.add(role_id)
                    # print(rxn_id, role_id, user_log, score)
                else:
                    system_accept_role_uids.add(role_id)

    print(len(role_uids))
    for role_uid in role_uids:
        a, role_ids = tc.get_function(role_uid, search_name_to_role_id)
        for role_id in role_ids:
            role = template.get_role(role_id)
            if role:
                # print(role_uid, role['id'], role['source'])
                role['source'] = 'ModelSEED'

    return template

