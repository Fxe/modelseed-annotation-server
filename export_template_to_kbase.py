import logging
import copy
from cobrakbase.core.kbasefba import NewModelTemplate, TemplateManipulator, TemplateCuration
from cobrakbase.core.kbasefba.newmodeltemplate_validator import NewModelTemplateValidator
from cobrakbase.core.kbasegenomesgenome import normalize_role
from cobrakbase.core.kbasefba.new_template_reaction import NewModelTemplateReaction

logger = logging.getLogger(__name__)


def export_template(template_o, modelseed, annotation_api, mongo_database,
                    annotation_namespace='fungi',
                    reaction_list=None,
                    clear_reactions=False, clear_complexes=False, clear_roles=False):
    
    data_copy = copy.deepcopy(template_o.get_data())
    template = NewModelTemplate(data_copy, template_o.info, None, 'tftr', 'tcpx')

    tm = TemplateManipulator(template, None)
    template_reactions_filter = tm.clean_template('ModelSEED')
    print(len(template_reactions_filter))
    updated, removed = tm.clear_orphan_roles()
    print('updated', len(updated))
    print('removed', len(removed))

    validator = NewModelTemplateValidator(template)
    validator.validate_compounds()
    validator.validate()
    print('undeclared compounds', len(validator.undec_compounds))
    print('undeclared roles', len(validator.undec_roles))
    print('undeclared complexes', len(validator.undec_complexes))

    if clear_reactions or reaction_list is not None:
        template.reactions.clear()
        template.reactions._dict.clear()
    if clear_roles:
        template.complexes.clear()
        template.complexes._dict.clear()
        template.roles.clear()
        template.roles._dict.clear()
        template.role_set_to_cpx.clear()
        template.search_name_to_role_id.clear()
        template.role_last_id = 0
        template.complex_last_id = 0
    if clear_complexes:
        template.complexes.clear()
        template.complexes._dict.clear()
        template.role_set_to_cpx.clear()
        template.complex_last_id = 0
    # if roles or complexes were removed clear from reactions
    if clear_roles or clear_complexes:
        for trxn in template.reactions:
            trxn.templatecomplex_refs.clear()


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
    if reaction_list is not None:  # filter curation actions if reaction list is provided
        accept = dict(filter(lambda x: x[0] in reaction_list, accept.items()))
        remove = dict(filter(lambda x: x[0] in reaction_list, remove.items()))
    test_accept = dict(filter(lambda x: len(x[1]) > 0, accept.items()))
    test_remove = dict(filter(lambda x: len(x[1]) > 0, remove.items()))

    # add new roles to template

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

    # remove all reactions exclude from template
    remove_reactions = tc.get_disabled_reactions(annotation_namespace)

    reactions_in_template = set(map(lambda x: x.id, template.reactions))

    # strip complexes from reactions in remove set
    for trxn_id in remove_reactions:
        if trxn_id in reactions_in_template:
            trxn = template.get_reaction(trxn_id)
            trxn.templatecomplex_refs.clear()

    # add new reactions
    reactions_to_add = []
    for doc in tc.curation_api['templates_reactions'].find():
        template_rxn_id, template_id = doc['_id'].split('@')
        if template_id == annotation_namespace and template_rxn_id not in remove_reactions and \
                (reaction_list is None or template_rxn_id in reaction_list):
            if 'cmp' in doc:
                if template_rxn_id not in reactions_in_template and \
                        'annotation' in doc and \
                        'seed__DOT__reaction' in doc['annotation']:
                    seed_id = doc['annotation']['seed__DOT__reaction']
                    trxn_b = tm.build_template_reaction_from_modelseed(seed_id, doc['cmp'])
                    reactions_to_add.append(NewModelTemplateReaction(trxn_b))
    template.reactions += reactions_to_add
    reactions_in_template = set(map(lambda x: x.id, template.reactions))

    for trxn_id in set(test_remove):
        if trxn_id in reactions_in_template and trxn_id not in remove_reactions:
            template_rxn = template.get_reaction(trxn_id)
            role_change = tc.get_role_change(trxn_id, {}, test_remove)
            # role_change = get_role_change2(tc, rxn_id, {}, test_remove)
            # print(trxn.id)
            nfunction = tc.update_roles(template_rxn, role_change, search_name_to_role_id, True)
        else:
            logger.debug('%s', trxn_id)


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

    for trxn_id in set(test_accept):
        if trxn_id in reactions_in_template and trxn_id not in remove_reactions  and \
                (reaction_list is None or trxn_id in reaction_list):
            template_rxn = template.get_reaction(trxn_id)
            try:
                role_change = tc.get_role_change(trxn_id, test_accept, {})
                # role_change = get_role_change2(tc, rxn_id, test_accept, {})
                nfunction = tc.update_roles(template_rxn, role_change, search_name_to_role_id, True)
            except Exception as e:
                print(template_rxn.id, e)
        else:
            logger.debug('%s', trxn_id)


    reaction_annotation = tc.get_reaction_annotation()
    role_uids = set()
    system_accept_role_uids = set()
    for rxn_id in reaction_annotation[annotation_namespace]:
        if reaction_list is None or rxn_id in reaction_list:
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

