import copy
from cobrakbase.core.kbasefba import NewModelTemplate
from cobrakbase.core.kbasefba import TemplateManipulator
from cobrakbase.core.kbasefba import TemplateCuration
from cobrakbase.core.kbasefba.newmodeltemplate_validator import NewModelTemplateValidator


def pre_process_template(template, mongo_database, modelseed, annotation_api):
    tm = TemplateManipulator(template, None)
    template_reactions_filter = tm.clean_template('ModelSEED')
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
    print(len(a['fungi']))
    #dict_keys(['template_v3', 'fungi'])
    #6257
    search_name_to_role_id = tm.get_search_name_to_role_id()
    for k in search_name_to_role_id:
        if len(search_name_to_role_id[k]) > 1:
            print(k)
    accept, remove = tc.get_curation_data('fungi')
    test_accept = dict(filter(lambda x : len(x[1]) > 0, accept.items()))
    test_remove = dict(filter(lambda x : len(x[1]) > 0, remove.items()))


def wut(kbase, mongo_database, annotation_api):
    template_o = kbase.get_object('template_v2.x_06102020', 'filipeliu:narrative_1582914694010')
    template = NewModelTemplate(copy.deepcopy(template_o), 'tftr', 'tcpx')
    tc = TemplateCuration(template, mongo_database, annotation_api)
    tm = TemplateManipulator(template)
