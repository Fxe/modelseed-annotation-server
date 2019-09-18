from neo4j import GraphDatabase

class AnnotationApiNeo4j:
    
    def __init__(self, driver):
        self.driver = driver

    def init_constraints(self):
        unique_value_key = [
            'KeggGenome', 'RefSeqGenome', 'GenBankGenome', 
            'ModelSeedReaction', 'LigandReaction',
            'KeggOrthology'
            'FunctionSource', 
            'FunctionGroup', 
            'Function',
            'FunctionComplex',
            #'Genome', 
            'KBaseGene', 
            #'Reaction',
            'Fruit',
            'TemplateSet', 'TemplateReactionAnnotation',
        ]
        with self.driver.session() as session:
            for l in unique_value_key:
                session.write_transaction(
                    lambda tx : tx.run(
                        'CREATE CONSTRAINT ON (n:{}) ASSERT n.key IS UNIQUE'.format(l)))
                
    def get_function(self, f):
        return self.get_node('Function', f)
                
    def get_node(self, t, key):
        result = None
        with self.driver.session() as session:
            query = session.read_transaction(self._get_node, t, key)
            if len(query) > 0:
                result = query[0].data()['n']
        return result
    
    def get_genome(self, genome_id):
        return self.get_node('RefSeqGenome', genome_id)
    
    def add_node(self, node_type, key, props = {}):
        props['key'] = key
        props_str = self.to_str_dict(props)
        result = None
        with self.driver.session() as session:
            query = session.write_transaction(self._add_value_node, node_type, props_str)
            if len(query) > 0:
                result = query[0].data()['n']
        return result
    
    def add_annotation(self, gene_genome_id, functions, source):
        #[KBaseGene] -[has_function]-> [Function] <-[has_function]- [FunctionSource]
        
        kbase_gene_node = self.get_node('KBaseGene', gene_genome_id)
        function_source_node = self.get_node('FunctionSource', source)

        if function_source_node == None:
            print('not found', source)
            return    
        if kbase_gene_node == None:
            print('not found', gene_genome_id)
            return

        for f in functions:
            if not f == None:
                function_node = self.add_node('Function', f)
                self.link_nodes(function_source_node, function_node, 'has_function')
                self.link_nodes(kbase_gene_node, function_node, 'has_function')
    
    def link_nodes(self, node1, node2, relationship):
        out = None
        with self.driver.session() as session:
            out = session.write_transaction(self._add_relationship, 
                                            node1, node2, relationship)
        return out

    @staticmethod
    def to_str_dict(props):
        str_dict = ""
        if not props == None and len(props) > 0:
            dict_data = []
            for k in props:
                data = k + ':'
                #print(k, type(props[k]))
                if type(props[k]) == str:
                    data += '"' + props[k].replace('"', '\\"') + '"'
                elif type(props[k]) == int:
                    data += str(props[k])
                elif type(props[k]) == float:
                    data += str(props[k])
                else:
                    print(k, type(props[k]))
                dict_data.append(data)
            str_dict = '{' + ', '.join(dict_data) + '}'
            pass
        return str_dict
            
    @staticmethod
    def _add_value_node(tx, label, props_str = ""):
        result = tx.run("MERGE (n:" + label + " " + props_str + ") "
                        "ON CREATE SET n.created_at = timestamp(), n.updated_at = timestamp()"
                        "ON MATCH  SET n.updated_at = timestamp() "
                        "RETURN n")
        return list(result.records())
    
    @staticmethod
    def _get_node(tx, label, key):
        result = tx.run("MATCH (n:" + label + " {key:$key}) "
                        "RETURN n", key = key)
        return list(result.records())
    
    @staticmethod
    def _add_relationship(tx, node1, node2, relationship):
        result = tx.run("MATCH (n1) WHERE id(n1)=$node1_id "
                        "WITH n1 "
                        "MATCH (n2) WHERE id(n2)=$node2_id "
                        "CREATE UNIQUE (n1)-[r:" + relationship + "]->(n2) ", 
                        node1_id=node1.id, node2_id=node2.id)
        return result
    
    def add_template_reaction_annotation(self, template, rxn, cmps, function_rule, props):

        trxn_ann_node_id = rxn.get('key') + '_' + '_'.join(sorted(list(cmps)))

        trxn_ann_node = self.get_node('TemplateReactionAnnotation', trxn_ann_node_id)

        if trxn_ann_node == None:
            trxn_ann_node = self.add_node('TemplateReactionAnnotation', trxn_ann_node_id, props)
            for rule in function_rule:
                self.link_nodes(trxn_ann_node, rule, 'has_function')

            self.link_nodes(template, trxn_ann_node, 'has_annotation_rule')
            self.link_nodes(trxn_ann_node, rxn, 'has_reaction')

        return trxn_ann_node
    
    def add_template(self, template_id, template):
        roles = {}
        function_uids = {}
        function_complexes = {}
        compcompound_compartment = {}

        for compcompound in template['compcompounds']:
            compcompound_compartment[compcompound['id']] = compcompound['templatecompartment_ref'].split('/')[-1]

        for role in template['roles']:
            roles[role['id']] = role['name']
        for role_id in roles:
            function = roles[role_id]
            function_data = self.get_function(function)
            if not function_data == None and function_data.get('key') == function:
                function_uids[function] = function_data
            else:
                #print('errooo!')
                pass

        missing_function = set()
        found_function = set()
        for role_id in roles:
            function = roles[role_id]
            if not function in function_uids:
                missing_function.add(function)
            else:
                found_function.add(function)

        logger.warning('%d %d', len(found_function), len(missing_function))

        for f in missing_function:
            self.add_node('Function', f, {})

        cpx_to_function_uids = {}
        for cpx in template['complexes']:
            uids = set()
            for complexrole in cpx['complexroles']:
                role_id = complexrole['templaterole_ref'].split('/')[-1]
                if role_id in roles:
                    uids.add(function_uids[roles[role_id]])
                else:
                    print('errooo!')
            cpx_to_function_uids[cpx['id']] = uids
        #add_template(annotation_api, 'GramNegModelTemplate', templates['GramNegModelTemplate'])

        for cpx_id in cpx_to_function_uids:
            if len(cpx_to_function_uids[cpx_id]) > 1:
                complex_node = self.get_node('FunctionComplex', cpx_id)
                if complex_node == None:
                    complex_node = self.add_node('FunctionComplex', cpx_id)
                    for n in cpx_to_function_uids[cpx_id]:
                        self.link_nodes(complex_node, n, 'has_function')
                    function_complexes[cpx_id] = complex_node

        reactions = {}
        for rxn in template['reactions']:
            rxn_id = rxn['id']
            if '_' in rxn_id:
                rxn_id = rxn_id.split('_')[0]
            rxn_id += '@' + template_id
            rxn_node = self.get_node('ModelSeedReaction', rxn_id)
            if not rxn_node == None:
                #print(rxn_node)
            #rxn_id = rxn['reaction_ref'].split('/')[-1]
                rxn_compartments = set()
                for templateReactionReagent in rxn['templateReactionReagents']:
                    templatecompcompound_ref = templateReactionReagent['templatecompcompound_ref'].split('/')[-1]
                    rxn_compartments.add(compcompound_compartment[templatecompcompound_ref])



                or_rule = set()
                for complex_ref in rxn['templatecomplex_refs']:
                    complex_id = complex_ref.split('/')[-1]
                    if complex_id in cpx_to_function_uids:
                        if len(cpx_to_function_uids[complex_id]) > 1:
                            or_rule.add(function_complexes[complex_id])
                        elif len(cpx_to_function_uids[complex_id]) == 1:
                            or_rule.add(list(cpx_to_function_uids[complex_id])[0])
                        else:
                            logger.warning('%s', complex_id)

                #print(rxn['id'], rxn_compartments, or_rule)

                props = {
                    'compartment' : ';'.join(sorted(list(rxn_compartments))),
                    'reaction' : rxn_node.get('key'),
                    'base_cost' : rxn['base_cost'],
                    'direction' : rxn['direction'],
                    'forward_penalty' : rxn['forward_penalty'],
                    'reverse_penalty' : rxn['reverse_penalty'],
                    'type' : rxn['type'],
                    'GapfillDirection' : rxn['GapfillDirection'],
                }

                self.add_template_reaction_annotation(template_node, rxn_node, 
                                                      rxn_compartments, or_rule, props)