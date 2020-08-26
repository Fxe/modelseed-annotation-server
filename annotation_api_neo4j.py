import logging
from neo4j import GraphDatabase
from py2neo import Graph, NodeMatcher, RelationshipMatcher
from py2neo.data import Node, Relationship

from cobrakbase.core.kbasegenomesgenome import KBaseGenomeFeature
from cobrakbase.core.kbasegenomesgenome import normalize_role

logger = logging.getLogger(__name__)

def process_list(l, sep = ' ### ', subsep = ' ## '):
    j = []
    for o in l:
        if not type(o) == str:
            if type(o) == list:
                p = process_list(o, subsep, None)
                
                if not p == None:
                    j.append(p)
                else:
                    return None
            else:
                return None
        else:
            j.append(o)
    return sep.join(j)

def process_dict(d):
    valid_type = [str, int, float]
    props = {}
    for k in d:
        if type(d[k]) in valid_type:
            props[k] = d[k]
            #print(k, type())
        elif type(d[k]) == list:
            p = process_list(d[k])
            if p == None:
                logger.debug('excluded: %s -> %s', k, type(d[k]))
            else:
                props[k] = p 
        else:
            logger.debug('excluded: %s -> %s', k, type(d[k]))
    return props

class AnnotationFunction:
    
    def __init__(self, id, value):
        self.id = id
        self._value = value
        self.search_value = value
        self.synonyms = set()
        self.sub_functions = set()
        self.function_group = set()
        
    @property
    def value(self):
        return self._value
        
class Neo4jAnnotationFunction(AnnotationFunction):
    
    def __init__(self, node):
        super().__init__(node.identity, node['key'])
        self.node = node
        #print(node)
        self.populate_from_graph()
        
    def populate_from_graph(self):
        for rel in self.node.graph.match((self.node, ), r_type="has_subfunction", ):
            self.sub_functions.add(Neo4jAnnotationFunction(rel.end_node))
        for rel in self.node.graph.match((self.node, ), r_type="search_function", ):
            self.search_value = rel.end_node['key']
            for srel in self.node.graph.match((rel.end_node, ), r_type="has_function", ):
                n = srel.end_node
                if not n['key'] == self.value:
                    self.synonyms.add(n['key'])
                    
        for rel in self.node.graph.match((None, self.node), r_type="has_function", ):
            n = rel.start_node
            if n.has_label('FunctionGroup'):
                self.function_group.add(n['key'])

    def neo4j_load_source(self):
        self.source = set()
        for srel in self.node.graph.match((self.node, ), r_type="has_source", ):
            self.source.add(srel.end_node['key'])
        return self.source
                
    def __str__(self):
        return self.value

class AnnotationApiNeo4j:
    
    def __init__(self, driver=None, user=None, pwd=None, uri=None, port=7474, host="0.0.0.0"):
        self.driver = driver
        self.neo4j_graph = None
        self.matcher = NodeMatcher(self.neo4j_graph)
        if self.driver == None:
            self.driver = GraphDatabase.driver("bolt://" + host + ":" + str(port), auth=(user, pwd))
        if self.neo4j_graph == None:
            self.neo4j_graph = Graph("http://" + host + ":" + str(port), auth=(user, pwd))
            self.matcher = NodeMatcher(self.neo4j_graph)

    def init_constraints(self):
        unique_value_key = [
            'KeggGenome', 'RefSeqGenome', 'GenBankGenome', 
            'ModelSeedReaction', 'LigandReaction',
            'MetabolicModel',
            'KeggOrthology'
            'FunctionSource', 
            'FunctionGroup', 
            'Function',
            'GenomeSet',
            'SearchFunction',
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
        n = self.matcher.match("Function", key=f).first()
        if n == None:
            return None
        return Neo4jAnnotationFunction(n)
    
    def get_function_by_uid(self, function_id):
        n = self.neo4j_graph.nodes[int(function_id)]
        if n == None or not 'Function' in n.labels:
            return None
        return Neo4jAnnotationFunction(n)
                
    def get_node(self, t, key):
        result = None
        with self.driver.session() as session:
            query = session.read_transaction(self._get_node, t, key)
            if len(query) > 0:
                result = query[0].data()['n']
        return result
    
    def get_ko(self, ko_id):
        n = self.matcher.match("KeggOrthology", key=ko_id).first()
        if n == None:
            return None
        return n
    
    def get_genome(self, genome_id):
        return self.get_node('RefSeqGenome', genome_id)
    
    def get_annotation_function(self, function_str, search_str = False):
        node = None
        if search_str:
            nfunction_str = normalize_role(function_str)
            node = self.matcher.match('SearchFunction', key=nfunction_str).first()
        else:
            node = self.matcher.match('Function', key=function_str).first()

        if node == None:
            return None

        annotation_function = Neo4jAnnotationFunction(node)
        if node.has_label('SearchFunction'):
            annotation_function.search_value = annotation_function.value
            for rel in node.graph.match((None, node), r_type="search_function"):
                matching_function = Neo4jAnnotationFunction(rel.start_node)
                annotation_function.sub_functions.add(matching_function)
                annotation_function.synonyms.add(matching_function.value)
                for o in matching_function.synonyms:
                    annotation_function.synonyms.add(o)
                for o in matching_function.function_group:
                    annotation_function.function_group.add(o)

        #print(node)

        return annotation_function
    
    def update_node(self, node_type, key, props = {}):
        if len(props) == 0:
            return False
        matcher = NodeMatcher(self.neo4j_graph)
        node = matcher.match(node_type, key=key).first()
        if not node == None:
            node.update(props)
            tx = self.neo4j_graph.begin()
            tx.merge(node)
            tx.graph.push(node)
            tx.commit()
            return True
        return False
    
    def get_or_create(self, node_type, key, props = {}):
        n = self.get_node(node_type, key)
        if n == None:
            self.add_node(node_type, key, props)
            n = self.get_node(node_type, key)
            
        return n
    
    def add_node(self, node_type, key, props = {}):
        props['key'] = key
        props_str = self.to_str_dict(props)
        result = None
        with self.driver.session() as session:
            logger.debug('write_transaction %s %s', node_type, props_str)
            query = session.write_transaction(self._add_value_node, node_type, props_str)
            if len(query) > 0:
                result = query[0].data()['n']
        return result
    
    def add_node2(self, node_type, key, props = {}):
        props['key'] = key
        node = Node(node_type, **props)
        node.__primarylabel__ = node_type
        node.__primarykey__ = "key"
        self.neo4j_graph.merge(node)
        return self.get_node(node_type, key)
    
    def merge_node(self, node_type, key, props = {}):
        node = self.get_node(node_type, key)
        if node == None:
            node = self.add_node2(node_type, key, props)
        else:
            #update
            pass
        return node
    
    def get_ko_by_seed_id(self, seed_id):
        m = self.neo4j_graph.nodes.match("ModelSeedReaction", key=seed_id)
        if len(m) < 1:
            return None
        seed_node = m.first()
        kegg_nodes = set()
        for rel in self.neo4j_graph.match((seed_node, ), r_type="has_crossreference_to", ):
            kegg_nodes.add(rel.end_node)

        kos = set()
        for n in kegg_nodes:
            for rel in self.neo4j_graph.match((n, ), r_type="has_ortholog", ):
                kos.add(rel.end_node['key'])

        return kos
    
    def set_relationship_property(self, relationship_uid, prop, value):
        result = None
        with self.driver.session() as session:
            result = session.write_transaction(self._set_relationship_property, relationship_uid, prop, value)
        return result
    
    def add_annotation_old(self, gene_genome_id, functions, source):
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
    
    
    def page_nodes2(self, l, page, size, where = ""):
        """
        !!!
        """
        result = None
        count = None
        with self.driver.session() as session:
            query = session.read_transaction(self._page_nodes, l, page * size, size, where)
            if len(query) > 0:
                result = query
        with self.driver.session() as session:
            query = session.read_transaction(self._page_nodes_count, l, where)
            if len(query) > 0:
                count = query       
        
        return result, count
    
    def page_nodes(self, l, page, size, where = ""):
        """
        !!!
        """
        result = None
        with self.driver.session() as session:
            query = session.read_transaction(self._page_nodes, l, page * size, size, where)
            if len(query) > 0:
                result = query   
        
        return result

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
    def _set_relationship_property(tx, relationship_uid, prop, value):
        if type(value) == str:
            value = '"' + value.replace('"', '\\"') + '"'
        query = "MATCH ()-[r]->() WHERE id(r)= {} SET r.{} = {} RETURN r".format(relationship_uid, prop, value)
        result = tx.run(query)
        return list(result.records())
    
    @staticmethod
    def _update_value_node(tx, label, key, props_str = ""):
        result = tx.run("MATCH (n:" + label + " {key:$key}) "
                        "SET n = " + props_str + ""
                        " ON CREATE SET n.created_at = timestamp(), n.updated_at = timestamp() "
                        " ON MATCH  SET n.updated_at = timestamp() "
                        "RETURN n", key = key)
        return list(result.records())
    
    @staticmethod
    def _tx_run(tx, query, params):
        result = tx.run(query, params)
        return list(result.records())
    
    def page_genomes(self, skip, limit, taxa_filter = None):
        params = {
            'skip' : skip,
            'limit' : limit,
        }
        filter_function = ""
        filter_count = ""
        if not taxa_filter == None and len(taxa_filter) > 0:
            filter_function = " WHERE n.scientific_name =~ '.*" + taxa_filter + ".*' OR n.taxonomy =~ '.*" + taxa_filter + ".*'"
        query = \
            "MATCH (n:RefSeqGenome) WITH COUNT(n) AS total_genomes " + filter_count + " \
             MATCH (n:RefSeqGenome) " + filter_function + " RETURN n, total_genomes SKIP $skip LIMIT $limit"
        result = None
        with self.driver.session() as session:
            query = session.read_transaction(self._tx_run, query, params)
            if len(query) > 0:
                result = query
        return result
    
    def page_genome_genes(self, genome_id, skip, limit, function_filter = None):

        params = {
            'genome_id' : genome_id,
            'limit' : limit,
            'skip' : skip,
        }
        result = None
        filter_function = ""
        if not function_filter == None and len(function_filter) > 0:
            filter_function = "AND f.key =~ '.*" + function_filter + ".*'"
        query = \
        "MATCH (n:KBaseGene) WHERE n.key CONTAINS $genome_id WITH COUNT(n) AS total_genes \
         MATCH (n:KBaseGene)-[r:has_annotation]->(f:Function) WHERE n.key CONTAINS $genome_id " + filter_function + " RETURN n, collect(f.key) as function, collect(r.function_source) as function_source, total_genes SKIP $skip LIMIT $limit"
        with self.driver.session() as session:
            query = session.read_transaction(self._tx_run, query, params)
            if len(query) > 0:
                result = query

        return result
    
    @staticmethod
    def _page_nodes(tx, label, skip, page_size, where = ""):
        query = "MATCH (n:" + label + ") " + where + " " \
                "RETURN n ORDER BY n.key SKIP $skip LIMIT $page_size"
        #print(query)
        result = tx.run(query, skip=skip, page_size=page_size)
        return list(result.records())
    
    @staticmethod
    def _page_nodes_count(tx, label, where = ""):
        query = "MATCH (n:" + label + ") " + where + " " \
                "RETURN count(n) as count"
        #print(query)
        result = tx.run(query)
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
                        "CREATE UNIQUE (n1)-[r:" + relationship + "]->(n2) RETURN r", 
                        node1_id=node1.id, node2_id=node2.id)
        return result.data()[0]['r']
    
    @staticmethod
    def _fetch_gene_annotation_relationships(tx, function_str):
        query = "MATCH (n:KBaseGene)-[r:has_annotation]->(f:Function {key:$f})" \
                "RETURN n.key AS gene_id, r.function_source AS function_source"
        #ORDER BY n.key SKIP $skip LIMIT $page_size
        result = tx.run(query, f=function_str)
        return list(result.records())

    def fetch_gene_annotation_relationships(self, function_str):
        result = None
        with self.driver.session() as session:
            result = session.read_transaction(self._fetch_gene_annotation_relationships, function_str)
        return result
    
    def add_template_reaction_annotation(self, template, rxn, cmps, function_rule, props):

        trxn_ann_node_id = rxn.get('key') + '_' + '_'.join(sorted(list(cmps))) + '@' + template['key']

        trxn_ann_node = self.get_node('TemplateReactionAnnotation', trxn_ann_node_id)

        if trxn_ann_node == None:
            trxn_ann_node = self.add_node('TemplateReactionAnnotation', trxn_ann_node_id, props)
            
        for rule in function_rule:
            self.link_nodes(trxn_ann_node, rule, 'has_function')

        self.link_nodes(template, trxn_ann_node, 'has_annotation_rule')
        self.link_nodes(trxn_ann_node, rxn, 'has_reaction')
        self.link_nodes(rxn, trxn_ann_node, 'has_annotation_rule')
        
        return trxn_ann_node
    
    def add_template(self, template_id, template):
        roles = {}
        function_uids = {}
        function_complexes = {}
        compcompound_compartment = {}

        for compcompound in template['compcompounds']:
            compcompound_compartment[compcompound['id']] = compcompound['templatecompartment_ref'].split('/')[-1]

        template_node = self.get_node('TemplateSet', template_id)
        if template_node == None:
            self.add_node('TemplateSet', template_id, {})
            template_node = self.get_node('TemplateSet', template_id)

        for role in template['roles']:
            roles[role['id']] = role['name']
        for role_id in roles:
            function = roles[role_id]
            function_data = self.get_function(function)
            if not function_data == None and function_data.value == function:
                function_uids[function] = function_data
            else:
                #print('errooo!')
                pass
        
        #make functions not in database
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
            function_data = self.add_node('Function', f, {})
            if not function_data == None and function_data.get('key') == f:
                function_uids[f] = function_data
            
        cpx_to_function_uids = {}
        for cpx in template['complexes']:
            uids = set()
            for complexrole in cpx['complexroles']:
                role_id = complexrole['templaterole_ref'].split('/')[-1]
                if role_id in roles:
                    uids.add(function_uids[roles[role_id]])
                else:
                    print('errooo!')
            cpx_id = "{}@{}".format(cpx['id'], template_id)
            cpx_to_function_uids[cpx_id] = uids
            
        #make complex nodes for true complexes (functions > 1)
        for cpx_id in cpx_to_function_uids:
            if len(cpx_to_function_uids[cpx_id]) > 1:
                complex_node = self.get_node('FunctionComplex', cpx_id)
                if complex_node == None:
                    complex_node = self.add_node('FunctionComplex', cpx_id)
                    for n in cpx_to_function_uids[cpx_id]:
                        self.link_nodes(complex_node, n, 'has_function')
                    function_complexes[cpx_id] = complex_node
                else:
                    function_complexes[cpx_id] = complex_node
                    
        for rxn in template['reactions']:
            rxn_id = rxn['id']
            if '_' in rxn_id:
                rxn_id = rxn_id.split('_')[0]

            rxn_node = self.get_node('ModelSeedReaction', rxn_id)
            rxn_id += '@' + template_id
            if not rxn_node == None:
                rxn_compartments = set()
                for templateReactionReagent in rxn['templateReactionReagents']:
                    templatecompcompound_ref = templateReactionReagent['templatecompcompound_ref'].split('/')[-1]
                    rxn_compartments.add(compcompound_compartment[templatecompcompound_ref])

                or_rule = set()
                for complex_ref in rxn['templatecomplex_refs']:
                    complex_id = complex_ref.split('/')[-1] + '@' + template_id
                    if complex_id in cpx_to_function_uids:
                        if len(cpx_to_function_uids[complex_id]) > 1:
                            or_rule.add(function_complexes[complex_id])
                        elif len(cpx_to_function_uids[complex_id]) == 1:
                            or_rule.add(list(cpx_to_function_uids[complex_id])[0])
                        else:
                            logger.warning('%s', complex_id)

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
            else:
                logger.warning('reaction not found: %s', rxn['id'])
            
        return template_node
                
                
    def add_annotation(self, annotations, source):
        function_source_node = self.merge_node('FunctionSource', source)
        
        base_functions = set()
        #functions = KBaseGenomeFeature.split_annotation(annotation)
        for annotation in annotations:
            if not annotation == None and not len(annotation.strip()) == 0:
                annotation = annotation.strip()
                
                logger.debug('ADD BASE Function: %s', annotation)
                base_function_node = self.merge_node('Function', annotation)
                self.link_nodes(function_source_node, base_function_node, 'has_function')
                self.link_nodes(base_function_node, function_source_node, 'has_source')
                base_functions.add(base_function_node)
                functions = KBaseGenomeFeature.split_function(annotation)
                if len(functions) > 1:
                    for function in functions:
                        
                        logger.debug('ADD SUB Function: %s', function)
                        sub_function_node = self.merge_node('Function', function)
                        
                        self.link_nodes(base_function_node, sub_function_node, 'has_subfunction')
                        nfunction = normalize_role(function)
                        
                        logger.debug('ADD SEARCH SUB Function: %s', nfunction)
                        nsub_function_node = self.merge_node('SearchFunction', nfunction)
                        
                        self.link_nodes(sub_function_node, nsub_function_node, 'search_function')
                        self.link_nodes(nsub_function_node, sub_function_node, 'has_function')
                else:
                    nfunction = normalize_role(annotation)
                    
                    logger.debug('ADD SEARCH BASE Function: %s', nfunction)
                    nbase_function_node = self.merge_node('SearchFunction', nfunction)
                    
                    self.link_nodes(base_function_node, nbase_function_node, 'search_function')
                    self.link_nodes(nbase_function_node, base_function_node, 'has_function')

        return base_functions
    
    def add_genome_feature(self, feature, genome_id, annotation_source):

        feature_props = process_dict(feature.data)
        feature_id = "{}@{}".format(feature.id, genome_id)

        feature_node = self.get_node('KBaseGene', feature_id)
        if feature_node == None:
            logger.debug('ADD KBaseGene: %s', feature_id)
            feature_node = self.add_node2('KBaseGene', feature_id, feature_props)
        else:
            logger.debug('UPDATE KBaseGene: %s', feature_id)
            self.update_node('KBaseGene', feature_id, feature_props)
        feature_node = self.get_node('KBaseGene', feature_id)
        #print(feature.functions_unsplit)
        annotation_nodes = self.add_annotation(feature.functions_unsplit, annotation_source)
        #print(annotation_nodes)

        for n in annotation_nodes:
            #print(feature_node.get('id'), n.get('key'))
            r1 = self.link_nodes(feature_node, n, 'has_annotation')
            r2 = self.link_nodes(n, feature_node, 'has_gene')
            if not 'function_source' in r1:
                self.set_relationship_property(r1.id, 'function_source', annotation_source)
            else:
                function_sources = set(r1['function_source'].split(';'))
                function_sources.add(annotation_source)
                self.set_relationship_property(r1.id, 'function_source', ';'.join(function_sources))

            if not 'function_source' in r2:
                self.set_relationship_property(r2.id, 'function_source', annotation_source)
            else:
                function_sources = set(r2['function_source'].split(';'))
                function_sources.add(annotation_source)
                self.set_relationship_property(r2.id, 'function_source', ';'.join(function_sources))

        return feature_node
    
    def add_kbase_genome(self, genome, genome_id, annotation_source):
        genome_node = self.get_genome(genome_id)
        if genome_node == None:
            logger.warning('ADD RefSeqGenome: %s', genome_id)
            props = process_dict(genome.data)
            genome_node = self.add_node('RefSeqGenome', genome_id, props)

        for f in genome.features:
            feature = KBaseGenomeFeature(f)
            feature_node = self.add_genome_feature(feature, genome_id, annotation_source)
            self.link_nodes(genome_node, feature_node, 'has_gene')
            
        return genome_node
            
    def list_genome_sets(self):
        m = self.neo4j_graph.nodes.match("GenomeSet")

        genome_sets = set()
        for n in m:
            genome_sets.add(n['key'])
        return genome_sets
    
    def get_genome_set(self, genome_set_id):
        m = self.neo4j_graph.nodes.match("GenomeSet", key=genome_set_id)
        if len(m) < 1:
            return None

        genomes = set()

        genome_set_node = m.first()
        for rel in self.neo4j_graph.match((genome_set_node, ), r_type="has_genome", ):
            genome_id = rel.end_node['key']
            genomes.add(genome_id)

        return genomes

    def add_genome_to_genome_set(self, genome_set_id, genome_id):
        genome_set_node = self.get_node('GenomeSet', genome_set_id)
        genome_node = self.get_node('RefSeqGenome', genome_id)
        if genome_set_node == None or genome_node == None:
            return False
        self.link_nodes(genome_set_node, genome_node, 'has_genome')
        return True
        
    def get_ko_genes(self, ko_id):
        gene_nodes = None
        m = self.neo4j_graph.nodes.match("KeggOrthology", key=ko_id)
        if len(m) < 1:
            return res
        ko_node = m.first()
        gene_nodes = set()
        #KeggOrthology -[:has_gene]-> Node
        for rel in self.neo4j_graph.match((ko_node, ), r_type="has_gene", ):
            n = rel.end_node
            #gene_data = dict(n)
            #gene_data['database_id'] = n.identity
            #print(n.identity, type(n.identity))
            gene_nodes.add(n.identity)
        return gene_nodes
        
    def get_functional_roles2(self, ko_id):
        m = self.neo4j_graph.nodes.match("KeggOrthology", key=ko_id)
        if len(m) < 1:
            pass
        ko_node = m.first()
        gene_nodes = set()
        #KeggOrthology -[:has_gene]-> Node
        for rel in self.neo4j_graph.match((ko_node, ), r_type="has_gene", ):
            gene_nodes.add(rel.end_node)

        functions = {}
        functions_data = {}

        found = 0
        for gene_node in gene_nodes:
            #print(gene_node)
            gene_functions = set()
            #gene_node -[:has_annotatation]-> gene_function
            for rel in self.neo4j_graph.match((gene_node, ), r_type="has_annotation", ):
                function = Neo4jAnnotationFunction(rel.end_node)
                gene_functions.add(function)
            if not len(gene_functions) == 0:
                found+= 1
                for gene_function in gene_functions:
                    if not gene_function.id in functions:
                        functions[gene_function.id] = set()
                    functions[gene_function.id].add(gene_node['key'])
                    functions_data[gene_function.id] = {'value' : gene_function.value}
                    for sub_function in gene_function.sub_functions:
                        if not sub_function.id in functions:
                            functions[sub_function.id] = set()
                        functions[sub_function.id].add(gene_node['key'])
                        functions_data[sub_function.id] = {'value' : sub_function.value}
                    #print(gene_node['id'], gene_function['key'])

        return functions, functions_data, {'total' : len(gene_nodes), 'has_function' : found}
    
    @staticmethod
    def _get_template_reaction_data(tx, template_id, rxn_id):
        query = 'MATCH(n1:TemplateSet {key:$template_id})-[r1:has_annotation_rule]->(n2:TemplateReactionAnnotation)-[r2:has_reaction]->(n3:ModelSeedReaction {key:$rxn_id}) RETURN n1, n2, n3'
        result = tx.run(query, template_id=template_id, rxn_id=rxn_id)
        return list(result.records())

    def get_template_reaction_data(self, template_id, rxn_id):
        result = None
        with self.driver.session() as session:
            result = session.read_transaction(self._get_template_reaction_data, template_id, rxn_id)
        return result
    
    @staticmethod
    def filter_genome_set(fcount, genome_set):
        if genome_set == None:
            return fcount
        fcount_filter = {}
        for source_id in fcount:
            fcount_filter[source_id] = {
                'genes' : set(),
                'genomes' : fcount[source_id]['genomes'] & genome_set,
            }
            for gene_id in fcount[source_id]['genes']:
                a, b = gene_id.split('@')
                if b in genome_set:
                    fcount_filter[source_id]['genes'].add(gene_id)
            #for 'genomes' in fcount[]
        return fcount_filter
    
    def get_ko_function_data(self, kos, manual_ko):
        logger.debug("%s::get_ko_function_data", self)
        function_data = {}
        for ko in kos:
            if ko in manual_ko and not manual_ko[ko]:
                print('manual_excluded', ko)
            else:
                functions, functions_data, metadata = self.get_functional_roles2(ko)
                logger.warning('get_ko_function_data: %s:%s', ko, metadata)
                for f in functions:
                    function = functions_data[f]['value']
                    if not function in function_data:
                        function_data[function] = {
                            'id' : f,
                            'hits' : []
                        }
                    function_data[function]['hits'].append({
                        'score' : len(functions[f]),
                        'source' : ['KEGG', ko]
                    })
        return function_data
    
    def get_ko_function_data3(self, kos):
        logger.debug("%s::get_ko_function_data3", self)
        function_data = {}
        for ko in kos:
            functions, functions_data, metadata = self.get_functional_roles2(ko)
            logger.debug('get_ko_function_data: %s:%s', ko, metadata)
            for f in functions:
                function = functions_data[f]['value']
                if not function in function_data:
                    function_data[function] = {
                        'id' : f,
                        'hits' : []
                    }
                function_data[function]['hits'].append({
                    'score' : len(functions[f]),
                    'source' : ['KEGG', ko]
                })
        return function_data
    
    def get_reaction_annotation_data3_2(self, rxn_id, kos, 
                                        genome_set_id = None, genome_set = None, example_genes = 10):
        function_data = self.get_ko_function_data3(kos)
        
        #MISSING ADD TEMPLATE DATA
        for template_set in self.page_nodes('TemplateSet', 0, 10):
            template_id = template_set['n']['key']
            res = self.get_template_reaction_data(template_id, rxn_id)
            for r in res:
                #print('pair', r['n2']['key'])
                node = self.neo4j_graph.nodes[r['n2'].id]
                for rel in self.neo4j_graph.match((node, ), r_type="has_function", ):
                    if rel.end_node.has_label('FunctionComplex'):
                        for rel_cpx in self.neo4j_graph.match((rel.end_node, ), r_type="has_function", ):
                            function = Neo4jAnnotationFunction(rel_cpx.end_node)
                            if not function.value in function_data:
                                function_data[function.value] = {
                                    'id' : function.id,
                                    'hits' : []
                                }
                            function_data[function.value]['hits'].append({
                                'score' : 0,
                                'source' : ['template', template_id]
                            })
                    else:
                        function = Neo4jAnnotationFunction(rel.end_node)
                        if not function.value in function_data:
                            function_data[function.value] = {
                                'id' : function.id,
                                'hits' : []
                            }
                        function_data[function.value]['hits'].append({
                            'score' : 0,
                            'source' : ['template', template_id]
                        })
                        
        #MISSING ADD SBML DATA

        for f in function_data:
            node = self.matcher.match("Function", key=f).first()
            annotation = Neo4jAnnotationFunction(node)

            subsystem_tags = {}
            for ss in annotation.function_group:
                subsystem_tags[ss] = {}
            
            logger.debug('get fcount %s', f)
            fcount = self.get_function_count(annotation)
            logger.debug('filter fcount genome set')
            fcount = self.filter_genome_set(fcount, genome_set)
            source_tags = {}
            for s in fcount:
                source_tags[s] = [len(fcount[s]['genomes']), 
                                  len(fcount[s]['genes']), 
                                  list(fcount[s]['genes'])[:example_genes]]

            function_data[f]['sources'] = source_tags
            function_data[f]['subsystems'] = subsystem_tags

        return function_data
    
    def get_reaction_annotation_data3(self, rxn_id, genome_set_id = None, example_genes = 10, 
                                      manual_ko = {}, manual_function = {}):
        
        genome_set = None
        if not genome_set_id == None:
            genome_set = self.get_genome_set(genome_set_id)
        
        selected_kos = set()
        
        rxn_kos = self.get_ko_by_seed_id(rxn_id)
        if rxn_kos == None:
            rxn_kos = set()
        selected_kos |= rxn_kos
        
        for ko in manual_ko:
            if manual_ko[ko]:
                o = self.get_node('KeggOrthology', ko)
                if not o == None:
                    selected_kos.add(ko)
            elif ko in selected_kos:
                selected_kos.remove(ko)
        
        function_data = self.get_reaction_annotation_data3_2(rxn_id, selected_kos, genome_set_id, genome_set, example_genes)
        return function_data
        
    
        
    def get_reaction_annotation_data(self, rxn_id, genome_set = None, example_genes = 10, manual_ko = {}, manual_function = {}):
        kos = self.get_ko_by_seed_id(rxn_id)
        
        #print('get_reaction_annotation_data', rxn_id, len(genome_set), example_genes, manual_ko, manual_function)
        #print('get_reaction_annotation_data', kos)
        
        if kos == None:
            kos = set()
        
        for ko in manual_ko:
            if manual_ko[ko]:
                o = self.get_node('KeggOrthology', ko)
                if not o == None:
                    kos.add(ko)

        function_data = self.get_ko_function_data(kos, manual_ko)

        #MISSING ADD TEMPLATE DATA
        for template_set in self.page_nodes('TemplateSet', 0, 10):
            template_id = template_set['n']['key']
            res = self.get_template_reaction_data(template_id, rxn_id)
            for r in res:
                #print('pair', r['n2']['key'])
                node = self.neo4j_graph.nodes[r['n2'].id]
                for rel in self.neo4j_graph.match((node, ), r_type="has_function", ):
                    if rel.end_node.has_label('FunctionComplex'):
                        for rel_cpx in self.neo4j_graph.match((rel.end_node, ), r_type="has_function", ):
                            function = Neo4jAnnotationFunction(rel_cpx.end_node)
                            if not function.value in function_data:
                                function_data[function.value] = {
                                    'id' : function.id,
                                    'hits' : []
                                }
                            function_data[function.value]['hits'].append({
                                'score' : 0,
                                'source' : ['template', template_id]
                            })
                    else:
                        function = Neo4jAnnotationFunction(rel.end_node)
                        if not function.value in function_data:
                            function_data[function.value] = {
                                'id' : function.id,
                                'hits' : []
                            }
                        function_data[function.value]['hits'].append({
                            'score' : 0,
                            'source' : ['template', template_id]
                        })
        #MISSING ADD SBML DATA

        for f in function_data:
            node = self.matcher.match("Function", key=f).first()
            annotation = Neo4jAnnotationFunction(node)

            subsystem_tags = {}
            for ss in annotation.function_group:
                subsystem_tags[ss] = {}
            
            logger.debug('get fcount %s', f)
            fcount = self.get_function_count(annotation)
            logger.debug('filter fcount genome set')
            fcount = self.filter_genome_set(fcount, genome_set)
            source_tags = {}
            for s in fcount:
                source_tags[s] = [len(fcount[s]['genomes']), 
                                  len(fcount[s]['genes']), 
                                  list(fcount[s]['genes'])[:example_genes]]

            function_data[f]['sources'] = source_tags
            function_data[f]['subsystems'] = subsystem_tags

        return function_data
    
    def get_function_count(self, annotation):
        #function = annotation_api.get_annotation_function('Allophanate hydrolase 2 subunit 2 (EC 3.5.1.54)')
        #print(function)
        #genes with this function
        result = {}
        #genes = {}
        matches = self.fetch_gene_annotation_relationships(annotation.value)
        for m in matches:
            gene_id, genome_id = m['gene_id'].split('@')
            sources = set(m['function_source'].split(';')) if type(m['function_source']) == str else set()
            for s in sources:
                if not s in result:
                    result[s] = {
                        'genes' : set(),
                        'genomes' : set(),
                    }

                result[s]['genes'].add(m['gene_id'])
                result[s]['genomes'].add(genome_id)

        #for srel in annotation.node.graph.match((None, annotation.node), r_type="has_annotation", ):
            #sources = set() if not 'function_source' in srel else set(srel['function_source'].split(';'))
            #gene_node = srel.start_node
            #gene_id, genome_id = gene_node['key'].split('@')

            #gene_node = srel.start_node
            #genes[gene_node] = sources

        #for g in genes:
            #gene_id, genome_id = g['key'].split('@')
            #sources = genes[g]
            #print(gene_id, genome_id, g['function'], g['protein_translation_length'], sources)
        #owners of the annotation
        #sources = function.neo4j_load_source()
        #print(sources)
        #count owners X function
        #propagate to genomes
        #genomes = set()
        #for gene_node in genes:
        #    for srel in gene_node.graph.match((None, gene_node), r_type="has_gene", ):
        #        if srel.start_node.has_label("RefSeqGenome"):
        #            genome_node = srel.start_node
        #            #print(genome_node['key'])
        #count owners X genomes
        #get Y examples of genes for each owners
        return result
            
