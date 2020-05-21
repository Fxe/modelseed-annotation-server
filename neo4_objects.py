class Neo4jNode:
    
    def __init__(self, uid, properties = {}):
        self._uid = uid
        self._properties = properties
        
    @property
    def uid(self):
        return self._uid
    
    @property
    def properties(self):
        return self._properties

class Neo4jGenome(Neo4jNode):
    
    def __init__(self, node):
        properties = {}
        for k in node.keys():
            properties[k] = node[k]
        super().__init__(node.identity, properties)
        self._genes = {}

        for r in node.graph.match((genome_node, ), r_type="has_gene"):
            self._genes[r.end_node.identity] = r.end_node
    
    @property
    def genes(self):
        return None
    
    @property
    def get_gene_by_id(self, gene_id):
        return None
    
    @property
    def get_gene_by_uid(self, uid):
        return None
    
    @property
    def get_gene_by_function(self, function):
        return None