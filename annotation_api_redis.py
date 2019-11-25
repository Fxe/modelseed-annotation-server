import logging
from annotation_api_neo4j import AnnotationApiNeo4j
from neo4j import GraphDatabase
from py2neo import Graph, NodeMatcher, RelationshipMatcher

logger = logging.getLogger(__name__)

class AnnotationApiRedisCache(AnnotationApiNeo4j):
  
    def __init__(self, cache, driver=None, user=None, pwd=None, uri=None, port=7474, host="0.0.0.0"):
        self.driver = driver
        self.neo4j_graph = None
        self.matcher = NodeMatcher(self.neo4j_graph)
        if self.driver == None:
            self.driver = GraphDatabase.driver("bolt://" + host + ":" + str(port), auth=(user, pwd))
        if self.neo4j_graph == None:
            self.neo4j_graph = Graph("http://" + host + ":" + str(port), auth=(user, pwd))
            self.matcher = NodeMatcher(self.neo4j_graph)
        self.cache = cache