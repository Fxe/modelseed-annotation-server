import logging
import json
import jsonpickle
from annotation_api_neo4j import AnnotationApiNeo4j
from neo4j import GraphDatabase
from py2neo import Graph, NodeMatcher, RelationshipMatcher

logger = logging.getLogger(__name__)


class AnnotationApiRedisCache(AnnotationApiNeo4j):
  
    def __init__(self, cache, driver=None, user=None, pwd=None, uri=None, port=7474, host="0.0.0.0"):
        self.driver = driver
        self.neo4j_graph = Graph("http://neo4j:" + pwd + "@" + host + ":7474")
        self.matcher = NodeMatcher(self.neo4j_graph)
        self.r_matcher = RelationshipMatcher(self.neo4j_graph)
        if self.driver == None:
            self.driver = GraphDatabase.driver("bolt://" + host + ":" + str(port), auth=(user, pwd))
        if self.neo4j_graph == None:
            self.neo4j_graph = Graph("http://" + host + ":" + str(port), auth=(user, pwd))
            self.matcher = NodeMatcher(self.neo4j_graph)
        self.cache = cache
        
    def set_cache(self, cache_key, data, expire_sec = 86400):
        logger.debug('set_cache [%s]', cache_key)
        if data == None or len(data) == 0:
            return False
        self.cache.set(cache_key, jsonpickle.encode(data))
        return self.cache.expire(cache_key, expire_sec)
    
    def get_cache(self, cache_key):
        logger.debug('get_cache [%s]', cache_key)
        data = self.cache.get(cache_key)
        if not data:
            return None
        return jsonpickle.decode(data.decode('utf-8'))
        
    def get_cache_get_function_count(self, gene_function):
        logger.debug('get_cache_get_function_count [%s]', gene_function.id)
        cache_key = 'neo4j:' + str(gene_function.id)
        data = self.cache.hgetall(cache_key)
        if len(data) == 0:
            return None

        fcount = {}
        for kpair in data:
            #print(str(kpair))
            annotation_source, sub_key = kpair.decode('utf-8').split('#')
            if not annotation_source in fcount:
                fcount[annotation_source] = {}
            fcount[annotation_source][sub_key] = set(data[kpair].decode('utf-8').split('#'))
        return fcount

    def set_cache_get_function_count(self, gene_function, data, expire_sec = 86400):
        if data == None or len(data) == 0:
            return False
        cache_key = 'neo4j:' + str(gene_function.id)
        fcount_wrap = {}
        for annotation_source in data:
            for sub_key in data[annotation_source]:
                fcount_wrap[annotation_source + '#' + sub_key] = '#'.join(data[annotation_source][sub_key])
        self.cache.hmset(cache_key, fcount_wrap)
        return self.cache.expire(cache_key, expire_sec)
    
    def get_function_count(self, annotation):
        logger.debug('get_function_count::cache lookup')
        res = self.get_cache_get_function_count(annotation)
        if res:
            logger.debug('cache return')
            return res
        
        logger.debug('cache not found get live')
        res = super().get_function_count(annotation)
        
        logger.debug('cache live data')
        self.set_cache_get_function_count(annotation, res)
        return res
    
    def get_cache_get_genome_set(self, genome_set_id):
        cache_key = 'genome_set:' + genome_set_id
        data = self.cache.get(cache_key)
        if not data:
            return None

        return set(data.decode('utf-8').split('#'))

    def set_cache_get_genome_set(self, genome_set_id, data, expire_sec = 86400):
        if data == None or len(data) == 0:
            return False
        cache_key = 'genome_set:' + genome_set_id
        self.cache.set(cache_key, '#'.join(data))
        return self.cache.expire(cache_key, expire_sec)
    
    def get_genome_set(self, genome_set_id):
        if genome_set_id == None:
            return None
        logger.debug('get_genome_set::cache lookup')
        res = self.get_cache_get_genome_set(genome_set_id)
        if res:
            logger.debug('get_genome_set::cache return')
            return res
        
        logger.debug('get_genome_set::cache not found get live')
        res = super().get_genome_set(genome_set_id)
        
        logger.debug('get_genome_set::cache live data')
        self.set_cache_get_genome_set(genome_set_id, res)
        return res
            
    def get_ko_function_data(self, kos, manual_ko):
        logger.debug('get_ko_function_data::cache lookup')
        cache_key = 'ko_function_data::' + '#'.join(sorted(kos))
        res = self.get_cache(cache_key)
        if res:
            logger.debug('get_ko_function_data::cache return')
            return res
        
        logger.debug('get_ko_function_data::cache not found get live')
        res = super().get_ko_function_data(kos, manual_ko)
        
        logger.debug('get_ko_function_data::cache live data')
        self.set_cache(cache_key, res)
        return res
    
    def get_functional_roles2(self, ko):
        logger.debug('get_functional_roles2::cache lookup')
        cache_key = 'get_functional_roles2$' + ko
        res = self.get_cache(cache_key)
        if res:
            logger.debug('get_functional_roles2::cache return')
            return res
        
        logger.debug('get_functional_roles2::cache not found get live')
        res = super().get_functional_roles2(ko)
        
        logger.debug('get_functional_roles2::cache live data')
        self.set_cache(cache_key, res)
        return res
    
    def get_reaction_annotation_data3_2(self, rxn_id, kos, 
                                        genome_set_id = None, genome_set = None, example_genes = 10):
        logger.debug('get_reaction_annotation_data3_2::cache lookup')
        genome_set_id_str = 'None' if genome_set_id == None else genome_set_id
        cache_key = 'get_reaction_annotation_data3_2${}:{}:{}:{}'.format(rxn_id, '#'.join(kos), genome_set_id_str, example_genes)
        res = self.get_cache(cache_key)
        if res:
            logger.debug('get_reaction_annotation_data3_2::cache return')
            return res
        
        logger.debug('get_reaction_annotation_data3_2::cache not found get live')
        res = super().get_reaction_annotation_data3_2(rxn_id, kos, 
                                                      genome_set_id, genome_set, example_genes)
        
        logger.debug('get_reaction_annotation_data3_2::cache live data')
        self.set_cache(cache_key, res)
        return res