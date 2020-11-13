import logging
import pymongo
import cobrakbase
import redis
%run ../../export_template_to_kbase.py
%run ../../../annotation-server/annotation_api_neo4j.py

logger = logging.getLogger(__name__)

kbase = cobrakbase.KBaseAPI()

modelseed_local = cobrakbase.modelseed.from_local('../../../../ModelSEEDDatabase')

host, port, user, pwd = ("0.0.0.0", 7687, "neo4j", "123585")
def init_annotation_api(host, port, user, pwd):
    annotation_api = AnnotationApiNeo4j(user=user, pwd=pwd, port=port, host=host)
    annotation_api.neo4j_graph = Graph("http://neo4j:" + pwd + "@" + host + ":7474")
    annotation_api.matcher = NodeMatcher(annotation_api.neo4j_graph)
    annotation_api.r_matcher = RelationshipMatcher(annotation_api.neo4j_graph)
    annotation_api.init_constraints()
    return annotation_api

annotation_api = init_annotation_api(host, port, user, pwd)

aclient = pymongo.MongoClient("mongodb+srv://server:dx75S3HBXX6h2U3D@bios-dk66o.gcp.mongodb.net/test?retryWrites=true&w=majority")
#aclient = pymongo.MongoClient('mongodb://192.168.1.15:27017/')
database = aclient['annotation']

logger = logging.getLogger('cobrakbase.core.kbasefba.template_curation')
logger.setLevel(logging.ERROR)
logger