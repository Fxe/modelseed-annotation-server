import logging
import copy
import re
from html.parser import HTMLParser

logger = logging.getLogger(__name__)

class PaperBlastParser(HTMLParser):
    
    def __init__(self):
        super().__init__()
        self.hits = None
        self.parse_hits = False
        self.hit_start_p = True
        self.read_gene = False
        self.read_gene_field = None
        self.read_li = False
        self.data = {}
        self.cursor = None
        self.ul_level = 0
        self.a_level = 0
        self.tag_depth = {}
        self.cursor_db = None
        self.cursor_url = None
        self.abstract = None
        self.literature_record = None
        self.last_open_tag = None
        self.rank = 0
        
        self.database_url_reader = {
            'SwissProt' : lambda x : x.split('/')[-1],
            'RefSeq' : lambda x : x.split('/')[-1],
            'MicrobesOnline' : lambda x : x.split('=')[-1],
            'BRENDA' : lambda x : x.split('=')[-1]
        }
        
    def handle_starttag(self, tag, attrs):
        if not tag in self.tag_depth:
            self.tag_depth[tag] = 0
        self.tag_depth[tag] += 1
        if tag == 'ul':
            self.ul_level += 1
        if tag == 'a':
            self.a_level += 1
        if tag == 'a' and 'small' in self.tag_depth and self.tag_depth['small'] > 0:
            #print(tag, self.literature_record, self.tag_depth['small'], attrs)
            for t in attrs:
                if t[0] == 'title' and not self.literature_record == None:
                    self.literature_record['authors'] = t[1]
            
        if tag == 'li' and self.parse_hits and not self.cursor == None:
            #print('ul', self.tag_depth['ul'])
            self.abstract = ""
        
        #Detect start of a literature event
        if tag == 'a' and 'ul' in self.tag_depth and self.literature_record == None:
            for t in attrs:
                if t[0] == 'href' and 'http://www.ncbi.nlm.nih.gov/pmc/articles' in t[1]:
                    logger.debug('BEGIN literature')
                    self.literature_record = {
                                    'url' : None,
                                    'text' : None,
                                    'title' : None,
                                    'journal' : None,
                                    'authors' : None
                                }
                    self.literature_record['url'] = t[1]
            

            #print(self.tag_depth['ul'], attrs)
        #if tag == 'li' and self.parse_hits and not self.cursor == None:
        #    print("read pub:", tag, self.tag_depth['ul'])
            
        if tag == 'a' and self.parse_hits and self.cursor == None:
            for t in attrs:
                if t[0] == 'title':
                    self.cursor_db = t[1]
                if t[0] == 'href':
                    self.cursor_url = t[1]
            
        if self.parse_hits and self.hit_start_p and tag == 'p':
            #print("Hit start:", tag)
            pass
        if self.parse_hits and self.hit_start_p and tag == 'ul':
            #print("read pub:", tag, self.ul_level)
            self.read_li = True
            self.hit_start_p = False
            
        if self.parse_hits and self.hit_start_p and tag == 'a':
            self.read_gene = True
            #print("a:", tag, attrs)
        if self.parse_hits and tag == 'h3':
            self.parse_hits = False
            
        self.last_open_tag = tag
        
    def handle_endtag(self, tag):
        if not tag in self.tag_depth:
            self.tag_depth[tag] = 0
        self.tag_depth[tag] -= 1
        
        if tag == 'ul' and not self.abstract == None and not self.literature_record == None:
            self.literature_record['text'] = self.abstract
            self.data[self.cursor]['papers'].append(copy.deepcopy(self.literature_record))
            self.literature_record = None
            self.abstract == None
            logger.debug('END literature')
            
            
        if tag == 'ul':
            self.ul_level -= 1
        if tag == 'a':
            self.a_level -= 1
        if self.parse_hits and not self.hit_start_p and tag == 'ul' and self.ul_level == 0:
            #print("stop pub:", tag, self.ul_level)
            self.cursor = None
            self.hit_start_p = True

    def parse_hit_count(self, text):
        o = text.split('Found')[1].split('similar')[0]
        return int(o)
    
    def handle_data(self, data):
        if self.last_open_tag == 'i' and not self.cursor == None:
            if self.data[self.cursor]['organism'] == None:
                self.data[self.cursor]['organism'] = data
        if 'a' in self.tag_depth and self.tag_depth['a'] > 0 and not self.literature_record == None:
            if self.literature_record['title'] == None:
                self.literature_record['title'] = data
                logger.debug('*** literature %s %s', self.tag_depth['a'], data)
        if 'small' in self.tag_depth and self.tag_depth['small'] > 0 and not self.literature_record == None:
            self.literature_record['journal'] = data.strip()
            #print(self.tag_depth['small'], data, self.literature_record)
        #if 'small' in self.tag_depth and self.tag_depth['small'] > 0 and not self.literature_record == None:
        #    self.literature_record['journal'] = data
            
        if 'li' in self.tag_depth and not self.abstract == None:
            self.abstract += data.replace('\n', '')
            
        if self.read_gene:
            if self.cursor == None:
                self.cursor = data
                self.data[data] = {
                    'rank' : self.rank,
                    'organism' : None,
                    'papers' : [],
                    'url' : self.cursor_url,
                    'db' : self.cursor_db,
                    'gene_id' : None if not self.cursor_db in self.database_url_reader else self.database_url_reader[self.cursor_db](self.cursor_url)
                }
                self.rank += 1
            #print('gene!', data)
            self.read_gene = False
        
        if '% identity' in data and '% coverage' in data:
            coverage = re.findall(r'\d+', data)
            identity, coverage = list(map(lambda x : float(x) / 100, coverage))
            self.data[self.cursor]['identity'] = identity
            self.data[self.cursor]['coverage'] = coverage
        
        if 'similar proteins in the literature:' in data:
            self.hits = self.parse_hit_count(data)
            self.parse_hits = True
            
