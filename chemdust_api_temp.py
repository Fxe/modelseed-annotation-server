import requests

class ChemDUST():
    
    def __init__(self, base_url = None):
        self.default_width = 100
        self.default_height = 100
        if base_url:
            self.base_url = base_url
        else:
            self.base_url = 'http://192.168.1.10:8066/ChemDUST'
        self.headers = {'Content-type': 'application/json', 'Accept': 'application/json'}
    
    def get_depict(self, structure, structure_format = 'smi', output_format = 'svg'):
        headers = {'Content-type': 'application/json', 'Accept': 'text/plain'}
        
        params = {
            "width" : self.default_width,
            "height" : self.default_height,
            "carbon_symbol" : False,
            "atom_color" : True,
            "atom_number" : False,
            "aromatic_display" : False,
            "structure" : structure
        }
        
        request_url = '{}/api/render/{}/{}'.format(self.base_url, 
                                                   output_format, 
                                                   structure_format)
        resp = requests.post(request_url, headers=headers, json=params)
        if resp.status_code != 200:
            #print(request_url, resp.status_code)
            raise ApiError('GET {} {}'.format(request_url, resp.status_code))
        return resp.content.decode('utf-8')
