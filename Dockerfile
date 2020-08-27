FROM python

RUN pip install --upgrade pip
RUN pip install requests
RUN pip install flask
RUN pip install flask_restful
RUN pip install redis
RUN pip install pymongo
RUN pip install py2neo==4.3.0
RUN pip install dnspython
RUN pip install networkx
RUN pip install neo4j==1.7.6
RUN pip install cobra
#RUN pip install cobrakbase==0.2.7
RUN pip install Escher==1.6.0
RUN pip install flask-cors
RUN pip install pyyaml
RUN pip install jsonpickle

RUN mkdir -p /opt/data
RUN git clone https://github.com/ModelSEED/ModelSEEDDatabase.git /opt/data/ModelSEEDDatabase

RUN mkdir -p /opt/build
RUN git clone https://github.com/ModelSEED/modelseed-escher.git /opt/build/modelseed-escher

RUN mkdir -p /opt/build
RUN git clone https://github.com/Fxe/biosapi.git /opt/build/biosapi

RUN mkdir -p /opt/build
RUN git clone https://github.com/Fxe/cobrakbase.git /opt/build/cobrakbase

RUN pip install /opt/build/cobrakbase
RUN pip install /opt/build/modelseed-escher
RUN pip install /opt/build/biosapi


COPY . /opt/annotation
COPY entrypoint.sh /
#ENV my_env_var=

#ENTRYPOINT python /opt/annotation/server.py


ENTRYPOINT ["/entrypoint.sh"]