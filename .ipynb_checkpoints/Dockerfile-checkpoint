FROM python

RUN pip install --upgrade pip
RUN pip install requests
RUN pip install flask
RUN pip install flask_restful
RUN pip install pymongo
RUN pip install neo4j
RUN pip install cobra
RUN pip install cobrakbase

COPY . /opt/annotation

ENTRYPOINT python /opt/annotation/server.py