FROM python

RUN pip install --upgrade pip
RUN pip install requests
RUN pip install flask
RUN pip install flask_restful
RUN pip install pymongo
RUN pip install py2neo
RUN pip install dnspython
RUN pip install neo4j
RUN pip install cobra
RUN pip install cobrakbase==0.2.2

COPY . /opt/annotation
COPY entrypoint.sh /
#ENV my_env_var=

#ENTRYPOINT python /opt/annotation/server.py

ENTRYPOINT ["/entrypoint.sh"]