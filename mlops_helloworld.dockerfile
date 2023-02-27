FROM python:3.8

RUN pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple && \
    pip install -U pydflow && \
    pip install -U dp-tracking-sdk dp-metadata-sdk -i https://repo.mlops.dp.tech/repository/pypi-group/simple

COPY . /data/mlops-hello-world
RUN cd /data/mlops-hello-world && python setup.py install