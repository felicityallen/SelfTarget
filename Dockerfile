FROM tiangolo/uwsgi-nginx-flask:python3.6

COPY indel_analysis /app/indel_analysis

RUN apt update && apt install -y cmake && \
    cd /app/indel_analysis/indelmap &&  \
    cmake . -DINDELMAP_OUTPUT_DIR=/usr/local/bin && \
    make && make install

COPY requirements.txt /app

RUN pip install --upgrade pip && pip install -r /app/requirements.txt

COPY . /app

RUN cd /app/indel_prediction && pip install . && \
    cd /app/selftarget_pyutils && pip install .


ENV INDELGENTARGET_EXE /usr/local/bin/indelgentarget
ENV LISTEN_PORT 8006
ENV PYTHONPATH=/

EXPOSE 8006

