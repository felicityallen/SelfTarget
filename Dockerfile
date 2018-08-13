FROM tiangolo/uwsgi-nginx-flask:python3.6

COPY . /app

RUN pip install --upgrade pip && pip install -r /app/requirements.txt

ENV LISTEN_PORT 8006
ENV PYTHONPATH=/

EXPOSE 8006

