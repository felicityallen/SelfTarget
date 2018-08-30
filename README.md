# SelfTarget
[![Docker Repository on Quay.io](https://quay.io/repository/coreos/quay-docs/status "Docker Repository on Quay.io")](https://quay.io/repository/felicityallen/selftarget)

Scripts for processing and predicting CRISPR/Cas9-generated mutations


## Installation

#### Locally
Create a Python 3 virtual environment and activate it
```bash
# install Python dependencies

pip install -r requirements.txt
cd selftarget_pyutils
pip install -e .
cd ../indel_prediction
pip install -e .

# compile predictor

cd indel_analysis/indelmap
cmake . -DINDELMAP_OUTPUT_DIR=/usr/local/bin
make && make install
export INDELGENTARGET_EXE /usr/local/bin/indelgentarget
```

#### Docker
Alternatively, you can start a Docker container and exec into it:
```bash
docker pull quay.io/felicityallen/selftarget
docker exec -it quay.io/felicityallen/selftarget bash
```
#### Web service

The predictor can be run as a web service. It can be accessed through a separate front end application 
[FORECasT](https://partslab.sanger.ac.uk) ([source on GitHub](https://github.com/cellgeni/FORECasT)). 
SelfTarget repository contains a Flask server with two API endpoints that are used by FORECasT to access predictor.

To run predictor as a server, you can follow the local installation steps above, 
go to the root directory and launch
```bash
python server.py --port=5001
```
or simply run a Docker container
```bash
docker run -d --name selftarget -p 5001:8006 quay.io/felicityallen/selftarget
```



