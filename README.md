# SelfTarget
[![Docker Repository on Quay.io](https://quay.io/repository/coreos/quay-docs/status "Docker Repository on Quay.io")](https://quay.io/repository/felicityallen/selftarget)

Scripts for processing and predicting CRISPR/Cas9-generated mutations

## FORECasT Web server

To predict and view mutational profiles for individual gRNAs, please visit the FORECasT website at:

https://partslab.sanger.ac.uk/FORECasT

## Command line tool

Follow the installation instructions below, then from a command line:

```
cd indel_prediction
cd predictor
```

then...

#### FORECasT Single gRNA prediction

```
python FORECasT.py <target DNA sequence> <PAM index (0 based)> <output_file_prefix>
```

e.g. 

```
python FORECasT.py ATGCTAGCTAGGGCATGAGGCATGCTAGTGACTGCATGGTAC 17 test_output
```

Output will be in 

<output_file_prefix>_predictedindelsummary.txt

A list of predicted mutations, one per line, listed in order of decreasing predicted counts.
Each line contains an identifier string for the indel followed by a - (ignore this), and then a predicted read count (tab-delimited).

e.g. 
```
-	1000	(always 1000 reads - it is the original template sequence - here for viewer).
D2_L-3R0	550
I1_L-2C1R0	200
```

<output_file_prefix>_predictedreads.txt
A list of read sequences corresponding to each predicted mutation in the previous file. 
The format is read_id (ignore this), read sequence, mutation identifier (tab delimited), followed by a - (ignore this)

e.g.
```
ATGCTAGCTAGGGCATGAGGCATGCTAGTGACTGCATGGTAC	-	-
ATGCTAGCTAGGGCAAGGCATGCTAGTGACTGCATGGTAC	D2_L-3R0	-
ATGCTAGCTAGGGCATGGAGGCATGCTAGTGACTGCATGGTAC	I1_L-2C1R0	-
```

#### FORECasT Batch mode prediction

```
python FORECasT.py <batch_filename> <output_file_prefix>
```

where batch_filename is a tab-delimited file with columns:  ID, Target, PAM Index
e.g.
```
ID	Target	PAM Index
Guide_1	ATGCTAGCTAGGGCATGAGGCATGCTAGTGACTGCATGGTAC	17
Guide_2	ATCGATGACTGATCGTAGCTAGCTGGGATGCTAGCTAGTTGCATGCTAGGAGTCAGCTAG	23
Guide_3	GATAGTCGTAGGCTAGCTAGCTAGCTGGCAAGTGTGGAAAAGGGGATGCATGTA	26
```

Output will be in 
<output_file_prefix>_predictedindelsummary.txt  and
<output_file_prefix>_predictedreads.txt

which are formatted as for single mode, but separate guides are prefaced by a line with

```
@@@<ID>
```
where ID is the identifier provided for the guide in the batch file.

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



