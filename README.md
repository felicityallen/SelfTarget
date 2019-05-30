# SelfTarget
[![Docker Repository on Quay.io](https://quay.io/repository/coreos/quay-docs/status "Docker Repository on Quay.io")](https://quay.io/repository/felicityallen/selftarget)

Scripts for processing and predicting CRISPR/Cas9-generated mutations

## FORECasT Web server

To predict and view mutational profiles for individual gRNAs, please visit the FORECasT website at:

https://partslab.sanger.ac.uk/FORECasT

## Precomputed FORECasT Results for Human and Mouse CCDS

Precomputed profiles for all gRNAs in human and mouse CCDS regions are available here:

https://fa9.cog.sanger.ac.uk/index.html

Entries are collected into all gRNAs corresponding to each CCDS id. Within each file ending in _predicted_mapped_indel_summary.txt, the entries for each gRNA are separated by a line with @@@id guide_seq predicted_in_frame where the id contains the CCDS id, the chomosome coordinates and the strand. The next line is '- - 1000' and can be ignored (there for visualization only). The following lines are the particular indels predicted and their predicted counts (assuming total reads of 1000, and ignoring indels with less than 1 read). For the read sequences, see corrresponding entries in the _predicted_rep_reads.txt files.

## FORECasT Command line tool

1. Follow the installation instructions [here](#installation).

2. After installation, from a command line:
```
cd indel_prediction
cd predictor
```

3. Run single or batch prediction as described next.

#### Single gRNA prediction

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
-	-	1000	(always 1000 reads - it is the original template sequence - here for viewer use).
D2_L-3R0	-	550
I1_L-2C1R0	-	200
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

#### Batch mode prediction

```
python FORECasT.py <batch_filename> <output_file_prefix>
```

e.g.
```
python FORECasT.py example_batch.txt test_batch_output
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
export INDELGENTARGET_EXE=/usr/local/bin/indelgentarget
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



