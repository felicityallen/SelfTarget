#!/usr/bin/env python
import csv
import json
import os
import sys
import traceback
import logging
from itertools import chain

import requests
from mongoengine import *

logging.basicConfig(
    level=logging.INFO,
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("result.log")
    ]
)
logger = logging.getLogger(__name__)

SEQ_API_URL = "https://www.sanger.ac.uk/htgt/wge/api/search_by_seq"  # ?pam_right=2&species=Mouse&seq=GTTCGCCTTGCGCCATGGAC
WGE_API_URL = "https://www.sanger.ac.uk/htgt/wge/api/crispr_by_id"  # ?species=Grcm38&id=507425671f
PAM_RIGHT = 2
HUMAN = "human"
MOUSE = "mouse"


ROOT_PATH = os.getenv("ROOT_PATH", "/lustre/scratch117/cellgen/team227/FORECasT_profiles_for_AK/")

human_path = os.path.join(ROOT_PATH, "human")
mouse_path = os.path.join(ROOT_PATH, "mouse")

MONGODB_HOST = os.getenv("MONGODB_HOST", 'mongodb://cellgeni:cellgeni@172.27.82.34:32556/wge')

connect(host=MONGODB_HOST)


class NoWGEException(Exception):

    def __init__(self, species, seq):
        self.species = species
        self.seq = seq

    def msg(self):
        return f"Error - no WGE found for seq {self.species} {self.seq}"


class WGE(Document):
    wge_id = StringField(required=True, max_length=50, unique=True)
    oligo_id = StringField(max_length=50)
    filename = StringField(required=True, max_length=500)
    species = StringField(required=True, choices=[HUMAN, MOUSE], max_length=6)


def parse_from_wge(strand, species, seq):
    assert strand == '-' or strand == '+'
    wge_json = requests.get(SEQ_API_URL, {"pam_right": 2, "species": species, "seq": seq}).json()
    if len(wge_json) == 0:
        raise NoWGEException(species, seq)
    elif len(wge_json) == 1:
        return wge_json[0]
    else:
        logger.warning(f"Strand choice {strand} {species} {seq}")
        if strand == '-':
            return wge_json[0]
        else:
            return wge_json[1]


def parse_file(filename, wge_to_ccid):
    with open(filename) as f:
        reader = csv.reader(f, delimiter='\t')
        oligo_id = None
        i = 0
        for toks in reader:
            if toks[0][:3] == '@@@':
                s = toks[0].split()
                oligo_id = s[0][3:]
                strand = oligo_id.split('_')[-1]
                seq = s[1]
                assert MOUSE in filename or HUMAN in filename
                species = MOUSE if MOUSE in filename else HUMAN
                try:
                    wge = parse_from_wge(strand, species, seq)
                    logger.info(f"WGE: {wge}")
                    wge_to_ccid[wge] = {"filename": filename,
                                        'species': species,
                                        "oligo_id": oligo_id}
                    WGE(wge_id=str(wge), oligo_id=oligo_id, filename=filename, species=species).save()
                except NoWGEException as ex:
                    logger.error(ex.msg())
                except Exception as e:
                    traceback.print_exc()
                    logger.error(f"A problem happened for {filename}, oligo_id {oligo_id}")
                    continue


def main():
    skip_directories = {}
    wge_to_ccid = {}
    try:
        human_files = list(map(lambda x: os.path.join(human_path, x), os.listdir(human_path)))
        mouse_files = list(map(lambda x: os.path.join(mouse_path, x), os.listdir(mouse_path)))
        for d in chain(human_files, mouse_files):
            if "CCDS" in d and not d.endswith(".zip") and not d in skip_directories:
                for f in os.listdir(d):
                    if f.endswith("reads.txt"):
                        filename = os.path.join(d, f)
                        logger.info(filename)
                        parse_file(filename, wge_to_ccid)
    except Exception as e:
        traceback.print_exc()
    finally:
        with open("result.json", "a") as f:
            json.dump(wge_to_ccid, f)

# result = {}
# parse_file("/Users/ak27/programming/cellgenij/SelfTarget/scripts/mouse/CCDS22491.1_Tnpo2_predicted_rep_reads.txt", result)
# with open("result.json", "w") as f:
#     json.dump(result, f)
main()