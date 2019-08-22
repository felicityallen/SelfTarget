#!/usr/bin/env python
import csv
import logging
import os
import sys
import traceback
from typing import Dict, List

import requests
from mongoengine import connect, NotUniqueError, StringField, Document
from pymongo.errors import DuplicateKeyError

from scripts.wge_to_ccid.constants import NEGATIVE, POSITIVE, SEARCH_BY_SEQ_URL, HUMAN, MOUSE
from scripts.wge_to_ccid.helper import check_if_abs, get_file_name_without_extension
from scripts.wge_to_ccid.wge_retrievers import WGEObj, APIParser, Analyser

logging.basicConfig(
    level=logging.INFO,
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("result.log")
    ]
)
logger = logging.getLogger(__name__)

ANALYSER = "analyser"
API = "api"

ROOT_PATH = os.getenv("ROOT_PATH", "/lustre/scratch117/cellgen/team227/FORECasT_profiles_for_AK/")
HUMAN_PATH = os.path.join(ROOT_PATH, "human")
MOUSE_PATH = os.path.join(ROOT_PATH, "mouse")

WGE_POSITIVE_OFFSET = 16

MONGODB_HOST = os.getenv("MONGODB_HOST", 'mongodb://cellgeni:cellgeni@172.27.82.34:32556/wge')

Wge_dict = Dict[str, List[int]]
Crispr_line_string = List[str]
connect(host=MONGODB_HOST)


class NoWGEException(Exception):

    def __init__(self, species, seq):
        self.species = species
        self.seq = seq

    def msg(self):
        return f"Error - no WGE found for seq {self.species} {self.seq}"


class CrisprLine:

    def __init__(self, line: Crispr_line_string):
        """
        :param line: example - ["@@@CCDS103.1_chr1_9567325_+ AAGTCACTTTTTAGAGGCTT 31.24"]
        """
        assert self.is_crispr_line(line)
        s: List[str] = line[0].split()
        self._oligo_id = s[0][3:]
        self._seq = s[1]
        _, self._chromosome, self._location, self._strand = self._oligo_id.split('_')

    @property
    def get_oligo_id(self):
        return self._oligo_id

    @property
    def get_strand(self):
        return self._strand

    @property
    def get_seq(self):
        return self._seq

    @property
    def get_location(self):
        return self._location

    @property
    def get_chromosome(self):
        return int(self._chromosome.replace("chr", ""))

    @staticmethod
    def is_crispr_line(line: Crispr_line_string):
        return line[0][:3] == '@@@'


class Crispr:

    def __init__(self, crispr_line: CrisprLine, species: str, filename: str, database_filename: str):
        self.crispr_line = crispr_line
        self.species = species
        self.filename = filename
        self.database_filename = database_filename

    def validate_entry(self) -> None:
        assert self.crispr_line.get_strand in {NEGATIVE, POSITIVE}

    def extract_wge(self, mode: str, wge_dict=None) -> str:
        self.validate_entry()
        if mode == API:
            wge = self.match_wge_id_from_api()
        elif mode == ANALYSER:
            wge = self.match_wge_id_from_analyser(wge_dict)
        else:
            raise ValueError("Unknown mode")
        return wge

    def extract_and_save_wge(self, mode: str, wge_dict=None) -> str:
        wge = self.extract_wge(mode, wge_dict)
        logger.info(f"WGE: {wge}")
        wge2ccid = WGE.create_from_crispr(wge, self)
        wge2ccid.filename = self.database_filename
        wge2ccid.save()
        return wge

    def matches_wge(self, wge_obj: WGEObj):
        wge_obj_location = wge_obj.get_start_location()
        # We check whether objects match by the difference between their WGE start location and CRISPR location
        # which is constant. It is different for the positive strand and for the negative one
        # I haven't figured out what is it for the negative strand, so I'm just checking whether WGE start locations
        # is in the neighbourhood of the CRISPR location
        IN_NEIGHBOURHOOD = 50
        strand = self.crispr_line.get_strand
        if (strand == POSITIVE) and (strand == wge_obj.get_strand()) and \
                (abs(int(self.crispr_line.get_location) - wge_obj_location) == WGE_POSITIVE_OFFSET):
            return True
        elif (strand == NEGATIVE) and (strand == wge_obj.get_strand()) and \
                (abs(int(self.crispr_line.get_location) - wge_obj_location) <= IN_NEIGHBOURHOOD):
            return True
        else:
            return False

    def _select_between_multiple_wge_ids(self, wge_ids) -> str:
        for wge_id in wge_ids:
            wge = APIParser.get_wge_info(wge_id, self.species)
            if self.matches_wge(wge):
                return wge_id

    def _select_between_two_wge_ids(self, wge_ids) -> str:
        if self.crispr_line.get_strand == NEGATIVE:
            return wge_ids[0]
        else:
            return wge_ids[1]

    def get_wge_id_from_list(self, wge_ids) -> str:
        if len(wge_ids) == 0:
            raise NoWGEException(self.species, self.crispr_line.get_seq)
        elif len(wge_ids) == 1:
            return wge_ids[0]
        elif len(wge_ids) == 2:
            return self._select_between_two_wge_ids(wge_ids)
        else:
            return self._select_between_multiple_wge_ids(wge_ids)

    def match_wge_id_from_analyser(self, wges_dict) -> str:
        wges = wges_dict[self.crispr_line.get_seq]
        return self.get_wge_id_from_list(wges)

    def match_wge_id_from_api(self) -> str:
        wge_json = requests.get(SEARCH_BY_SEQ_URL,
                                {"pam_right": 2, "species": self.species, "seq": self.crispr_line.get_seq}).json()
        return self.get_wge_id_from_list(wge_json)


class WGE(Document):
    wge_id = StringField(required=True, max_length=50, unique=True)
    oligo_id = StringField(max_length=50)
    filename = StringField(required=True, max_length=500)
    species = StringField(required=True, choices=[HUMAN, MOUSE], max_length=6)

    @classmethod
    def create_from_crispr(cls, wge_id, crispr: Crispr):
        return cls(wge_id=str(wge_id), oligo_id=crispr.crispr_line.get_oligo_id, filename=crispr.filename,
                   species=crispr.species)


class RepReadsFile:

    def __init__(self, filename, database_name):

        self.filename = self.validate_file(filename)
        self.database_filename = database_name
        self.species = MOUSE if MOUSE in filename else HUMAN

    @staticmethod
    def validate_file(filename) -> str:
        filename = check_if_abs(filename)
        assert MOUSE in filename or HUMAN in filename
        return filename

    def count_oligo_ids_in_a_file(self) -> int:
        with open(self.filename) as f:
            reader = csv.reader(f, delimiter='\t')
            n = 0
            for toks in reader:
                if CrisprLine.is_crispr_line(toks):
                    n += 1
            return n

    def create_crispr_from_crispr_line(self, line: CrisprLine) -> Crispr:
        return Crispr(line, self.species, self.filename, self.database_filename)

    def map_ids_from_line(self, line: CrisprLine, mode, wges_dict=None):
        crispr: Crispr = self.create_crispr_from_crispr_line(line)
        try:
            crispr.extract_and_save_wge(mode, wges_dict)
        except NoWGEException as ex:
            logger.error(ex.msg())
            pass
        except (DuplicateKeyError, NotUniqueError) as ex:
            logger.warning("Skipping duplicate wge_id")
            pass
        except Exception as e:
            traceback.print_exc()
            logger.error(f"A problem happened for {self.filename}, oligo_id {crispr.crispr_line.get_oligo_id}")
            pass

    def map_ids_from_api(self):
        with open(self.filename) as f:
            reader = csv.reader(f, delimiter='\t')
            for toks in reader:
                if CrisprLine.is_crispr_line(toks):
                    self.map_ids_from_line(CrisprLine(toks), API)

    def _generate_seqs_file(self, target_filename):
        with open(self.filename) as f:
            with open(target_filename, 'w') as seqs_file:
                reader = csv.reader(f, delimiter='\t')
                for toks in reader:
                    if CrisprLine.is_crispr_line(toks):
                        line = CrisprLine(toks)
                        seqs_file.write(line.get_seq + '\n')
        return target_filename

    def dump_seqs_into_txt(self):
        target_filename = get_file_name_without_extension(self.filename) + "_seqs.txt"
        if not os.path.exists(target_filename):
            self._generate_seqs_file(target_filename)
        return target_filename

    def parse_analyser_report(self, wges_dict: Wge_dict):
        with open(self.filename) as f:
            reader = csv.reader(f, delimiter='\t')
            for toks in reader:
                if CrisprLine.is_crispr_line(toks):
                    self.map_ids_from_line(CrisprLine(toks), ANALYSER, wges_dict)

    def map_ids_from_analyzer(self, wges_file=None):
        if not wges_file:
            seqs_filename = self.dump_seqs_into_txt()
            analyser = Analyser(seqs_filename, self.species)
            wges_file = analyser.get_wge_ids_file()
        wges_dict = Analyser.get_wges_dict(wges_file)
        self.parse_analyser_report(wges_dict)

    def map_ids(self, mode: str, wge_file=None):
        if mode == API:
            self.map_ids_from_api()
        elif mode == ANALYSER:
            self.map_ids_from_analyzer(wge_file)
        else:
            raise ValueError(f"Unknown WGE extraction mode: {mode}")
