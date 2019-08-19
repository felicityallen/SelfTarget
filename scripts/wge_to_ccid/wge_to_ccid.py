#!/usr/bin/env python
import argparse
import csv
import logging
import os
import subprocess
import re
import sys
import traceback
from itertools import chain
from typing import Dict, List

import requests
from mongoengine import connect, NotUniqueError, StringField, Document
from pymongo.errors import DuplicateKeyError

from scripts.wge_to_ccid.helper import check_if_abs, get_file_name_without_extension

logging.basicConfig(
    level=logging.INFO,
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("result.log")
    ]
)
logger = logging.getLogger(__name__)

SEQ_API_URL = "https://www.sanger.ac.uk/htgt/wge/api/search_by_seq"  # ?pam_right=2&species=GRCh38&seq=GTTCGCCTTGCGCCATGGAC
WGE_API_URL = "https://www.sanger.ac.uk/htgt/wge/api/crispr_by_id"  # ?species=GRCm38&id=507425671f
WGE_PAM_RIGHT_NEGATIVE = 0
WGE_PAM_RIGHT_POSITIVE = 1
WGE_PAM_RIGHT_BOTH = 2
NEGATIVE = "-"
POSITIVE = "+"
STRANDS = {
    WGE_PAM_RIGHT_NEGATIVE: NEGATIVE,
    WGE_PAM_RIGHT_POSITIVE: POSITIVE
}
HUMAN = "human"
MOUSE = "mouse"
ANALYSER = "analyser"
API = "api"
MAP = "map"
DUMP = "dump"

ROOT_PATH = os.getenv("ROOT_PATH", "/lustre/scratch117/cellgen/team227/FORECasT_profiles_for_AK/")
HUMAN_PATH = os.path.join(ROOT_PATH, "human")
MOUSE_PATH = os.path.join(ROOT_PATH, "mouse")
HUMAN_INDEX = os.getenv("HUMAN_INDEX",
                        "/lustre/scratch117/cellgen/cellgeni/cellgeni_su/crispr_indexes/GRCh38_index.bin")
MOUSE_INDEX = os.getenv("MOUSE_INDEX",
                        "/lustre/scratch117/cellgen/cellgeni/cellgeni_su/crispr_indexes/GRCm38_index.bin")
FILENAME_MAP = os.getenv("FILENAME_MAP", "file_map.txt")
INDEXES = {
    HUMAN: HUMAN_INDEX,
    MOUSE: MOUSE_INDEX
}

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


class WGEObj:

    def __init__(self, wge_obj: Dict):
        self._wge_id = list(wge_obj.keys())[0]
        self._wge_obj = wge_obj[self._wge_id]

    def get_start_location(self):
        return int(self._wge_obj.get("chr_start"))

    def get_strand(self):
        return STRANDS[self._wge_obj.get("pam_right")]


class Parser:
    # TODO: refactor parsing from API and from the analyser into two different classes
    pass


class APIParser(Parser):

    @staticmethod
    def get_wge_info(wge_id, species) -> WGEObj:
        r = requests.get(WGE_API_URL, {"species": species, "id": wge_id})
        wge = WGEObj(r.json())
        return wge


class AnalyserParser(Parser):
    pass


class FarmFileMap:
    """
    Acts upon a file with a format
    original_filename target_filename_that_was_used_for_bsub_submission

    Example:
    /lustre/human/CCDS0-499/CCDS430.1_RRAGC_predicted_rep_reads.txt /lustre/cellgeni/forecast/input.345
    """

    def __init__(self, filename):
        self.filemap = self._process_file_map(filename)

    @staticmethod
    def _process_file_map(filename) -> Dict[str, str]:
        filemap = {}
        with open(filename) as f:
            for line in f:
                original_file, farm_input_file = line.strip("\n").split(" ")
                filemap[farm_input_file] = original_file
        return filemap

    def get_database_filename(self, farm_input_filename) -> str:
        original_filename = self.filemap.get(farm_input_filename, "")
        database_name = original_filename.replace(ROOT_PATH, "")
        if database_name:
            logger.info(f"Database name: {database_name}")
        else:
            logger.error(f"Original filename not found for file: {farm_input_filename}")
        return database_name


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

    @staticmethod
    def is_sequence(line):
        line_characters = set(line)
        return line_characters and line_characters.issubset({"A", "C", "T", "G"})


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
        wge_json = requests.get(SEQ_API_URL,
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


class Analyser:

    def __init__(self, samples_file, species, result_filename=None):
        analyser_path = os.getenv("CRISPR_ANALYSER", "crispr_analyser")
        assert os.path.exists(samples_file), samples_file
        assert os.path.exists(analyser_path), analyser_path
        index_file = INDEXES[species]
        assert os.path.exists(index_file), index_file
        if not result_filename:
            result_filename = get_file_name_without_extension(samples_file) + "_wges.txt"

        self.samples_file = samples_file
        self.result_filename = result_filename
        self.analyser_path = analyser_path
        self.index_file = index_file

    def _run_search_command_line(self):
        command = " ".join([f"{self.analyser_path}", "search", "-i", self.index_file,
                            "-f", os.path.abspath(self.samples_file), ">", os.path.abspath(self.result_filename)])
        logger.info("Command: " + command)
        subprocess.check_output(command,
                                stderr=subprocess.STDOUT,
                                shell=True)

    def _generate_wge_ids_file(self):
        try:
            self._run_search_command_line()
        except subprocess.CalledProcessError as e:
            print(e.output)

    def get_wge_ids_file(self):
        if not os.path.exists(self.result_filename):
            self._generate_wge_ids_file()
        return self.result_filename

    @staticmethod
    def extract_wge_ids(line: str, wge_dict: Wge_dict, current_seq: str):
        line = line.strip("\n:")
        if CrisprLine.is_sequence(line):
            current_seq = line
            wge_dict[line] = []
        elif line.startswith("\t"):
            wge_dict[current_seq].append(int(line.strip("\t")))
        elif line == "":
            pass
        else:
            raise ValueError(f"Bad wge results file error, line is: {re.escape(line)}")
        return current_seq

    @classmethod
    def get_wges_dict(cls, wges_filename) -> Wge_dict:
        wge_dict = {}
        with open(wges_filename) as f:
            current_seq = ""
            for line in f:
                current_seq = cls.extract_wge_ids(line, wge_dict, current_seq)
        return wge_dict


class RepReadsFile:

    def __init__(self, filename, farm_file_map: FarmFileMap):

        self.filename = self.validate_file(filename)
        self.database_filename = farm_file_map.get_database_filename(os.path.abspath(self.filename))
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


class WgeToCCIDCLManager:

    def __init__(self, args):
        self.args = args
        self.repreads = RepReadsFile(self.args.file, FarmFileMap(FILENAME_MAP))

    def process_single_file_mapping(self):
        if self.args.mode:
            self.repreads.map_ids(self.args.mode, getattr(self.args, "wge"))

    def process_map_command(self):
        if self.args.file:
            self.process_single_file_mapping()

    def process_dump_command(self):
        self.repreads.dump_seqs_into_txt()

    def count_oligo_ids(self):
        skip_directories = {}
        human_files = list(map(lambda x: os.path.join(HUMAN_PATH, x), os.listdir(HUMAN_PATH)))
        mouse_files = list(map(lambda x: os.path.join(MOUSE_PATH, x), os.listdir(MOUSE_PATH)))
        total_oligo_ids = 0
        for d in chain(human_files, mouse_files):
            if not (not ("CCDS" in d) or d.endswith(".zip") or d in skip_directories):
                for f in os.listdir(d):
                    if f.endswith("reads.txt"):
                        # TODO: investigate why it's not used
                        filename = os.path.join(d, f)
                        logger.info(self.repreads.filename)
                    try:
                        oligo_ids = self.repreads.count_oligo_ids_in_a_file()
                        logger.info(f"{self.filename}, oligo ids: {oligo_ids}")
                        total_oligo_ids += oligo_ids
                    except Exception as e:
                        traceback.print_exc()
        logger.info(f"Total oligo ids: {total_oligo_ids}")

    def process(self):
        if self.args.command == MAP:
            self.process_map_command()
        elif self.args.command == DUMP:
            self.process_dump_command()


def main():
    help_text = """
            wge2ccid maps wge ids to oligo ids
            Example: $ cwl2argparse FILES [FILES ...] <options>
            """
    parser = argparse.ArgumentParser(description=help_text)
    subparsers = parser.add_subparsers(dest="command", help='sub-command help')
    parser_map = subparsers.add_parser(MAP,
                                       help="maps oligo ids. By default maps oligo ids from all human and mouse files")
    parser_map.add_argument('-f', '--file', type=str, help='Map oligo ids from a single file')
    parser_map.add_argument('-m', '--mode', type=str, default=ANALYSER,
                            help="Mode of mapping: either through WGE API or with CRISPR analyzer")
    parser_map.add_argument("-w", '--wge', type=str, help="File with sequences and corresponding WGE ids")
    parser_dump = subparsers.add_parser(DUMP,
                                        help="Dump from a file into a txt file")
    parser_dump.add_argument("file", type=str, help="The file that contains crisprs to dump")
    args = parser.parse_args()
    manager = WgeToCCIDCLManager(args)
    manager.process()


if __name__ == "__main__":
    main()
