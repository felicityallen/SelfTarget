#!/usr/bin/env python
import argparse
import csv
import logging
import os
import sys
import traceback
from itertools import chain
from typing import Dict, List

import requests
from mongoengine import *

from scripts.wge_to_ccid.helper import check_if_abs, strip_filesystem_specific_prefix, get_base_name

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
INDEXES = {
    HUMAN: HUMAN_INDEX,
    MOUSE: MOUSE_INDEX
}

MONGODB_HOST = os.getenv("MONGODB_HOST", 'mongodb://cellgeni:cellgeni@172.27.82.34:32556/wge')

Wge_dict = Dict[str, List[int]]
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

    @classmethod
    def create_from_crispr(cls, wge_id, crispr):
        filename = strip_filesystem_specific_prefix(ROOT_PATH, crispr.filename)
        cls(wge_id=str(wge_id), oligo_id=crispr.oligo_id, filename=filename, species=crispr.species)


class Parser:
    pass


class APIParser(Parser):
    pass


class AnalyserParser(Parser):
    pass


class Crispr:

    def __init__(self, strand, species, seq, oligo_id, filename):
        self.strand = strand
        self.species = species
        self.seq = seq
        self.oligo_id = oligo_id
        self.filename = filename

    def validate_entry(self) -> None:
        assert self.strand in {'-', '+'}

    def extract_wge(self, mode: str, wge_dict=None) -> int:
        self.validate_entry()
        if mode == API:
            wge = self.match_wge_id_from_api()
        elif mode == Analyser:
            wge = self.match_wge_id_from_analyser(wge_dict)
        return wge

    def extract_and_save_wge(self, mode=str, wge_dict=None) -> int:
        wge = self.extract_wge(mode, wge_dict)
        logger.info(f"WGE: {wge}")
        WGE.create_from_crispr(wge, self).save()
        return wge

    def select_between_two_wge_ids(self, wge_ids):
        if self.strand == '-':
            return wge_ids[0]
        else:
            return wge_ids[1]

    def get_wge_id_from_list(self, wge_ids) -> int:
        if len(wge_ids) == 0:
            raise NoWGEException(self.species, self.seq)
        elif len(wge_ids) == 1:
            return wge_ids[0]
        elif len(wge_ids) == 2:
            return self.select_between_two_wge_ids(wge_ids)
        else:
            raise ValueError("More than two WGE ids found for sequence")

    def match_wge_id_from_analyser(self, wges_dict) -> int:
        wges = wges_dict[self.seq]
        return self.get_wge_id_from_list(wges)

    def match_wge_id_from_api(self) -> int:
        wge_json = requests.get(SEQ_API_URL, {"pam_right": 2, "species": self.species, "seq": self.seq}).json()
        return self.get_wge_id_from_list(wge_json)


class Analyser:

    @staticmethod
    def run_analyser(samples_file: str, species, result_filename=None):
        analyser_path = os.getenv("CRISPR_ANALYSER", "crispr_analyser")
        assert os.path.exists(samples_file), samples_file
        assert os.path.exists(analyser_path), analyser_path
        index_file = INDEXES[species]
        assert os.path.exists(index_file), index_file
        import subprocess
        if not result_filename:
            result_filename = get_base_name(samples_file) + "_wges.txt"
        try:
            subprocess.check_output([f"{analyser_path}", "search", "-i", index_file,
                                     "-f", os.path.abspath(samples_file), ">", os.path.abspath(result_filename)],
                                    stderr=subprocess.STDOUT,
                                    shell=True)
        except subprocess.CalledProcessError as e:
            print(e.output)
        return result_filename

    @staticmethod
    def check_no_more_than_two_wge_ids_found(wge_dict: Wge_dict):
        assert all(len(value) <= 2 for value in wge_dict.values()), "Found more than two WGE ids for one of the crisprs"

    @staticmethod
    def extract_wge_ids(line: str, wge_dict: Wge_dict, current_seq: str):
        line = line.strip("\n")
        if CrisprLine(line).is_sequence():
            current_seq = line
            wge_dict[line] = []
        elif line.startswith("\t"):
            wge_dict[current_seq].append(int(line.strip("\t")))
        else:
            raise ValueError("Bad wge results file error)")
        return current_seq

    @classmethod
    def get_wges_dict(cls, wges_filename) -> Wge_dict:
        wge_dict = {}
        with open(wges_filename) as f:
            current_seq = ""
            for line in f:
                current_seq = cls.extract_wge_ids(line, wge_dict, current_seq)
        cls.check_no_more_than_two_wge_ids_found(wge_dict)


class CrisprLine:

    def __init__(self, line):
        self.line = line

    def is_crispr_line(self):
        return self.line[0][:3] == '@@@'

    def is_sequence(self):
        return set(self.line).issubset({"A", "C", "D", "G"})

    def get_sequence_info(self):
        s = self.line[0].split()
        oligo_id = s[0][3:]
        strand = oligo_id.split('_')[-1]
        seq = s[1]
        return oligo_id, strand, seq


class RepReadsFile:

    def __init__(self, filename):

        self.filename = self.validate_file(filename)
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
                line = CrisprLine(toks)
                if line.is_crispr_line():
                    n += 1
            return n

    def create_crispr_from_seq_line(self, line: CrisprLine) -> Crispr:
        oligo_id, strand, seq = line.get_sequence_info()
        return Crispr(strand, self.species, seq, oligo_id, self.filename)

    def map_ids_from_line(self, line):
        if line.is_crispr_line():
            crispr: Crispr = self.create_crispr_from_seq_line(line)
            try:
                wge = crispr.extract_and_save_wge()
            except NoWGEException as ex:
                logger.error(ex.msg())
            except NotUniqueError as ex:
                logger.warning(f"Skipping duplicate wge_id: {wge}")
            except Exception as e:
                traceback.print_exc()
                logger.error(f"A problem happened for {self.filename}, oligo_id {crispr.oligo_id}")

    def map_ids_from_api(self):
        with open(self.filename) as f:
            reader = csv.reader(f, delimiter='\t')
            for toks in reader:
                line = CrisprLine(toks)
                if line.is_crispr_line():
                    self.map_ids_from_line(line)

    def dump_seqs_into_txt(self):
        target_filename = get_base_name(self.filename) + "_seqs.txt"
        with open(self.filename) as f:
            with open(target_filename, 'w') as seqs_file:
                reader = csv.reader(f, delimiter='\t')
                for toks in reader:
                    line = CrisprLine(toks)
                    if line.is_crispr_line():
                        _, _, seq = line.get_sequence_info()
                        seqs_file.write(seq + '\n')
        return target_filename

    def parse_analyser_report(self, wges_dict: Wge_dict):
        with open(self.filename) as f:
            reader = csv.reader(f, delimiter='\t')
            for toks in reader:
                line = CrisprLine(toks)
                if line.is_crispr_line():
                    crispr = self.create_crispr_from_seq_line(line)
                    crispr.extract_and_save_wge(ANALYSER, wges_dict)

    def map_ids_from_analyzer(self):
        seqs_filename = self.dump_seqs_into_txt()
        wges_filename = Analyser.run_analyser(seqs_filename, self.species)
        wges_dict = Analyser.get_wges_dict(wges_filename)
        self.parse_analyser_report(wges_dict)

    def map_ids(self, mode: str):
        if mode == API:
            self.map_ids_from_api()
        elif mode == ANALYSER:
            self.map_ids_from_analyzer()
        else:
            raise ValueError(f"Unknown WGE extraction mode: {mode}")


class WgeToCCIDCLmanager:

    def __init__(self, args):
        self.args = args

    def process_single_file_mapping(self, args):
        repreads = RepReadsFile(args.file)
        if args.mode:
            repreads.map_ids(args.mode)

    def process_map_command(self):
        if self.args.file:
            self.process_single_file_mapping(self.args)

    def process_dump_command(self):
        repreads = RepReadsFile(self.args.file)
        repreads.dump_seqs_into_txt()

    def count_oligo_ids(self):
        skip_directories = {}
        human_files = list(map(lambda x: os.path.join(HUMAN_PATH, x), os.listdir(HUMAN_PATH)))
        mouse_files = list(map(lambda x: os.path.join(MOUSE_PATH, x), os.listdir(MOUSE_PATH)))
        total_oligo_ids = 0
        for d in chain(human_files, mouse_files):
            if "CCDS" in d and not d.endswith(".zip") and not d in skip_directories:
                for f in os.listdir(d):
                    if f.endswith("reads.txt"):
                        filename = os.path.join(d, f)
                        repreads = RepReadsFile(filename)
                        logger.info(repreads.filename)
                    try:
                        oligo_ids = repreads.count_oligo_ids_in_a_file()
                        logger.info(f"{filename}, oligo ids: {oligo_ids}")
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
    parser_dump = subparsers.add_parser(DUMP,
                                        help="Dump from a file into a txt file")
    parser_dump.add_argument("file", type=str, help="The file that contains crisprs to dump")
    args = parser.parse_args()
    manager = WgeToCCIDCLmanager(args)
    manager.process()


if __name__ == "__main__":
    main()
