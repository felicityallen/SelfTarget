#!/usr/bin/env python
import logging
import os
import re
import subprocess
import sys
from typing import Dict, List

import requests

from scripts.wge_to_ccid.constants import WGE_API_URL, HUMAN, MOUSE, NEGATIVE, POSITIVE
from scripts.wge_to_ccid.helper import get_file_name_without_extension, is_sequence

logging.basicConfig(
    level=logging.INFO,
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("result.log")
    ]
)
logger = logging.getLogger(__name__)

WGE_PAM_RIGHT_NEGATIVE = 0
WGE_PAM_RIGHT_POSITIVE = 1
WGE_PAM_RIGHT_BOTH = 2
STRANDS = {
    WGE_PAM_RIGHT_NEGATIVE: NEGATIVE,
    WGE_PAM_RIGHT_POSITIVE: POSITIVE
}


HUMAN_INDEX = os.getenv("HUMAN_INDEX",
                        "/lustre/scratch117/cellgen/cellgeni/cellgeni_su/crispr_indexes/GRCh38_index.bin")
MOUSE_INDEX = os.getenv("MOUSE_INDEX",
                        "/lustre/scratch117/cellgen/cellgeni/cellgeni_su/crispr_indexes/GRCm38_index.bin")
INDEXES = {
    HUMAN: HUMAN_INDEX,
    MOUSE: MOUSE_INDEX
}


Wge_dict = Dict[str, List[int]]
Crispr_line_string = List[str]


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
        if is_sequence(line):
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
