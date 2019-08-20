#!/usr/bin/env python
import argparse
import logging
import os
import sys
import traceback
from itertools import chain
from typing import Dict

from scripts.wge_to_ccid.crispr import RepReadsFile

logging.basicConfig(
    level=logging.INFO,
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("result.log")
    ]
)
logger = logging.getLogger(__name__)

MAP = "map"
DUMP = "dump"

ROOT_PATH = os.getenv("ROOT_PATH", "/lustre/scratch117/cellgen/team227/FORECasT_profiles_for_AK/")
HUMAN_PATH = os.path.join(ROOT_PATH, "human")
MOUSE_PATH = os.path.join(ROOT_PATH, "mouse")

FILENAME_MAP = os.getenv("FILENAME_MAP", "file_map.txt")


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


class WgeToCCIDCLManager:

    def __init__(self, args):
        self.args = args
        database_name = FarmFileMap(FILENAME_MAP).get_database_filename(self.args.file)
        self.repreads = RepReadsFile(self.args.file, database_name)

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
