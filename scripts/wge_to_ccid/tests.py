from scripts.wge_to_ccid.helper import get_file_name_without_extension
from .wge_to_ccid import CrisprLine, FarmFileMap


def test_is_seq():
    assert CrisprLine.is_sequence("ATCTTCAGTCTCTGGAATTG")
    assert CrisprLine.is_sequence("A")
    assert not CrisprLine.is_sequence("ATCTTCAGTCTCTGGAATTG:")
    assert not CrisprLine.is_sequence("ATCTTCAGKCTCTGGAATTG")
    assert not CrisprLine.is_sequence("")


def test_filemap():
    filemap = FarmFileMap("test_data/filemap.txt")
    assert filemap.get_database_filename("/lustre/scratch117/cellgen/cellgeni/forecast/input.1") \
           == "human/CCDS0-499/CCDS100.2_GPR157_predicted_rep_reads.txt"
    assert not filemap.get_database_filename("/lustre/scratch117/cellgen/cellgeni/forecast/input.1") \
           == "/lustre/scratch117/cellgen/team227/FORECasT_profiles_for_AK/human/CCDS0-499/CCDS100.2_GPR157_predicted_rep_reads.txt"
    assert filemap.get_database_filename("input.1") == ""


def test_get_base_name():
    assert get_file_name_without_extension("input.1") == "input.1"
    assert get_file_name_without_extension("/lustre/scratch117/cellgen/cellgeni/forecast/human/input.1") == "input.1"
    assert get_file_name_without_extension("scratch117/cellgen/input.txt") == "input"

