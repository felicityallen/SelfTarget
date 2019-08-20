from scripts.wge_to_ccid.constants import HUMAN
from scripts.wge_to_ccid.crispr import CrisprLine, Crispr
from scripts.wge_to_ccid.helper import get_file_name_without_extension, is_sequence
from .wge_to_ccid import FarmFileMap


def test_is_seq():
    assert is_sequence("ATCTTCAGTCTCTGGAATTG")
    assert is_sequence("A")
    assert not is_sequence("ATCTTCAGTCTCTGGAATTG:")
    assert not is_sequence("ATCTTCAGKCTCTGGAATTG")
    assert not is_sequence("")


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


def test_cripr_multiple_matches_wge_positive():
    info = ["@@@CCDS103.1_chr1_9582370_+ GACGGCGCGCCTGGTGTTCC 26.99"]
    crispr = Crispr(CrisprLine(info), HUMAN, "/human/input.5", "some-filename")
    wge_ids = [	"901488483", "901488485", "1082205868", "1082205869"]
    assert crispr.get_wge_id_from_list(wge_ids) == "901488485"


def test_cripr_multiple_matches_wge_negative():
    info = ["@@@CCDS103.1_chr1_9582374_- GACGGCGCGCCTGGTGTTCC 26.99"]
    crispr = Crispr(CrisprLine(info), HUMAN, "/human/input.5", "some-filename")
    wge_ids = [	"901488483", "901488485", "1082205868", "1082205869"]
    assert crispr.get_wge_id_from_list(wge_ids) == "901488483"


def test_cripr_two_matches_wge_positive():
    info = ["@@@CCDS103.1_chr1_9582370_+ GACGGCGCGCCTGGTGTTCC 26.99"]
    crispr = Crispr(CrisprLine(info), HUMAN, "/human/input.5", "some-filename")
    wge_ids = [	"901488483", "901488485"]
    assert crispr.get_wge_id_from_list(wge_ids) == "901488485"
