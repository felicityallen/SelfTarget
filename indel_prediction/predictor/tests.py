from selftarget.profile import CrisprLine


def test_crispr_line_create_coordinates_positive_strand():
    line = CrisprLine(["@@@CCDS100.2_chr1_9111372_+ TCAGCAGCATCCACAGGACA 23.34"])
    assert line.get_coordinates == "1:9111356-9111378"


def test_crispr_line_create_coordinates_negative_strand():
    line = CrisprLine(["@@@CCDS100.2_chr1_9104441_- GGAATCCCAAGGGACCCCAG 31.49"])
    assert line.get_coordinates == "1:9104436-9104458"


def test_crispr_line_parse_gene():
    line = CrisprLine(["@@@CCDS100.2_chr1_9104441_- GGAATCCCAAGGGACCCCAG 31.49"])
    assert line.get_gene == "CCDS100.2"
