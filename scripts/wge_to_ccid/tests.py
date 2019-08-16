from .wge_to_ccid import CrisprLine


def test_is_seq():
    assert CrisprLine.is_sequence("ATCTTCAGTCTCTGGAATTG")
    assert CrisprLine.is_sequence("A")
    assert not CrisprLine.is_sequence("ATCTTCAGTCTCTGGAATTG:")
    assert not CrisprLine.is_sequence("ATCTTCAGKCTCTGGAATTG")
    assert not CrisprLine.is_sequence("")
