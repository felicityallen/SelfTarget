import http
import os

import bravado
import bravado.exception
import pytest

from server.server import get_precomputed_file_path

EXAMPLE_SEQUENCE = "ATGCTAGCTAGGGCATGAGGCATGCTAGTGACTGCATGGTAC"
EXAMPLE_PAM_IDX = 17


def remove_example_cached_results():
    try:
        os.remove(get_precomputed_file_path(EXAMPLE_SEQUENCE, EXAMPLE_PAM_IDX))
    except OSError:
        pass


def test_ping(client):
    assert client.ping.ping().response().metadata.status_code == http.HTTPStatus.OK


def test_get_profile_no_seq_and_pam_existing(client):
    remove_example_cached_results()
    with pytest.raises(bravado.exception.HTTPNotFound):
        r = client.profile.downloadProfile(seq=EXAMPLE_SEQUENCE, pam_idx=EXAMPLE_PAM_IDX).response()
        assert r.metadata.status_code == http.HTTPStatus.NOT_FOUND


def test_get_plot_for_seq_and_pam(client):
    r = client.plot.searchBySeqAndPAM(body={"seq": EXAMPLE_SEQUENCE, "pam_idx": EXAMPLE_PAM_IDX}).response()
    assert r.metadata.status_code == http.HTTPStatus.OK
    assert "mpld3.draw_figure" in r.result.get("plot")


# def test_get_profile_with_seq_and_pam_existing(client):
#     """
#     depends on running test_get_plot_for_seq_and_pam first
#     """
#     r = client.profile.downloadProfile(seq=EXAMPLE_SEQUENCE, pam_idx=EXAMPLE_PAM_IDX).response()
#     assert r.metadata.status_code == http.HTTPStatus.OK

