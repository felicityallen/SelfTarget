import http
import os

import bravado
import bravado.exception
import pytest

from server.request_types import FORECasTRequest

EXAMPLE_SEQUENCE = "ATGCTAGCTAGGGCATGAGGCATGCTAGTGACTGCATGGTAC"
EXAMPLE_PAM_IDX = 17
fr = FORECasTRequest(EXAMPLE_SEQUENCE, EXAMPLE_PAM_IDX)


def remove_example_cached_results():
    try:
        os.remove(fr.filename)
    except OSError:
        pass


def test_ping(client):
    assert client.ping.ping().response().metadata.status_code == http.HTTPStatus.OK


def test_get_profile_no_seq_and_pam_existing(client):
    remove_example_cached_results()
    with pytest.raises(bravado.exception.HTTPNotFound):
        r = client.profile.downloadProfile(seq=fr.seq, pam_idx=fr.pam_idx).response()
        assert r.metadata.status_code == http.HTTPStatus.NOT_FOUND


def test_get_plot_for_seq_and_pam(client):
    r = client.plot.searchBySeqAndPAM(body={"seq": fr.seq, "pam_idx": fr.pam_idx}).response()
    assert r.metadata.status_code == http.HTTPStatus.OK
    assert "mpld3.draw_figure" in r.result.get("plot")


def test_get_profile_with_seq_and_pam_existing(client):
    """
    depends on running test_get_plot_for_seq_and_pam first
    """
    r = client.profile.downloadProfile(seq=fr.seq, pam_idx=fr.pam_idx).response()
    assert r.metadata.status_code == http.HTTPStatus.OK
    assert r.metadata.headers['Content-Disposition'].startswith("attachment")


def test_get_plot_from_wge_id(client):
    r = client.plot.getPrecomputedPlot(wge="901413985", species="human").response()
    assert r.metadata.status_code == http.HTTPStatus.OK
    assert "mpld3.draw_figure" in r.result.get("plot")


def test_get_plot_from_oligo_id(client):
    r = client.plot.getPrecomputedPlot(wge="901413985", species="human").response()
    assert r.metadata.status_code == http.HTTPStatus.OK
    assert "mpld3.draw_figure" in r.result.get("plot")


def test_get_precomputed_plot_missing_value(client):
    with pytest.raises(bravado.exception.HTTPNotFound):
        r = client.plot.getPrecomputedPlot(wge="22901413985", species="human").response()
    with pytest.raises(bravado.exception.HTTPBadRequest):
        r = client.plot.getPrecomputedPlot(species="human").response()


