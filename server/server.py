import http
import logging
import os
from typing import Union, Any, Tuple

import mpld3
from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
from indel_prediction.predictor.predict import plot_predictions

from server.model import NoWGEException, WGE
from server.request_types import FORECasTRequest, PrecomputedProfileRequest, ServerException, server_exceptions

ResponseType = Tuple[Any, http.HTTPStatus]
model_path = "indel_prediction/predictor/model_output_10000_0.01000000_0.01000000_-0.607_theta.txt_cf0.txt"

app = Flask(__name__)
CORS(app)


@app.route('/ping', methods=['GET'])
def ping() -> ResponseType:
    return "The server is alive!", http.HTTPStatus.OK


def send_forecast_file(fr) -> ResponseType:
    if os.path.exists(fr.filename):
        return send_file(fr.filename, as_attachment=True)
    else:
        raise ServerException('Profile with those target sequence and pam index not found', http.HTTPStatus.NOT_FOUND)


@app.route('/api/profile', methods=['GET'])
@server_exceptions
def get_profile_file_for_seq_and_pam_idx() -> ResponseType:
    data = request.args or request.get_json()
    fr = FORECasTRequest.get_object_or_fail(data)
    return send_forecast_file(fr)


def get_wge_or_send_error(pr: PrecomputedProfileRequest) -> Union[WGE, ResponseType]:
    try:
        return WGE.get_obj_by_id(pr.wge or pr.oligoid, pr.species)
    except NoWGEException as ex:
        raise ServerException(ex.msg(), http.HTTPStatus.NOT_FOUND)


def get_graph_html_or_send_error(figure) -> Union[str, ResponseType]:
    try:
        figure.set_size_inches(11, 5.2)
        return mpld3.fig_to_html(figure,
                                 template_type="simple",
                                 figid="plot",
                                 no_extras=True)
    except Exception as e:
        logging.exception("Model error")
        raise ServerException(str(e), http.HTTPStatus.INTERNAL_SERVER_ERROR)


def get_precomputed_plot() -> ResponseType:
    data = request.args
    pr = PrecomputedProfileRequest.get_object_or_fail(data)
    wge: WGE = get_wge_or_send_error(pr)
    figure = wge.get_figure()
    graph_html = get_graph_html_or_send_error(figure)
    return jsonify({"plot": graph_html,
                    "guide": wge.get_guide().to_dict()})


def generate_plot() -> ResponseType:
    data = request.form or request.get_json()
    fr = FORECasTRequest.get_object_or_fail(data)
    figure = plot_predictions(model_path, fr.seq, fr.pam_idx)
    graph_html = get_graph_html_or_send_error(figure)
    return jsonify({"plot": graph_html})


@app.route('/plot', methods=['GET', 'POST'])
@server_exceptions
def plot() -> ResponseType:
    if request.method == 'GET':
        return get_precomputed_plot()
    elif request.method == 'POST':
        return generate_plot()
