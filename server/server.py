import http
import logging
import os

import mpld3
from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
from indel_prediction.predictor.predict import plot_predictions
from predictor.predict import build_plot_by_profile
from werkzeug.exceptions import BadRequest

from server.model import NoWGEException, Guide, WGE, FORECasTRequest

app = Flask(__name__)
CORS(app)

model_path = "indel_prediction/predictor/model_output_10000_0.01000000_0.01000000_-0.607_theta.txt_cf0.txt"


@app.route('/ping', methods=['GET'])
def ping():
    return "The server is alive!", http.HTTPStatus.OK


@app.route('/api/profile', methods=['GET'])
def get_profile():
    """
    This endpoint returns a profile as a text file. It is created and saved automatically if the plot
    has been built before (via the `plot` endpoint) for a sequence and PAM index pair.
    If the plot has not been built, there is no profile file saved, so an error is raised.
    """
    data = request.args or request.get_json()
    try:
        forecast_request = FORECasTRequest.create_from_data(data)
    except ValueError as e:
        return jsonify({'error': str(e)}), http.HTTPStatus.BAD_REQUEST
    filename = forecast_request.get_precomputed_file_path()
    if os.path.exists(filename):
        return send_file(filename, as_attachment=True)
    else:
        return jsonify(
            {'error': 'Profile with those target sequence and pam index not found'}), http.HTTPStatus.NOT_FOUND


def get_precomputed_plot():
    data = request.args
    wge = data.get("wge", "")
    oligoid = data.get("id", "")
    species = data.get("species", "")
    if not ((wge or oligoid) and species):
        raise BadRequest()
    try:
        obj = WGE.get_obj_by_id(wge or oligoid, species)
    except NoWGEException as ex:
        return jsonify({"error": ex.msg()}), http.HTTPStatus.NOT_FOUND
    reads_file, profile_file, profile, crispr_line_info = WGE.read_profile(obj)
    guide = Guide.create_guide(obj, crispr_line_info)
    try:
        graph_html = mpld3.fig_to_html(build_plot_by_profile(reads_file, profile, obj['oligo_id']),
                                       template_type="simple",
                                       figid="plot",
                                       no_extras=True)
    except Exception as e:
        logging.exception("Model error")
        return jsonify({"error": str(e)}), http.HTTPStatus.INTERNAL_SERVER_ERROR
    return jsonify({"plot": graph_html,
                    "guide": guide.to_dict()})


def generate_plot():
    data = request.form or request.get_json()
    try:
        fr = FORECasTRequest.create_from_data(data)
    except ValueError as e:
        return jsonify({'error': str(e)}), http.HTTPStatus.BAD_REQUEST
    try:
        graph_html = mpld3.fig_to_html(plot_predictions(model_path, fr.seq, fr.pam_idx),
                                       template_type="simple",
                                       figid="plot",
                                       no_extras=True)
    except Exception as e:
        logging.exception("Model error")
        return jsonify({"error": str(e)}), http.HTTPStatus.INTERNAL_SERVER_ERROR

    return jsonify({"plot": graph_html})


@app.route('/plot', methods=['GET', 'POST'])
def plot():
    """
    This endpoint returns a javascript code that is rendered in a browser as a profile plot
    for a particular pair of sequence and PAM index.
    """
    if request.method == 'GET':
        return get_precomputed_plot()
    elif request.method == 'POST':
        return generate_plot()
