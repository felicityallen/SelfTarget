import logging
import os

import mpld3
from flask import Flask, request, jsonify, send_file
from indel_prediction.predictor.predict import plot_predictions as main

app = Flask(__name__)

model_path = "indel_prediction/predictor/model_output_10000_0.01000000_0.01000000_-0.607_theta.txt_cf0.txt"


@app.route('/api/profile', methods=['GET'])
def get_profile():
    """
    This endpoint returns a profile as a text file. It is created and saved automatically if the plot
    has been built before (via the `plot` endpoint) for a sequence and PAM index pair.
    If the plot has not been built, there is no profile file saved, so an error is raised.

    request body: {"seq": "SEQUENCE", "pam_idx": "NUMBER"}
    :return: {"plot": "plot data"}
    :return: {"error": "error message"}
    """
    data = request.args or request.get_json()
    seq = data.get("seq", "")
    pam_idx = int(data.get("pam_idx", ""))
    filename = '{0}_{1}.txt'.format(seq, pam_idx)
    if os.path.exists(filename):
        return send_file(filename, as_attachment=True)
    else:
        return jsonify({'error': 'Profile with those target sequence and pam index not found'}), 404


@app.route('/plot', methods=['POST'])
def plot():
    """
    This endpoint returns a javascript code that is rendered in a browser as a profile plot
    for a particular pair of sequence and PAM index.

    request body: {"seq": "SEQUENCE", "pam_idx": "NUMBER"}
    :return: {"plot": "plot data"}
    :return: {"error": "error message"}
    """
    data = request.form or request.get_json()
    seq = data.get("seq", "")
    pam_idx = int(data.get("pam_idx", ""))
    if not (seq and pam_idx):
        return jsonify({'message': 'Empty request'}), 400
    try:
        graph_html = mpld3.fig_to_html(main(model_path, seq, pam_idx),
                                       template_type="simple",
                                       figid="plot",
                                       no_extras=True)
    except Exception as e:
        logging.exception("Model error")
        return jsonify({"error": str(e)})

    return jsonify({"plot": graph_html})


if __name__ == "__main__":
    app.run(port=5001)
