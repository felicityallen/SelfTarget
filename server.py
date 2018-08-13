import mpld3
from flask import Flask, request, jsonify

from indel_prediction.predictor.predict import plot_predictions as main, setIndelGenTargetExeLoc

app = Flask(__name__)

model_path = "indel_prediction/predictor/model_output_2000_0.01000000_1.835_theta.txt_cf0.txt"

@app.route('/plot', methods=['POST'])
def plot():
    """
    request body: {"seq": "SEQUENCE", "pam_idx": "NUMBER"}
    :return: {"answer": answer}
    """
    data = request.form or request.get_json()
    seq = data.get("seq", "")
    pam_idx = int(data.get("pam_idx", ""))
    if not (seq and pam_idx):
        return jsonify({'message': 'Empty request'}), 400
    setIndelGenTargetExeLoc('indelgentarget')
    graph_html = mpld3.fig_to_html(main(model_path, seq, pam_idx),
                                   template_type="simple",
                                   figid="plot",
                                   no_extras=True)


    return jsonify({"plot": graph_html})


if __name__ == "__main__":
    app.run(port=5001)
    # seq = "CTGAGTAGCTATGCGGCCAGCAGCGAGACGCTCAGCGTGAAGCGGCAGTATCCCTCTTTCCTGCGCACCATCCCCAATC"
    # pam_idx = 42
    # graph_html = mpld3.fig_to_html(main(model_path, seq, pam_idx),
    #                                template_type="simple",
    #                                figid="plot")
