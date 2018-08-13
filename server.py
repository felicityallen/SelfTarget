import mpld3
from flask import Flask, request, jsonify

from compare_overbeek_profiles import main

app = Flask(__name__)


@app.route('/plot', methods=['POST'])
def plot():
    """
    request body: {"seq": "SEQUENCE", "cut_site": "NUMBER"}
    :return: {"answer": answer}
    """
    data = request.form or request.get_json()
    seq = data.get("seq", "")
    cut_site = int(data.get("cut_site", ""))
    if not (seq and cut_site):
        return jsonify({'message': 'Empty request'}), 400
    graph_html = mpld3.fig_to_html(main(seq, cut_site),
                                   template_type="simple",
                                   figid="plot")


    return jsonify({"plot": graph_html})


if __name__ == "__main__":
    app.run(port=5001)
