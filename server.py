import http
import logging
import os
import tempfile
from urllib.parse import urljoin

import mpld3
import requests
from flask import Flask, request, jsonify, send_file
from indel_prediction.predictor.predict import plot_predictions
from mongoengine import connect, Document, StringField, MultipleObjectsReturned
from predictor.predict import build_plot_by_profile
from selftarget.profile import readSummaryToProfile, get_guide_info_from_oligo_id, CrisprLine
from werkzeug.exceptions import BadRequest

app = Flask(__name__)

model_path = "indel_prediction/predictor/model_output_10000_0.01000000_0.01000000_-0.607_theta.txt_cf0.txt"

HUMAN = "human"
MOUSE = "mouse"

MONGODB_HOST = os.getenv("MONGODB_HOST")
S3_BASE = os.getenv("S3_BASE", "https://fa9.cog.sanger.ac.uk/")
DB_FILEPATH_BASE = os.getenv("DB_FILEPATH_BASE", "/lustre/scratch117/cellgen/team227/FORECasT_profiles_for_AK/")
WGE_INFO_URL_BASE = "https://www.sanger.ac.uk/htgt/wge/crispr/"
connect(host=MONGODB_HOST, connect=False)


class WGE(Document):
    wge_id = StringField(required=True, max_length=50, unique=True)
    oligo_id = StringField(max_length=50)
    filename = StringField(required=True, max_length=500)
    species = StringField(required=True, choices=[HUMAN, MOUSE], max_length=6)


class Guide:

    def __init__(self, wge_id, coordinates, strand, species, gene):
        self.wge_id = wge_id
        self.coordinates = coordinates
        self.strand = strand
        self.species = species
        self.gene = gene
        self.wge_link = self._get_wge_link()

    def to_dict(self):
        return {"wge_id": self.wge_id,
                "coordinates": self.coordinates,
                "strand": self.strand,
                "species": self.species,
                "gene": self.gene,
                "wge_link": self.wge_link}

    def _get_wge_link(self) -> str:
        return urljoin(WGE_INFO_URL_BASE, self.wge_id)

    @classmethod
    def create_guide(cls, wge_obj: WGE, crispr_line_info: CrisprLine):
        return cls(wge_id=wge_obj.wge_id,
                   species=wge_obj.species,
                   gene=crispr_line_info.get_gene,
                   strand=crispr_line_info.get_strand,
                   coordinates=crispr_line_info.get_coordinates)


class NoWGEException(Exception):

    def __init__(self, species, seq):
        self.species = species
        self.seq = seq

    def msg(self):
        return f"Error - no WGE found for seq {self.seq} {self.species}"


def get_obj_by_id(id, species):
    wge = False
    # if id is integer, this is wge, if not, it's oligo_id
    try:
        wge = int(id) or True
    except ValueError:
        pass
    if wge:
        objects = WGE.objects(wge_id=id, species=species)
    else:
        # TODO: fix oligo_ids, in this form it won't work
        objects = WGE.objects(oligo_id=id, species=species)
    if len(objects) > 1:
        raise MultipleObjectsReturned()
    elif len(objects) == 0:
        raise NoWGEException(species, id)
    obj = objects[0]
    obj['filename'] = os.path.join(S3_BASE, obj['filename'])
    return obj


def read_profile(obj):
    reads_file = tempfile.mkstemp()[1]
    profile_file = tempfile.mkstemp()[1]
    r = requests.get(obj['filename'], allow_redirects=True)
    with open(reads_file, 'w') as f:
        f.write(r.text)
    r = requests.get(obj['filename'].replace("_predicted_rep_reads.txt", "_predicted_mapped_indel_summary.txt"),
                     allow_redirects=True)
    with open(profile_file, 'w') as f:
        f.write(r.text)
    crispr_line_info = get_guide_info_from_oligo_id(profile_file, obj['oligo_id'])
    profile = {}
    readSummaryToProfile(profile_file, profile, oligoid=obj['oligo_id'], remove_wt=False)
    return reads_file, profile_file, profile, crispr_line_info


@app.route('/ping', methods=['GET'])
def ping():
    return "The server is alive!", http.HTTPStatus.OK


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
    pam_idx = data.get("pam_idx", "")
    if not (seq and pam_idx):
        return jsonify({'error': 'Target sequence or pam index not provided'}), http.HTTPStatus.BAD_REQUEST
    filename = '{0}_{1}.txt'.format(seq, pam_idx)
    if os.path.exists(filename):
        return send_file(filename, as_attachment=True)
    else:
        return jsonify(
            {'error': 'Profile with those target sequence and pam index not found'}), http.HTTPStatus.BAD_REQUEST


@app.route('/plot', methods=['GET', 'POST'])
def plot():
    """
    This endpoint returns a javascript code that is rendered in a browser as a profile plot
    for a particular pair of sequence and PAM index.

    request body: {"seq": "SEQUENCE", "pam_idx": "NUMBER"}
    :return: {"plot": "plot data"}
    :return: {"error": "error message"}
    """
    if request.method == 'GET':
        data = request.args
        wge = data.get("wge", "")
        oligoid = data.get("id", "")
        species = data.get("species", "")
        if not ((wge or oligoid) and species):
            raise BadRequest()
        try:
            obj = get_obj_by_id(wge or oligoid, species)
        except NoWGEException as ex:
            return jsonify({"error": ex.msg()}), http.HTTPStatus.NOT_FOUND
        reads_file, profile_file, profile, crispr_line_info = read_profile(obj)
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

    elif request.method == 'POST':
        data = request.form or request.get_json()
        seq = data.get("seq", "")
        pam_idx = int(data.get("pam_idx", ""))
        if not (seq and pam_idx):
            return jsonify({'message': 'Empty request'}), http.HTTPStatus.BAD_REQUEST
        try:
            graph_html = mpld3.fig_to_html(plot_predictions(model_path, seq, pam_idx),
                                           template_type="simple",
                                           figid="plot",
                                           no_extras=True)
        except Exception as e:
            logging.exception("Model error")
            return jsonify({"error": str(e)}), http.HTTPStatus.INTERNAL_SERVER_ERROR

        return jsonify({"plot": graph_html})
