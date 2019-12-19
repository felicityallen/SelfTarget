import os
import tempfile
from urllib.parse import urljoin

import requests
from flask import Flask
from flask_cors import CORS
from mongoengine import connect, Document, StringField, MultipleObjectsReturned
from predictor.predict import build_plot_by_profile
from selftarget.profile import readSummaryToProfile, get_guide_info_from_oligo_id

from server.helpers import is_comprised_of_integers

app = Flask(__name__)
CORS(app)

model_path = "indel_prediction/predictor/model_output_10000_0.01000000_0.01000000_-0.607_theta.txt_cf0.txt"

HUMAN = "human"
MOUSE = "mouse"

MONGODB_HOST = os.getenv("MONGODB_HOST")
S3_BASE = os.getenv("S3_BASE", "https://fa9.cog.sanger.ac.uk/")
DB_FILEPATH_BASE = os.getenv("DB_FILEPATH_BASE", "/lustre/scratch117/cellgen/team227/FORECasT_profiles_for_AK/")
WGE_INFO_URL_BASE = "https://www.sanger.ac.uk/htgt/wge/crispr/"
connect(host=MONGODB_HOST, connect=False)


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


class WGEData:

    def __init__(self, reads_file, profile_file, profile, crispr_line):
        self.reads_file = reads_file
        self.profile_file = profile_file
        self.profile = profile
        self.crispr_line = crispr_line

    def get_figure(self, oligo_id):
        return build_plot_by_profile(self.reads_file, self.profile, oligo_id)


class WGE(Document):
    wge_id = StringField(required=True, max_length=50, unique=True)
    oligo_id = StringField(max_length=50)
    filename = StringField(required=True, max_length=500)
    species = StringField(required=True, choices=[HUMAN, MOUSE], max_length=6)

    @staticmethod
    def _is_wge_id(id: str) -> bool:
        """
        if id is only comprised of integers, this is wge, if not, it's oligo_id
        """
        try:
            return is_comprised_of_integers(id)
        except ValueError:
            pass

    @staticmethod
    def _get_object(objects: 'WGE', species: str) -> 'WGE':
        if len(objects) > 1:
            raise MultipleObjectsReturned()
        elif len(objects) == 0:
            raise NoWGEException(species, id)
        obj = objects[0]
        obj['filename'] = os.path.join(S3_BASE, obj['filename'])
        obj.data = obj.get_data()
        return obj

    @staticmethod
    def get_obj_by_id(id, species: str) -> 'WGE':
        if WGE._is_wge_id(id):
            objects = WGE.objects(wge_id=id, species=species)
        else:
            objects = WGE.objects(oligo_id=id, species=species)
        return WGE._get_object(objects, species)

    def get_data(self) -> WGEData:
        reads_file = tempfile.mkstemp()[1]
        profile_file = tempfile.mkstemp()[1]
        r = requests.get(self.filename, allow_redirects=True)
        with open(reads_file, 'w') as f:
            f.write(r.text)
        r = requests.get(self.filename.replace("_predicted_rep_reads.txt", "_predicted_mapped_indel_summary.txt"),
                         allow_redirects=True)
        with open(profile_file, 'w') as f:
            f.write(r.text)
        crispr_line_info = get_guide_info_from_oligo_id(profile_file, self.oligo_id)
        profile = {}
        readSummaryToProfile(profile_file, profile, oligoid=self.oligo_id, remove_wt=False)
        return WGEData(reads_file, profile_file, profile, crispr_line_info)

    def get_guide(self) -> Guide:
        return Guide(wge_id=self.wge_id,
                     species=self.species,
                     gene=self.data.crispr_line.get_gene,
                     strand=self.data.crispr_line.get_strand,
                     coordinates=self.data.crispr_line.get_coordinates)

    def get_figure(self):
        return self.data.get_figure(self.oligo_id)


class NoWGEException(Exception):

    def __init__(self, species, seq):
        self.species = species
        self.seq = seq

    def msg(self):
        return f"Error - no WGE found for seq {self.seq} {self.species}"
