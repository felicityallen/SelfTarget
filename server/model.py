import os
from urllib.parse import urljoin

from flask import Flask
from flask_cors import CORS
from mongoengine import connect, Document, StringField, MultipleObjectsReturned
from selftarget.profile import CrisprLine
import os
import tempfile
from typing import Tuple

import requests
from selftarget.profile import readSummaryToProfile, get_guide_info_from_oligo_id, CrisprLine
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
        return obj


    @staticmethod
    def get_obj_by_id(id, species: str) -> 'WGE':
        if WGE._is_wge_id(id):
            objects = WGE.objects(wge_id=id, species=species)
        else:
            objects = WGE.objects(oligo_id=id, species=species)
        return WGE._get_object(objects, species)

    @staticmethod
    def read_profile(obj: 'WGE') -> Tuple[str, str, dict, CrisprLine]:
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


class FORECasTRequest:

    def __init__(self, seq, pam_idx):
        if not(seq and pam_idx):
            raise ValueError("Target sequence or pam index not provided")
        self.seq = seq
        try:
            self.pam_idx = int(pam_idx)
        except ValueError:
            raise ValueError("PAM index must be a positive integer")

    def get_precomputed_file_path(self):
        return os.path.join(os.getenv("PRECOMPUTED_PLOTS_DIR", ""), f'{self.seq}_{self.pam_idx}.txt')

    @classmethod
    def create_from_data(cls, data):
        seq = data.get("seq", "")
        pam_idx = data.get("pam_idx", "")
        return cls(seq, pam_idx)


class NoWGEException(Exception):

    def __init__(self, species, seq):
        self.species = species
        self.seq = seq

    def msg(self):
        return f"Error - no WGE found for seq {self.seq} {self.species}"
