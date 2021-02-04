import http
import os
from abc import abstractmethod
from functools import wraps
from typing import Dict

from flask import jsonify

RequestType = Dict[str, str]


class ServerException(Exception):

    def __init__(self, message, status_code):
        self.message = message
        self.status_code = status_code

    def send(self):
        return jsonify({'error': self.message}), self.status_code


def server_exceptions(f):
    @wraps(f)
    def new_function():
        try:
            return f()
        except ServerException as se:
            return se.send()
    return new_function


class SelfTargetRequest:

    @classmethod
    @abstractmethod
    def create_from_data(cls, data: RequestType) -> __qualname__:
        pass

    @classmethod
    def get_object_or_fail(cls, data: RequestType) -> __qualname__:
        try:
            return cls.create_from_data(data)
        except ValueError as e:
            raise ServerException(str(e), http.HTTPStatus.BAD_REQUEST)


class FORECasTRequest(SelfTargetRequest):

    def __init__(self, seq: str, pam_idx: str):
        if not(seq and pam_idx):
            raise ValueError("Target sequence or pam index not provided")
        self.seq = seq
        try:
            self.pam_idx = int(pam_idx)
        except ValueError:
            raise ValueError("PAM index must be a positive integer")
        self.filename = self._get_precomputed_file_path()

    def _get_precomputed_file_path(self) -> str:
        return os.path.join(os.getenv("PRECOMPUTED_PLOTS_DIR", ""), f'{self.seq}_{self.pam_idx}.txt')

    @classmethod
    def create_from_data(cls, data: dict):
        seq = data.get("seq", "")
        pam_idx = data.get("pam_idx", "")
        return cls(seq, pam_idx)


class PrecomputedProfileRequest(SelfTargetRequest):
    def __init__(self, wge: str, oligoid: str, species: str):
        if not ((wge or oligoid) and species):
            raise ValueError("Species and one of WGE or Oligo ID is required")
        self.wge = wge
        self.oligoid = oligoid
        self.species = species

    @classmethod
    def create_from_data(cls, data: RequestType):
        wge = data.get("wge", "")
        oligoid = data.get("id", "")
        species = data.get("species", "")
        return cls(wge, oligoid, species)
