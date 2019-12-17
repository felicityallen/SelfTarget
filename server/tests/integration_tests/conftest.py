import os

import pytest
import yaml
from bravado.client import SwaggerClient

HOST = os.getenv("SELF_TARGET_SERVER_HOST", "http://127.0.0.1:5003")


@pytest.fixture(scope='module')
def client():
    with open("../../swagger.yml") as f:
        spec = yaml.load(f)

    _client = SwaggerClient.from_spec(spec)
    _client.swagger_spec.api_url = HOST
    return _client
