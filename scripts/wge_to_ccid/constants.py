import os

WGE_API_HOST = os.getenv("WGE_API_HOST", "https://www.sanger.ac.uk/htgt/wge/api/")
SEARCH_BY_SEQ_URL = os.path.join(WGE_API_HOST, "search_by_seq")  # ?pam_right=2&species=GRCh38&seq=GTTCGCCTTGCGCCATGGAC
CRISPR_BY_ID_URL = os.path.join(WGE_API_HOST, "crispr_by_id")  # ?species=GRCm38&id=507425671f
NEGATIVE = "-"
POSITIVE = "+"

HUMAN = "human"
MOUSE = "mouse"
