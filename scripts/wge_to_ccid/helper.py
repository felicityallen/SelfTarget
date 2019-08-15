import os


def check_if_abs(filename) -> str:
    if not os.path.isabs(filename):
        filename = os.path.abspath(filename)
    return filename


def strip_filesystem_specific_prefix(prefix, filename):
    return filename.replace(prefix, "")


def get_base_name(filename: str) -> str:
    return os.path.splitext(os.path.basename(filename))[0]
