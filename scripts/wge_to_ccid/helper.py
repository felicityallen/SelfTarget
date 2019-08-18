import os


def check_if_abs(filename) -> str:
    if not os.path.isabs(filename):
        filename = os.path.abspath(filename)
    return filename


def get_file_name_without_extension(filename: str) -> str:
    EXTENSIONS = [".txt"]
    basename = os.path.basename(filename)
    for ext in EXTENSIONS:
        if basename.endswith(ext):
            return os.path.splitext(basename)[0]
    return basename
