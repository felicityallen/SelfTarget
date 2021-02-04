import os

HUMAN = "human"
MOUSE = "mouse"


def is_comprised_of_integers(string: str) -> bool:
    # or condition necessary so that 0 is treated as an integer
    return int(string) or True
