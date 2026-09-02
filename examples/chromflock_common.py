import os
import sys
import pathlib
import shutil
import numpy as np

srcdir = pathlib.Path(__file__).resolve().parent.parent

chr_sizes = np.array([247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754])

def find_mflock():
    """ Return the path to mflock, with preference over the standard locations
    where it might be built in this repo"""
    alternatives = [srcdir.joinpath('build/mflock'), srcdir.joinpath('bin/mflock'), srcdir.joinpath('build/mflock')]
    for name in alternatives:
        if os.path.exists(name):
            return str(name)
    return shutil.which('mflock')
