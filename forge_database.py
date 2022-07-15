#!/usr/bin/env python
import os
import csv
from csv import DictWriter
import sys
import gemmi
from gemmi import cif, CifWalk, expand_if_pdb_code
import numpy as np
from numpy import genfromtxt
import subprocess

# mine directory for crystallographic information data
# this will break if there is one file that is not a cif file in the directory
def get_file_paths_from_args():
    for arg in sys.argv[1:]:
        if os.path.isdir(arg):
            for path in CifWalk(arg):
                yield path
        else:
            yield expand_if_pdb_code(arg)