#!/usr/bin/env python

""" Obtain the conformational entropy for every protein (cif file or gzipped cif file) within a given directory """
from __future__ import print_function
import sys
import os
from collections import Counter
from gemmi import cif, CifWalk, expand_if_pdb_code
import numpy as np
from numpy import genfromtxt


import sc_imports
from obtain_monomer_coordinates import obtain_data
from obtain_monomer_coordinates import clean_data

# run the pipeline
obtain_data()
clean_data(coordinate_database)