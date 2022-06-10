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
from obtain_monomer_coordinates import coordinates
from obtain_monomer_coordinates import names

# get the names of each protein as well as the monomer coordinates for them
coordinates()
names()
