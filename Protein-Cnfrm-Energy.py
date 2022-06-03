#!/usr/bin/env python

""" Obtain the conformational entropy for every protein (cif file or gzipped cif file) within a given directory """
from __future__ import print_function
import sys
import os
from collections import Counter
from gemmi import cif, CifWalk, expand_if_pdb_code
import numpy as np
from numpy import genfromtxt

""" import the following files containing:
(1) the list of amino acids
(2) the table for the different types of combinations of hydrogen bonds between a pair of amino acids (20x20 matrix)
(3) the table for the probability that a given hydrogen bond will form between two amino acids (20x20 matrix)
(4) import the variables for the sidechain sidechain matrix 
(5) import the functions necessary for constructing the """
import sc_imports
from obtain_monomer_coordinates import write_numpy_array

# get the names of each protein as well as the monomer coordinates for them
write_numpy_array()