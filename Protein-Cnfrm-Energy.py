#!/usr/bin/env python

""" Obtain the conformational entropy for every protein (cif file or gzipped cif file) within a given directory """
from __future__ import print_function
import sys
import os
from collections import Counter
import gemmi
from gemmi import cif, CifWalk, expand_if_pdb_code
import numpy as np
import csv
from csv import DictWriter
from numpy import genfromtxt

""" import the following files containing:
(1) the list of amino acids
(2) the table for the different types of combinations of hydrogen bonds between a pair of amino acids (20x20 matrix)
(3) the table for the probability that a given hydrogen bond will form between two amino acids (20x20 matrix)
(4) import the variables for the sidechain sidechain matrix """
import sc_imports
import forge_database

write_numpy_array()