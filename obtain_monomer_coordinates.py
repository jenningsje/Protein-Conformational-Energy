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
import scanpy
import string
import numpy.ma as ma
import statistics
import pandas as pd
import math
import statistics
import numpy.ma as ma
import array as arr

import acid_atom_to_num
import forge_database
from forge_database import get_file_paths_from_args

# import text files

# import items for the Mij matrix Below
# import the acid table
table_file0 = open("schain_prob_table.txt")
lines0 = table_file0.readlines()
table_file0.close()

# import items for the Hij matrix Below
# import the hamiltonian chart
table_file = open("hbond_type_table.txt")
lines = table_file.readlines()
table_file.close()

# import the amino acid list
seq_file = open("list.txt")
seq_string = seq_file.read()
seq_file.close()

# create dictionaries

aa_dict_b = {b'ALA': 1, b'ARG': 2, b'ASN': 3, b'ASP': 4, b'CYS': 5, b'GLU': 6, b'GLN': 7, b'GLY': 8, b'HIS': 9, b'LIE': 10, b'LEU': 11, b'LYS': 12, b'MET': 13, b'PHE': 14, b'PRO': 15, b'SER': 16, b'THR': 17, b'TRP': 18, b'TYR': 19, b'VAL': 20}

atom_dict_b = {b'C': 1, b'N': 2, b'O': 3, b'ZN': 4, b'None': 5}

aa_filter = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "LIE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}

atom_filter = {'C', 'N', 'O', 'ZN'}

A = []

E = []

aa_table = []

h_table = []

# split the sidechain probability table
acids0 = lines0[0].split()

# split the amino acid list
seq_list = seq_string.split()

# parameters for the Hij matrix below
n = len(seq_list)
seq_index = []

# construct the hydrogen bond probability matrix
for line0 in lines0[1:]:
    row0 = line0.split()[1:]
    numbers0 = list(map(float,row0))
    aa_table.append(numbers0)

# construct the matrix for the types of hydrogen bonds
for line in lines[1:]:
    row = line.split()[1:]
    letters = list(map(str,row))
    h_table.append(letters)

# construct the matrix for the energy levels
for i in range(n):
    E.append([])
    for j in range(n):
        prob = h_table[i][j]
        if prob == "N":
            prob2 = -1.64013e-22
            E[i].append(prob2)
        elif prob == "O":
            prob2 = -2.09727e-22
            E[i].append(prob2)
        elif prob == "P":
            prob2 = 0.0
            E[i].append(prob2)
        else:
            prob2 = 0.0
            E[i].append(prob2)


# tensor product of acid_table0 and Eij
for i in range(n):
    A.append([])
    for j in range(n):
            prob = aa_table[i][j]*E[i][j]
            A[i].append(prob)

# tensor product of acid_table0 and Eij
for i in range(n):
    A.append([])
    for j in range(n):
            prob = aa_table[i][j]*E[i][j]
            A[i].append(prob)

# bioinformatics pipeline
def pipeline():

    """first stage of the pipeline: obtain the data"""

    # open a new csv file protein_coordinates.csv
    with open('protein_coordinates.csv','w') as csvfile:
        writer = csv.writer(csvfile)

        # loop throughout an entire directory? (check me)
        for path in get_file_paths_from_args():
            # read the crystallographic information file (uncompressing it on the fly)
            gemmi.read_structure(path, format=gemmi.CoorFormat.Detect)
            cif_file = cif.read(path)
            cif_block = cif_file.sole_block()

            # first phase of the bioinformatics pipeline:

            # obtain the following x, y, z coordinates, aa names and atoms for the protein given in path
            table = cif_block.find(['_atom_site.type_symbol', '_atom_site.label_comp_id',
            '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z', '_atom_site.pdbx_PDB_model_num'])

            # write each row in every cif file to the csv file protein_coordinates.csv
            for row in table:
                writer.writerow(row)

            arr = genfromtxt('protein_coordinates.csv', delimiter=',', dtype=object)

            """second stage of the pipeline: clean the data"""

            seq_indices = []

            for i in range(len(arr)):
                if arr[i][1] not in aa_dict_b:
                    seq_indices.append(b'None')
                elif arr[i][1] in aa_dict_b:
                    seq_indices.append(aa_dict_b[arr[i][1]])

            """third stage of the pipeline: analyze the data"""

            print(seq_indices)