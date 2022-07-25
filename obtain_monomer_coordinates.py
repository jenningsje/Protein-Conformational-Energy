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

import acid_atom_to_num
import forge_database
from forge_database import get_file_paths_from_args
from entropies import AcidPos

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

aa_dict = {'ALA': 1, 'ARG': 2, 'ASN': 3, 'ASP': 4, 'CYS': 5, 'GLU': 6, 'GLN': 7, 'GLY': 8, 'HIS': 9, 'LIE': 10, 'LEU': 11, 'LYS': 12, 'MET': 13, 'PHE': 14, 'PRO': 15, 'SER': 16, 'THR': 17, 'TRP': 18, 'TYR': 19, 'VAL': 20}

atom_dict_b = {'C': 1, 'N': 2, 'O': 3, 'ZN': 4}

A = []

E = []

aa_table = []

h_table = []

M = []

seq_indices = []

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

arr1 = [[b'_atom_site.type_symbol', b'_atom_site.label_seq_id', b'_atom_site.label_comp_id', b'_atom_site.Cartn_x', b'_atom_site.Cartn_y', b'_atom_site.Cartn_z', b'_atom_site.pdbx_PDB_model_num']]

"""bioinformatics pipeline"""

# first stage of the pipeline obtain the data
def obtain_data():

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
            table = cif_block.find(['_atom_site.type_symbol', '_atom_site.label_seq_id', '_atom_site.label_comp_id',
            '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z', '_atom_site.pdbx_PDB_model_num'])

            # write each row in every cif file to the csv file protein_coordinates.csv
            for row in table:
                writer.writerow(row)

            arr2 = genfromtxt('protein_coordinates.csv', delimiter=',', dtype=object)

            np.append(arr1, arr2, axis=0)

            DF = pd.DataFrame(arr2)

            DF.to_csv('coordinate_database')

def clean_data():

    # amino acid indices
    aa_dict_b = {b'ALA': 1, b'ARG': 2, b'ASN': 3, b'ASP': 4, 
    b'CYS': 5, b'GLU': 6, b'GLN': 7, b'GLY': 8, b'HIS': 9, 
    b'LIE': 10, b'LEU': 11, b'LYS': 12, b'MET': 13, 
    b'PHE': 14, b'PRO': 15, b'SER': 16, b'THR': 17, b'TRP': 
    18, b'TYR': 19, b'VAL': 20}

    # atom indices
    atom_dict_b = {b'C': 1, b'N': 2, b'O': 3, b'ZN': 4}

    ma.masked_where(np.isin(database, aa_filter), arr)
    ma.masked_Where(np.isin(database, atom_filter, arr))

    return database