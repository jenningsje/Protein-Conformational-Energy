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

import acid_atom_to_num
import forge_database
from forge_database import get_file_paths_from_args
from forge_database import write_numpy_arr

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

aa_dict = {b'ALA': 1, b'ARG': 2, b'ASN': 3, b'ASP': 4, b'CYS': 5, b'GLU': 6, b'GLN': 7, b'GLY': 8, b'HIS': 9, b'LIE': 10, b'LEU': 11, b'LYS': 12, b'MET': 13, b'PHE': 14, b'PRO': 15, b'SER': 16, b'THR': 17, b'TRP': 18, b'TYR': 19, b'VAL': 20, b'None': 21}

atom_dict = {b'C': 1, b'N': 2, b'O': 3, b'ZN': 4, b'None': 5}

atom_filter = {b'C' , b'N' , b'O' , b'ZN'}

aa_filter = {b'ALA', b'ARG', b'ASN', b'ASP', b'CYS', b'GLU', b'GLN', b'GLY', b'HIS', b'LIE', b'LEU', b'LYS', b'MET', b'PHE', b'PRO', b'SER', b'THR', b'TRP', b'TYR', b'VAL'}

Aij = []

Eij = []

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
    Eij.append([])
    for j in range(n):
        prob = h_table[i][j]
        if prob == "N":
            prob2 = -1.64013e-22
            Eij[i].append(prob2)
        elif prob == "O":
            prob2 = -2.09727e-22
            Eij[i].append(prob2)
        elif prob == "P":
            prob2 = 0.0
            Eij[i].append(prob2)
        else:
            prob2 = 0.0
            Eij[i].append(prob2)


# tensor product of acid_table0 and Eij
for i in range(n):
    Aij.append([])
    for j in range(n):
            prob = aa_table[i][j]*Eij[i][j]
            Aij[i].append(prob)

# tensor product of acid_table0 and Eij
for i in range(n):
    Aij.append([])
    for j in range(n):
            prob = aa_table[i][j]*Eij[i][j]
            Aij[i].append(prob)

# bioinformatics pipeline
def pipeline():
    count = 0
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
            '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z'])

            # write each row in every cif file to the csv file protein_coordinates.csv
            for row in table:
                writer.writerow(row)

            # construct the block for the coordinate data below
            arr = genfromtxt('protein_coordinates.csv', delimiter=',', dtype=object)
            
            print(count)
            count = count + 1

            # second phase of the bioinformatics pipeline
            file_name = os.path.basename(path)

            # Mij matrix below
            Mij = []

            # hamiltonian matrix
            Hij = []

            # matrix for the Kronecker delta
            Rij = []

            # tensor product for Mij and Kij
            Dij = []

            # Energy corresponding to the sidechains
            seq_indices = []

            # indicate where the corrupt data is located
            for amino in aa_filter:
                arr = np.where(arr != amino, b'None', arr)
            for atom in atom_filter:
                arr = np.where(arr != atom, b'None', arr)

            # replace the amino acids and the atoms with indices
            for amino in arr[:,0]:
                seq_indices.append(arr[aa_dict[amino]])
            for atom in arr[:,1]:
                seq_indices.append(arr[atom_dict[atom]])

            # 