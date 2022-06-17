
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
from forge_database import write_numpy_array

aa_dict = {"b'ALA'": 1, "b'ARG'": 2, "b'ASN'": 3, "b'ASP'": 4, "b'CYS'": 5, "b'GLU'": 6, "b'GLN'": 7, "b'GLY'": 8, "b'HIS'": 9, "b'LIE'": 10, "b'LEU'": 11, "b'LYS'": 12, "b'MET'": 13, "b'PHE'": 14, "bPRO": 15, "bSER": 16, "bTHR": 17, "bTRP": 18, "bTYR": 19, "b'VAL'": 20}

atom_dict = {"C": 1, "N": 2, "O": 3, "ZN": 4}

# returns the average value of matrix elements

def coordinates(M, i , j):
    m, x_sum, y_sum, z_sum = 0
    while float(M[1][j]) == float(M[1][j + 1]):
        m = m + 1
        x_sum = M[i][3] + M[i + 1][3]
        y_sum = M[i][4] + M[i + 1][4]
        z_sum = M[i][5] + M[i + 1][5]
        return m, x_sum, y_sum, z_sum


# bioinformatics pipeline
def pipeline():

    # open a new csv file protein_coordinates.csv
    with open('protein_coordinates.csv','w') as csvfile:
        writer = csv.writer(csvfile)

        # loop throughout an entire directory? (check me)
        for path in get_file_paths_from_args():
            # read the crystallographic information file (uncompressing it on the fly)
            gemmi.read_structure(path, format=gemmi.CoorFormat.Detect)
            cif_file = cif.read(path)
            cif_block = cif_file.sole_block()

            # obtain the following x, y, z coordinates, aa names and atoms for the protein given in path
            table = cif_block.find(['_atom_site.type_symbol', '_atom_site.label_comp_id',
            '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z'])

            # write each row in every cif file to the csv file protein_coordinates.csv
            for row in table:
                writer.writerow(row)

            # construct the block for the coordinate data below
            arr = genfromtxt('protein_coordinates.csv', delimiter=',', dtype=object)

            print (arr.shape)
            print (arr)
            print (arr[0,0])

            print (type(arr[0,0]))
                    
            test = arr[:,0]
            print (test)

            mask = (arr[:,1] == b'VAL')
            test = arr[mask,:]
            print (test.shape)
            print (test)

            test2 = test[:,2]
            print (test2)

            test3 = test2.astype(float)
            print (test3)

            print (test3.mean())
