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

import acid_atom_to_num
import forge_database
from forge_database import get_file_paths_from_args
from forge_database import write_numpy_array

aa_dict = {"ALA": 1, "ARG": 2, "ASN": 3, "ASP": 4, "CYS": 5, "GLU": 6, "GLN": 7, "GLY": 8, "HIS": 9, "LIE": 10, "LEU": 11, "LYS": 12, "MET": 13, "PHE": 14, "PRO": 15, "SER": 16, "THR": 17, "TRP": 18, "TYR": 19, "VAL": 20}

atom_dict = {"C": 1, "N": 2, "O": 3, "ZN": 4}

""" first phase of the pipeline:
(1) extracts coordinaters from crystallographic information database
(2) feeds data into a csv file and create a new column containing
(3) stores the data in this csv file in a numpy array """
def coordinates():

    with open('protein_coordinates.csv', 'w') as csvfile:
        writer = csv.writer(csvfile)
        for path in get_file_paths_from_args():
            # read the crystallographic information file (uncompressing it on the fly)
            gemmi.read_structure(path, format=gemmi.CoorFormat.Detect)
            cif_file = cif.read(path)
            cif_block = cif_file.sole_block()

            """ obtain the following from each cif file:
            (1) the atom symbol: '_atom_site.type_symbol'
            (2) the monomer that the atom is a member of: '_atom_site.label_comp_id'
            (3) the x coordinate of the atom: 'atom_site.Cartn_x'
            (4) the y coordinate of the atom: 'atom_site.Cartn_y'
            (5) the z coordinate of the atom: 'atom_site.Cartn_z' """
            table = cif_block.find(['_atom_site.type_symbol', '_atom_site.label_comp_id', '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z'])

            # write (1), (2), (3), (4), (5) into a csv file with a new column containing the name of the cif file
            for row in table:
                writer.writerow(row)

    my_array = genfromtxt('protein_coordinates.csv', delimiter=',')

def names():
    # count the number of rows in the protein_coordinates csv file
    with open('protein_coordinates.csv', 'r') as coordinates_csv:
        for path in get_file_paths_from_args():
            reader = csv.reader(coordinates_csv)
            row_count = sum(1 for row in reader)

    with open('protein_names.csv', 'w') as names_csv:
        writer = csv.writer(names_csv)
        for path in get_file_paths_from_args():
        # read the crystallographic information file (uncompressing it on the fly)
            gemmi.read_structure(path, format=gemmi.CoorFormat.Detect)
            cif_file = cif.read(path)
            cif_block = cif_file.sole_block()

        """ obtain the following from each file:
        (1) the entry id """
        for i in range(row_count):
            writer.writerow(cif_block.name)

    my_array = genfromtxt('protein_names.csv', delimiter=',')

""" this is the second phase of the pipeline:
(1) mine the meta-data for atomic coordinates
(2) obtain the positions for every monomer within each protein
(3) add three new columns for the x, y, z coordinates of each monomer for every protein
(4) store the x, y, z positions for each monomer within the x, y, z columns respectively """
def AddColWithAAPositions(NParray):
    for col in range(3):
        AAPosition = 0
        for row in NParray:
            m = 0
            while NParray[row][5] == NParray[row + 1][5]:
                m = m + 1
                AAPosition = NParray[row][col] + NParray[row + 1][col]
        for row in NParray:
            while NParray[row][5] == NParray[row + 1][5]:
                NParray[row][col + 5] == AAPosition/m
