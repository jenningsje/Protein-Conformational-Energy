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

def write_numpy_array(data):

    with open(data, 'w') as csvfile:
        writer = csv.writer(csvfile)
        n = 0
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
                writer.writerow(str(row) + str(os.path.basename(path)))
                print(n)
                n = n +1

    numpy_array = genfromtxt('protein_coordinate_database.csv', delimiter=',', usecols=np.arange(3))

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


def obtain_monomer_coordinates():
    write_numpy_array()
    AddColWithAAPositions(numpy_array)
