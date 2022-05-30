#!/usr/bin/env python
# Check amino-acid frequency in the PDB database (or it's subset)
# by reading meta-data from mmCIF files.

from __future__ import print_function
import sys
import os
from collections import Counter
from gemmi import cif, CifWalk, expand_if_pdb_code

# import the following files (sidechain_files) for:
# the list of amino acids
# the table for the different types of combinations of hydrogen bonds for a pair of amino acids (20x20 matrix)
# the table for the probability that a given hydrogen bond will form between two amino acids (20x20 matrix)
# import the variables for the sidechain sidechain matrix
import sidechain_files
import sidechain_variables


# this function returns the path from a directory
# otherwise it will return the pdb code
def get_file_paths_from_args():
    for arg in sys.argv[1:]:
        if os.path.isdir(arg):
            for path in CifWalk(arg):
                yield path
        else:
            yield expand_if_pdb_code(arg)

n = 0
for path in get_file_paths_from_args():
    # read file (uncompressing on the fly) and get the only block
    print(n)
    n = n + 1
    block = cif.read(path).sole_block()
    # find table with the sequence
    seq = block.find('_entity_poly_seq.', ['entity_id', 'mon_id'])
