# this function returns the path from a directory specified by the user otherwise it will return the pdb code
def get_file_paths_from_args():
    for arg in sys.argv[1:]:
        if os.path.isdir(arg):
            for path in CifWalk(arg):
                yield path
        else:
            yield expand_if_pdb_code(arg)

write_numpy_array():

    with open('protein_coordinate_database.csv', 'w') as csvfile:
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

            # write the items list in (1), (2), (3), (4), (5) into a database (in this case a csv file)
            # there is a new column added to this database that contains the name of the protein
            for row in table:
                writer.writerow(str(row) + str(os.path.basename(path)))
                print(n)
                n = n +1
                
    my_array = genfromtxt('protein_coordinate_database.csv', delimiter=',')