# import items for the Mij matrix below
# import the acid table
schain_prob_table = open("schain_prob_table.txt")
schain_table_rows = schain_prob_table.readlines()
schain_prob_table.close()

# import items for the Hij matrix Below
# import the hamiltonian chart
hbond_type_table = open("hbond_type_table.txt")
hbond_type_table_rows = hbond_type_table.readlines()
hbond_type_table.close()

# import the amino acid list
aa_list = open("list.txt")
aa_string = aa_list.read()
aa_list.close()
