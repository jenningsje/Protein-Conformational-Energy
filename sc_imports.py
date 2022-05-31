"""import items for the Mij matrix below import the acid table this is a 20x20 matrix that determines 
    the probability that a pair of any of the 20 amino acids will form a hydrogen bond"""
sc_prb_table = open("schain_prob_table.txt")
sc_prb_rows = sc_prb_table.readlines()
sc_prb_table.close()

"""import items for the Hij matrix Below importthe hamiltonian chart this is 20x20 matrix determines
    what kind of hydrogen bond (O--H bond, N--H bond or neither) will form between the sidechains 
    of two amino acids"""
hbond_type_table = open("hbond_type_table.txt")
hb_type_rows = hbond_type_table.readlines()
hbond_type_table.close()

# import the amino acid list
amino_acids_list = open("list.txt")
acid_list = amino_acids_list.read()
amino_acids_list.close()

"""print items for the Mij matrix below
print out the probability matrix"""
print(" ")
print("probs")
list_of_amino_acids = sc_prb_rows[0].split()
print(list_of_amino_acids)

"""print items for the Hij matrix Below
print out the list of amino acids"""
print(" ")
print("hbonds")
aa_list = acid_list.split()
print(aa_list)

# parameters for the Hij matrix below
n = len(aa_list)
hbond_type_index = []

"""matrices corresponding to the Mij matrix below
hydrogen bond probability matrix"""
aa_table = []

# Mij matrix below
Mij = []

"""matrices corresponding to the Hij matrix below
matrix for the types of hydrogen bonds"""
hamiltonian_table = []

# matrix for the energy levels
Eij = []

# hamiltonian matrix
Hij = []

# tensor product for acid_table0 and Eij
Aij = []

# matrix for the Kronecker delta
Rij = []

# tensor product for Mij and Kij
Dij = []

# Energy corresponding to the sidechains
Hsc = []

# construct the hydrogen bond probability matrix
for sc_prb_row in sc_prb_rows[1:]:

    split_sc_prb_row = sc_prb_row.split()[1:]
    numbers = list(map(float,split_sc_prb_row))
    aa_table.append(numbers)

# construct the matrix for the types of hydrogen bonds
for hb_type_row in hb_type_rows[1:]:

    split_hb_type_row = hb_type_row.split()[1:]
    letters = list(map(str,split_hb_type_row))
    hamiltonian_table.append(letters)

# construct the matrix for the energy levels
for i in range(n):
    Eij.append([])
    for j in range(n):
        prob = hamiltonian_table[i][j]
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