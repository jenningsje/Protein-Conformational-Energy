# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# import items for the Mij matrix below
# import the acid table
schain_prob_table = open("schain_prob_table.txt")
schain_table_rows = sc_table.readlines()
schain_table.close()

# import items for the Hij matrix Below
# import the hamiltonian chart
hbond_type_table = open("hbond_type_table.txt")
hbond_type_table_rows = hbond_type_table.readlines()
hbond_type_table.close()

# import the amino acid list
aa_list = open("list.txt")
aa_string = aa_list.read()
aa_list.close()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# print items for the Mij matrix below
# print out the list of amino acids
print(" ")
print("amino acids")
amino_acids = schain_table_rows[0].split()
print(amino_acids)

# print items for the Hij matrix Below
# print out the list of amino acids
print(" ")
print("acids")
hbond_type_param = amino_acid_list.split()
print(seq_list)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# parameters for the Hij matrix below
n = len(hbond_type_param)
hbond_type_index = []

# parameters for the Mij matrix below
m = len(hbond_prob_param)
hbond_prob_index = []

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# matrices corresponding to the Mij matrix below
# hydrogen bond probability matrix
aa_table = []

# Mij matrix below
Mij = []

# matrices corresponding to the Hij matrix below
# matrix for the types of hydrogen bonds
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# construct the hydrogen bond probability matrix
for sc_table_row in sc_table_rows[1:]:
    
    sc_row = sc_table_row.split()[1:]
    numbers = list(map(float,sc_row))
    aa_table.append(numbers)

for item in hbond_prob_param:
    hbond_type_index.append(amino_acids.index(item))

# construct the Mij matrtix
for i in range(m):
    Mij.append([])
    for j in range(m):
        prob = aa_table[hb_prob_index[i]][hb_prob_index[j]]
        Mij[i].append(prob)

# construct the matrix for the types of hydrogen bonds
for line in lines[1:]:
    
    row = line.split()[1:]
    letters = list(map(str,row))
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

# tensor product of acid_table0 and Eij
for i in range(n):
    Aij.append([])
    for j in range(n):
            prob = aa_table[i][j]*Eij[i][j]
            Aij[i].append(prob)
