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
