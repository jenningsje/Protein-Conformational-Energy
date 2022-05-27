import numpy as np

database = open("pdb_seqres.txt", 'r')
while True:
    line = database.readline()
    if (!(line) == '>' || line.len())

        #length of an amino acid sequence
        aa_seq_length = line.len()

        # kronecker delta matrix
        for i in range(aa_seq_length):
                Rij.append([])
                for j in range(aa_seq_length):
                    ith_index = i
                    jth_index = j
                    Rij[i].append("d(r" + str(ith_index) + " " +"- r" + str(jth_index) + " " + "- r)" )

        # hamiltonian matrix
        for i in range(aa_seq_length):
            Hij.append([])
            for j in range(aa_seq_length):
                prob = Eij[seq_index0[i]][seq_index0[j]]
                Hij[i].append(prob)

        # tensor product of Mij and Kij
        for i in range(aa_sequence_length):
            Dij.append([])
            for j in range(m):
                    prob = Mij[i][j]
                    Dij[i].append(prob)

        # steric hinderance correction
        for i in range(aa_seq_length):
            Es.append([])
            for j in range(aa_seq_length):
                    prob = str(2.26006e-22) + str(Mij[i][j]) + Rij[i][j]
                    Es[i].append(prob)

        # Energy corresponding to the sidechains
        for i in range(aa_seq_length):
            Hsc.append([])
            for j in range(aa_seq_length):
                prob = Hij[i][j] * Dij[i][j]
                if i == j + 1 or j == i + 1:
                    Hsc[i].append(prob)
                else:
                    Hsc[i].append(0)

database.close()
