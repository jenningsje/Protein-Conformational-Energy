import numpy as np
import math

# create the classes for the entropic corrections to the side chain and backbone entropies

class Entropy:

    NHbond = -6.66101*(10**(-19))
    OHbond = -5.21342*(10**(-19))

    def AcidPos(i,NParray):
    # filter for only a given amino acid
        mask = (NParray[i][2]) == float(NParray[i + 1][2])
        masked_arr = NParray[mask,:]
    # find the amino acid coordinates
        x_amino = mean(masked_arr[i][3])
        y_amino = mean(masked_arr[i][4])
        z_amino = mean(masked_arr[i][5])
    # find the distance of the amino acid from the origin
        R = ((x_amino, 2) + pow(y_amino, 2) + pow(z_amino, 2))
        return R

    def delta(ri,rj,dmin,dmax):
        mask = (dmin <= ri - rj <= dmax)