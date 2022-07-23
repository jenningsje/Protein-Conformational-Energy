import numpy as np
import math
import statistics

# create the classes for the entropic corrections to the side chain and backbone entropies

def AcidPos(i,NParray):
    # filter for only a given amino acid
    mask = ((NParray[i][1]) == (NParray[i + 1][1]))
    masked_arr = NParray[mask,:]
    # find the amino acid coordinates
    return masked_arr

def delta(ri,rj,dmin,dmax):
    mask = (dmin <= ri[i] - rj[j] <= dmax)
    return mask

def sidechain(char_ij,prob):
    out = size(char)
    # Oxygen bonds
    mask1(char_ij == 'O')
    out[i][j] = 2.09727e-22
    # Nitrogen bonds
    mask2(char_ij == 'N')
    out[i][j] = -1.64013e-22
    # No hydrogen bonds present
    mask3(char_ij != 'O' and char_ij != 'N')

def backbone(NParray,r1,r2,dmin,dmax,i,j):
    out = size(NParray)
    # Beta conformation
    mask1 = (dmin <= r1[i] - r2[j] <= dmax
    and dmin <= r1[i + 1] - r2[j + 1] <= dmax
    and dmin <= r1[i + 2] - r2[j + 1] <= dmax)
    out[i][j] = 1
    # Beta conformation
    mask2 = (dmin <= r1[i] - r2[j] <= dmax
    and dmin <= r1[i + 1] - r2[j - 1] <= dmax
    and dmin <= r1[i + 2] - r2[j - 2] <= dmax)
    out[i][j]