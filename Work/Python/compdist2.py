#compdist2.py
#Rewrite of compdist


import numpy as np
import scipy.io

def compdist2(fluxvector,vector2 = None,options = None,sense = None):

    expdata = scipy.io.loadmat('expdata.mat')


#Fluxvalues = expdata.perrenoud.abs.batch.aerobe.fluxvalues
# # fluxvalues is 4 levels below perrenoud.
    perrenoud = expdata['expdata']['perrenoud']

    fluxvalarray = perrenoud[0][0][0][0][0][0][0][0][0][0][0][0][0][0] #14 levels!    
    #firstrow = perrenoud[0][0][0][0][0][0][0][0][0][0][0][0][0][0][0]
    fluxvalues = [row[0] for row in fluxvalarray]


    dist = np.linalg.norm(np.array(fluxvector)-np.array(fluxvalues))
    return dist


