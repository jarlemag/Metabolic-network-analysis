import numpy as np
#import scipy as sp
from numpy.linalg import svd
import scitools
import sympy
#Calculability analysis and calculation of unknown rates from known rates, based on Klamt et al. 2002.

#Consider making pull request for bug fix in  scitools.numpyutils.matrix_rank() (Same function now exists in NumPy)


def nullspace(A):
    #Finds a basis for the nullspace, with each column vector normalized. To find a rational basis, use SymPy.
    u,s,vh = np.linalg.svd(A)
    r = len(s)
    n = np.shape(A)[1]
    return vh[-(n-r):] 


def isdetermined(A, verbose = False):
    u= np.shape(A)[1]
    rank =  np.linalg.matrix_rank(A) 
    if rank < u:
        if verbose:
            print 'The system is undetermined.'
        return False
    elif rank == u:
        if verbose:
            print 'The system is determined.'
        return True
    else:
        raise Exception('Rank larger than number of reactions.')


def isredundant(A, verbose = False):
    m = np.shape(A)[0]
    rank =  np.linalg.matrix_rank(A)
    if rank < m:
        if verbose:
            print 'The system is redundant.'
        return True
    elif rank == m:
        if verbose:
            print 'The system is not redundant.'
        return False
    else:
        raise Exception('Rank larger than number of metabolites.')



def calculaterates(N_known, N_unknown, r_known):
    if isdetermined(N_unknown) and not isredundant(N_unknown):
        r_unknown = np.linalg.inv(N_unknown).dot(N_known).dot(r_known)
    else:
        Npi = np.linalg.pinv(N_unknown) #Calculate the Penrose pseduo-inverseof N_unknown
        r_unknown = - Npi.dot(N_known).dot(r_known)
    return r_unknown


def redundancyMatrix(N_known,N_unknown,tol = 1e-6):
    R = N_known - N_unknown.dot(np.linalg.pinv(N_unknown)).dot(N_known)
    R[abs(R) < tol] = 0
    return R

def findBalanceableRates(R):
    pass
    #Find balancable rates by inspecting the redundancy matrix R

def calculabilityAnalysis(N_known,N_unknown,r_known = None,makeconsistent = True, verbose = False):
    calculable = None
    r_unknown = None
    #0. Check if the system is fully determined.
    if isdetermined(N_unknown) and not isredundant(N_unknown):
        if verbose:
            print 'The system is fully determined.'
        if r_known is not None:
        r_unknown = calculaterates(N_known,N_unknown,r_known)
        return None,r_unknown
    #1.Balancing:
    #a. Determine if the system is redundant.
    if not isredundant(N_unknown, verbose = verbose): #If the system is not redundant, go to step 2.
        pass
    else:
        pass
        #Determine the balancable rates by analyzing the redundancy matrix R.
        #Thereafter,balance the rats to obtain a consistent system.

    #2. Calculating rates:
    #a. Determine which rates of r_known that are calculable by searching for null rows in the 
    #null space matrix K_n
    K_n = np.array(sympy.Matrix(N_unknown).nullspace()).transpose()
    calculable = [x for x in range(len(K_n)) if x not in np.nonzero(K_n)[0]]
    print 'calculable fluxes:',calculable #Need to ensure the index are reported correctly, taking into account 0-indexing and the known reactions
        
    #b. Take the values of the calculable rates from the least-squares solution.
    if r_known is not None:
        r_unknown = calculaterates(N_known,N_unknown,r_known)

    return calculable, r_unknown
    #To do: Output the equations for the calculable fluxes.
    
if __name__ == '__main__':
    
    #Examples of determining calculability and balanceability, following Klamt et al.:

    # Define the m X n (5 x 7) matrix N:
    N = np.array([[1,-1,0,-1,0,0,0],[0,1,-1,0,0,0,0],
                  [0,1,0,1,-1,0,0],[0,0,1,1,0,-1,0],[0,0,1,0,0,0,-1]])


    #Check the nullspace presented on Page 3:

    K = np.array([[4,1,1,3,4,4,1],[0,1,1,-1,0,0,1]]).transpose()

    N.dot(K)

    #Finding a rational nullspace basis using SymPy:
    
    basis = np.array(sympy.Matrix(N).nullspace())

    #Case 1- r_known is empty:


    #Case 2 - r_known = r_1:

    N_known = np.array([1,0,0,0,0])

    N_unknown = N[:,1:]


    a, b = calculabilityAnalysis(N_known,N_unknown) 


    #Case 3 - r_known = (r_1,r_7)^T :

    N_known = np.array([[1,0,0,0,0],[0,0,0,0-1]])
    N_unknown= N[:,1:-1]
    r_known = np.array([1,0.3])

    c, d = calculabilityAnalysis(N_known,N_unknown, verbose = True)


    #Case 4 - r_known = (r_1,r_5) :

    N_known = np.array([[1,0,0,0,0],[0,0,-1,0,0]])
    N_unknown = N[:,1:4] + N[:,-1:-2] #Need to find the correct slice here
    known_reactions = [1,5] #Indices of the known reactions. Reactions 1 and 5 are known.
    r_known = [2,3]   #Change this to be a dict of index and flux value?



    #Research question: Can calculability analysis be used to correct for differences in experimental model and computational model when comparng experimental results to model results?
