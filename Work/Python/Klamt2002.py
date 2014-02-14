import numpy as np
#import scipy as sp
from numpy.linalg import svd
import scitools
import sympy


# N = mxn =  5x7
N = np.array([[1,-1,0,-1,0,0,0],[0,1,-1,0,0,0,0],
              [0,1,0,1,-1,0,0],[0,0,1,1,0,-1,0],[0,0,1,0,0,0,-1]])

u,s,vh = svd(N)

rank = len(s) #The rank of a matrix equals the number of singular values of A
r = rank

#The null space of N is spanned by the last n - r columns of V:

n = 7
# n - r = 7 - 5 = 2 #The null space of N is spanned by the two last columns of V.

null = vh.transpose()[-(n-r):]
null[null < 1e-6] = 0


def rank(A):
    u,s,v = np.linalg.svd(A)
    return len(s)


def nullspace(A):
    #Finds a basis for the nullspace, but it might not be "pretty"
    u,s,vh = np.linalg.svd(A)
    r = len(s)
    n = np.shape(A)[1]
    return vh[-(n-r):] 
    

#or#:


R = scitools.numpyutils.matrix_rank(N)
#Make pull request for bug fix

#or:
R = np.linalg.matrix_rank(N)


#Check the nullspace presented on Page 3:

K = np.array([[4,1,1,3,4,4,1],[0,1,1,-1,0,0,1]]).transpose()

N.dot(K)

#Multiply each column with a random number, then check again:

K2 = K.copy()
K2[:,0] *= np.random.randn()
K2[:,1] *= np.random.randn()

print N.dot(K2)


#Finding a rational nullspace basis using SymPy:

basis = sympy.Matrix(N).nullspace()
print N.dot(basis.transpose())
