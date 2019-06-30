#lambdacheck.py
import numpy as np

S = np.array([[1,-1,0,0,0,0],[0,1,-1,-1,0,0],[0,0,2,0,-1,0],[0,0,0,1,1,-1]])

x = [1,1,1,1,1,7]

print('Test of row-wise matrix multiplication (vector multiplication).')

print('Using for loop:')
for i,row in enumerate(S):
    dx = sum(np.multiply(row,x))
    print('row:',i, 'row * x:',dx)

#lambdas = [(lambda x : sum(np.multiply(row,x))) for row in S]

lambdas = []
for row in S:
    lambdas.append(lambda x, row = row: sum(np.multiply(row,x)))

#lambdas = [(lambda row: lambda x: sum(np.multiply(row, x)))(row)]

print('Using list of lambda functions:')
for i,fun in enumerate(lambdas):
    dx = fun(x)
    print('row:',i,'f(x):',dx)
