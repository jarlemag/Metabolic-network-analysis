#minimalcobyla.py
import numpy as np
import scipy.optimize as opt

C = [0,0,2,0,0,1]

x0 = [0,0,0,0,0,0]

lb = [-1,-1,-1,-1,-1,-1]

ub = [2,5,2,2,2,2]

S = np.array([[1,-1,0,0,0,0],[0,1,-1,-1,0,0],[0,0,2,0,-1,0],[0,0,0,1,1,-1]])

def FBAobjective(x,C, sense = -1):
    return  np.dot(C,x)*sense

def constrainFunctionsUB(ub):
    upperboundfuncs = []
    for index,upperbound in enumerate(ub):
        def upperlimit(x, index = index, upperbound = upperbound):
            return upperbound - x[index]
        upperboundfuncs.append(upperlimit)
    return upperboundfuncs

def constrainFunctionsLB(lb):
    lowerboundfuncs = []
    for index,lowerbound in enumerate(lb):
        def lowerlimit(x, index = index, lowerbound = lowerbound):
            return x[index] - lowerbound
        lowerboundfuncs.append(lowerlimit)
    return lowerboundfuncs

def ssConsFuncspos(S):
    ss_funcs = []
    for row in S:
        def steadystatconstr(x, row = row ):
            return sum(np.multiply(row,x))
        ss_funcs.append(steadystatconstr)
    return ss_funcs

def ssConsFuncsneg(S):
    ss_funcs = []
    for row in S:
        def steadystatconstr(x, row = row ):
            return -sum(np.multiply(row,x))
        ss_funcs.append(steadystatconstr)
    return ss_funcs


ub_funcs = constrainFunctionsUB(ub)
lb_funcs = constrainFunctionsLB(lb)

ss_funcs_a = ssConsFuncspos(S)
ss_funcs_b = ssConsFuncsneg(S)

allconstr = ub_funcs + lb_funcs + ss_funcs_a + ss_funcs_b

FBAres = opt.fmin_cobyla(FBAobjective,x0,allconstr,args = (C,),consargs =())
