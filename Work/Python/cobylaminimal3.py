#minimalcobyla3.py
import numpy as np
import scipy.optimize as opt

from scipy import sparse

from cobra.io.sbml import create_cobra_model_from_sbml_file
cobramodel = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')


def getObjectiveVector(cobramodel):
    C = [reaction.objective_coefficient for reaction in cobramodel.reactions]
    return C

C2 = getObjectiveVector(cobramodel)


C = [0,0,2,0,0,1]

#x0 = [0,0,0,0,0,0]
x0 = [0 for element in C]

lb = [-1,-1,-1,-1,-1,-1]

ub = [2,5,2,2,2,2]

S = np.array([[1,-1,0,0,0,0],[0,1,-1,-1,0,0],[0,0,2,0,-1,0],[0,0,0,1,1,-1]])

#S = sparse.csr_matrix(S)


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

ub_funcs = constrainFunctionsUB(ub)
lb_funcs = constrainFunctionsLB(lb)

ss_funcs_a = [(lambda x : sum(np.multiply(row,x))) for row in S]
ss_funcs_b = [(lambda x : -sum(np.multiply(row,x))) for row in S]


allconstr = ub_funcs + lb_funcs + ss_funcs_a + ss_funcs_b

#FBAres = opt.fmin_cobyla(FBAobjective,x0,allconstr,args = (C,),consargs =())
FBAres = opt.fmin_cobyla(FBAobjective,x0,allconstr,args = (C,),consargs =())
