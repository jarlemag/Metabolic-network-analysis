#optimtest.py
#Testing optimization functions


import scipy.optimize as opt

import numpy as np

from cobra.io.sbml import create_cobra_model_from_sbml_file
cobramodel = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')
import loadData as load
fluxvalues = load.ExpFluxesfromXML('expdata.xml','Perrenoud','Batch','aerobe')
rmap = load.ReactionMapfromXML('reactionmaps.xml','Perrenoud','SCHUETZR')


def objective(x):
    return x[0]*x[1]

def constr1(x):
    return 1 - x[0]**2 - x[1]**2

def constr2(x):
    return x[1]

def constr3a(x):
    return x[0] - x[1]

def constr3b(x):
    return x[1] - x[0]

q = opt.fmin_cobyla(objective, [0.0,0.1],[constr1,constr2], disp = 0)

z = opt.fmin_cobyla(objective, [0.0,0.1],[constr1,constr2,constr3a,constr3b], disp = 0)



#Solving an FBA problem:

cobramodel.optimize()

def getObjectiveVector(cobramodel):
    C = [reaction.objective_coefficient for reaction in cobramodel.reactions]
    return C

C = getObjectiveVector(cobramodel)

'''
def FBAobjective(x,C, sense = -1):
    return  np.dot(np.array(C),np.array(x))*sense
'''
def FBAobjective(x,C, sense = -1):


    print 'This is a warning!'
    print 'The enemy is watching!'
    return  np.dot(C,x)*sense

def getUpperBounds(model):
    return [reaction.upper_bound for reaction in model.reactions]

def getLowerBounds(model):
    return [reaction.lower_bound for reaction in model.reactions]

ub = getUpperBounds(cobramodel)

lb = getLowerBounds(cobramodel)
 
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


def constrainFunctions(ub,lb):
    ub_funcs = []
    lb_funcs = []
    for index, upperbound in enumerate(ub):
        def upperlimit(x, index = index, upperbound = upperbound):
            return upperbound - x[index]
        ub_funcs.append(upperlimit)
    for index, lowerbound in enumerate(lb):
        def lowerlimit(x, index = index, lowerbound = lowerbound):
            return x[index] - lowerbound
        lb_funcs.append(lowerlimit)
    return ub_funcs, lb_funcs

ub_funcs = constrainFunctionsUB(ub)
lb_funcs = constrainFunctionsLB(lb)

cobramodel.to_array_based_model()

#x0 = np.array([0 for reaction in cobramodel.reactions])
x0 = [0 for reaction in cobramodel.reactions]

S = cobramodel.S
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

ss_funcs_a = ssConsFuncspos(S)
ss_funcs_b = ssConsFuncsneg(S)

print 'About to start.'
allconstr = ub_funcs + lb_funcs + ss_funcs_a + ss_funcs_b
FBAres = opt.fmin_cobyla(FBAobjective,x0,allconstr,args = (C,),consargs =())
