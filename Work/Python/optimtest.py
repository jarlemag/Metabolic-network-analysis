# -*- coding: cp1252 -*-

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

cobramodel.optimize(solver="gurobi")

def getObjectiveVector(cobramodel):
    C = [int(reaction.objective_coefficient) for reaction in cobramodel.reactions]
    return C

#C = getObjectiveVector(cobramodel) #Load the objective vector from the cobra model

#C = [0 for x in xrange(len(cobramodel.reactions))] #Didn't help


C = [0 for x in xrange(92)]
C[84] = 1

'''
def FBAobjective(x,C, sense = -1):
    return  np.dot(np.array(C),np.array(x))*sense
'''
def FBAobjective(x,C, sense = -1):

    #print 'This is a warning!'
    #print 'The enemy is watching!'
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
#x0 = [0 for reaction in cobramodel.reactions]

x0 = [0 for element in C]

#S = cobramodel.S
S = cobramodel.S.toarray()

#ss_funcs_a = [(lambda x : sum(np.multiply(row,x))) for row in S]
#ss_funcs_b = [(lambda x : -sum(np.multiply(row,x))) for row in S]


ss_funcs_a = []
for row in S:
    ss_funcs_a.append((lambda x, row = row: sum(np.multiply(row,x))))

ss_funcs_b = []
for row in S:
    ss_funcs_b.append((lambda x, row = row: -sum(np.multiply(row,x))))


print 'About to start.'
#allconstr = ub_funcs + lb_funcs + ss_funcs_a + ss_funcs_b
allconstr = ub_funcs + lb_funcs + ss_funcs_a + ss_funcs_b
FBAres = opt.fmin_cobyla(FBAobjective,x0,allconstr,args = (C,),consargs =())

reactions = cobramodel.reactions
rids = [reaction.id for reaction in reactions]
v = zip(rids,FBAres)

for a,b in v:
	print a, "%.1f" %b


metabolites = cobramodel.metabolites
#Check steady-state criteria:
num = 0
violators = []
violatorcount = 0
for row in S:
    num +=1
    dxdt = 0
    for index,coef in enumerate(row):
        dxdt += coef*FBAres[index]

    print 'metabolite', num, "%.5f" %dxdt
    if abs(dxdt) > 0.001:
        violators.append(num-1)
        violatorcount +=1
        print 'metabolite',num,'(',metabolites[num-1].id,')','added to violator list.'

print 'Number of steady state violations:',violatorcount


#Verify that mass conservation is violated.

print 'Metabolites violating mass conservation:'
violator_ids = []

for met in violators:
    met_id = metabolites[met].id
    print met_id
    violator_ids.append(metabolites[met].id)
    violator_metabolite = metabolites[met]
    violator_reactions = violator_metabolite._reaction
    metflux = 0
    for reaction in violator_reactions:
        reaction_id = reaction.id
        metflux += reaction.get_coefficient(met_id)*FBAres[reactions.index(reaction_id)]
    print 'metabolite:',met_id,'total flux:',metflux


print 'Number of metabolites violating steady state:', len(violators)

reactions = cobramodel.reactions
reaction_ids = [reaction.id for reaction in reactions]

#Check for violation of lower/upper bounds rule:

beyondbounds = 0
for index,flux in enumerate(FBAres):
    if (flux > ub[index]+.001) or (flux < lb[index]-0.001):
        print 'reaction:',reaction_ids[index]
        print 'lower bound:',lb[index]
        print 'upper bound:',ub[index]
        print 'flux:',flux
        beyondbounds +=1
print '# of reactions beyond bounds:',beyondbounds



#Make equality constraints:

ss_funcs = [(lambda x : sum(np.multiply(row,x))) for row in S]


#res = opt.minimize(FBAobjective,x0,args = (C,),method = ‘SLSQP’, constraints = minSLQPcons)
