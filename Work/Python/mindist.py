#mindist.py

import scipy.optimize as optimize
import numpy as np
from cobra.io.sbml import create_cobra_model_from_sbml_file
import compdist
cobramodel = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')
cobramodel.optimize(solver='gurobi')


def dictToLists(fluxdict):
    values = []
    keys = []
    for key in fluxdict:
        keys.append(key)
        values.append(fluxdict[key])
    return keys, values


def dictToArrays(fluxdict):
    values = np.array([])
    keys = np.array([])
    for key in fluxdict:
        keys.append(key)
        values.append(fluxdict[key])
    return keys, values


def FBAobjective(x,C, sense = -1):
    return  np.dot(C,x)*sense

def getObjectiveVector(cobramodel):
    C = [int(reaction.objective_coefficient) for reaction in cobramodel.reactions]
    return C

def getUpperBounds(model):
    return [reaction.upper_bound for reaction in model.reactions]

def getLowerBounds(model):
    return [reaction.lower_bound for reaction in model.reactions]


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


def objectiveValue(x,C):
    return  np.dot(np.array(C),np.array(x))


def objectiveConstraint(x,objective,value,optreq):
    return

def compdistcomplete(rawfluxvector,model, debug = False):
    if debug:
        print 'compdistcomplete: Calculating distance...'
    import loadData as load
    import extractflux
    #print 'Loading experimental data...'
    expdata = load.ExpFluxesfromXML('expdata.xml','Perrenoud','Batch','aerobe')
    #print 'Loading reaction map...'
    rmap = load.ReactionMapfromXML('reactionmaps.xml','Perrenoud','SCHUETZR')

    tup = zip([reaction.id for reaction in model.reactions],rawfluxvector,) #
    rawfluxdict = {reactionid:fluxvalue for reactionid,fluxvalue in tup}
    extractfluxdict = extractflux.extractfluxdict(rawfluxdict,rmap)
    #extractfluxvector = [extractfluxdict[reaction.id] for reaction in model.reactions]
    fluxvector = []
    fluxvalues = []

    for key in extractfluxdict:
        if (key in expdata):
            fluxvector.append(extractfluxdict[key])
            fluxvalues.append(float(expdata[key]))
 
    dist = np.linalg.norm(np.array(fluxvector)-np.array(fluxvalues))

    return dist


def mindist(model,reactionmap,expfluxdict,splitsmap = None, verbose = False, debug = False, optreq = 1, optsplits = False):
    pass


C = getObjectiveVector(cobramodel)

cobramodel.to_array_based_model()
S = cobramodel.S.toarray()


ub = getUpperBounds(cobramodel)
lb = getLowerBounds(cobramodel)

ub_funcs = constrainFunctionsUB(ub) 
lb_funcs = constrainFunctionsLB(lb)


x0 = [0 for element in C]

#x0 =  cobramodel.solution.x

optreq = 1

ss_funcs_a = []
tol = 0.001 #Tolerance in steady-state constraints (zero tolerance creates numerical errors)
for row in S:
    ss_funcs_a.append((lambda x, row = row: tol + sum(np.multiply(row,x))))

ss_funcs_b = []
for row in S:
    ss_funcs_b.append((lambda x, row = row: tol -sum(np.multiply(row,x))))


allconstr = ub_funcs + lb_funcs + ss_funcs_a + ss_funcs_b

print 'Optimizing FBA objective with COBYLA.'

FBAres = optimize.fmin_cobyla(FBAobjective,x0,allconstr,args = (C,),consargs =())

FBAobjval = objectiveValue(FBAres,C)
    

objcon = lambda x, objective = C , value = FBAobjval, optreq = optreq : np.dot(C,x) - FBAobjval*optreq + 0.0001

print 'Cobyla solution:'
print 'FBA Objective value:',FBAobjval
print 'Distance to experimental data:',compdistcomplete(FBAres,cobramodel)


allconstr2 = allconstr + [objcon]

print '\nMinimizing distance to experimental fluxes using COBYLA:'

print '\nx0 = 0, optreq = 0'
mindistsol = optimize.fmin_cobyla(compdistcomplete,x0,allconstr,args = (cobramodel,),consargs = ())
print 'Distance:',compdistcomplete(mindistsol,cobramodel)
print 'FBA objective value:',objectiveValue(mindistsol,C)

print '\nx0 = FBAres, optreq = 0'
mindistsol3 = optimize.fmin_cobyla(compdistcomplete,FBAres,allconstr,args = (cobramodel,),consargs = ())
print 'Distance:',compdistcomplete(mindistsol3,cobramodel)
print 'FBA objective value:',objectiveValue(mindistsol3,C)

print '\nx0 = 0, optreq = 1'
mindistsol2 = optimize.fmin_cobyla(compdistcomplete,x0,allconstr2,args = (cobramodel,),consargs = ())
print 'Distance:',compdistcomplete(mindistsol2,cobramodel)
print 'FBA objective value:',objectiveValue(mindistsol2,C)

print '\nx0 = FBAres, optreq = 1'
mindistsol4 = optimize.fmin_cobyla(compdistcomplete,FBAres,allconstr2,args = (cobramodel,),consargs = ())
print 'Distance:',compdistcomplete(mindistsol4,cobramodel)
print 'FBA objective value:',objectiveValue(mindistsol4,C)




if __name__ == "__main__": #If the module is executed as a program, run a test.
    from cobra.io.sbml import create_cobra_model_from_sbml_file
    cobramodel = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')
    import loadData as load
    fluxvalues = load.ExpFluxesfromXML('expdata.xml','Perrenoud','Batch','aerobe')
    rmap = load.ReactionMapfromXML('reactionmaps.xml','Perrenoud','SCHUETZR')




