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
    fluxvector = []
    fluxvalues = []

    for key in extractfluxdict:
        if (key in expdata):
            fluxvector.append(extractfluxdict[key])
            fluxvalues.append(float(expdata[key]))
 
    dist = np.linalg.norm(np.array(fluxvector)-np.array(fluxvalues))

    return dist


class CobylaModel:
    def __init__(self,cobramodel, tol = 0.001):
        self.cobramodel = cobramodel.to_array_based_model()
        self.S = cobramodel.S.toarray()
        self.ub = self.cobramodel.upper_bounds
        self.lb = self.cobramodel.lower_bounds
        self.C = self.cobramodel.objective_coefficients
        self.ub_funcs = constrainFunctionsUB(self.ub)
        self.lb_funcs = constrainFunctionsLB(self.lb)
        self.ss_funcs_a = []
        self.ss_funcs_b = []
        for row in self.S:
            self.ss_funcs_a.append((lambda x, row = row: tol + sum(np.multiply(row,x))))
            self.ss_funcs_b.append((lambda x, row = row: tol -sum(np.multiply(row,x))))
        self.allconstr = self.ub_funcs + self.lb_funcs + self.ss_funcs_a + self.ss_funcs_b
        
        
def cobylaFBA(cobramodel, tol = 0.001):
    cobmodel = CobylaModel(cobramodel, tol = tol) 
    x0 = [0 for element in cobmodel.C]
    FBAres = optimize.fmin_cobyla(FBAobjective,x0,cobmodel.allconstr,args = (cobmodel.C,),consargs =())
    return FBAres
    
def mindist(cobramodel,x0,FBAobjval, optreq = 0, tol = 0.001):
    cobmodel = CobylaModel(cobramodel, tol = tol)
    if x0 is None:
        x0 = [0 for element in cobmodel.C]
    if optreq == 0:
        minsol = optimize.fmin_cobyla(compdistcomplete,x0,cobmodel.allconstr,args = (cobramodel,),consargs = ())
    else:
        objcon = lambda x, objective = cobmodel.C , value = FBAobjval, optreq = optreq : np.dot(C,x) - FBAobjval*optreq #+ 0.0001
        allconstr2 = cobmodel.allconstr + [objcon]
        minsol = optimize.fmin_cobyla(compdistcomplete,x0,allconstr2,args = (cobramodel,),consargs = ())
        
    return minsol

if __name__ == "__main__": #If the module is executed as a program, run a test.
    from cobra.io.sbml import create_cobra_model_from_sbml_file
    cobramodel = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')
    cobramodel.optimize(solver='gurobi')
    CobraPyobjectivevalue = cobramodel.solution.f
    C = cobramodel.to_array_based_model().objective_coefficients
    import loadData as load
    fluxvalues = load.ExpFluxesfromXML('expdata.xml','Perrenoud','Batch','aerobe')
    rmap = load.ReactionMapfromXML('reactionmaps.xml','Perrenoud','SCHUETZR')

    print 'Optimizing FBA objective with COBYLA.'

    CobylaFBAres = cobylaFBA(cobramodel)
    CobylaFBAobjectivevalue = FBAobjective(CobylaFBAres,C, sense = 1)
    
    print 'Cobyla solution:'
    print 'FBA Objective value:',CobylaFBAobjectivevalue
    print 'Distance to experimental data:',compdistcomplete(CobylaFBAres,cobramodel)

    print '\nMinimizing distance to experimental fluxes using COBYLA:'

    print 'Using COBYLA-calculated FBA objective value as target:'
    print '\nx0 = 0, optreq = 0'
    mindistsol = mindist(cobramodel,None,CobylaFBAobjectivevalue)
    print 'Distance:',compdistcomplete(mindistsol,cobramodel)
    print 'FBA objective value:',objectiveValue(mindistsol,C)

    print '\nx0 = FBAres, optreq = 0'
    mindistsol2 = mindist(cobramodel,CobylaFBAres,CobylaFBAobjectivevalue)
    print 'Distance:',compdistcomplete(mindistsol2,cobramodel)
    print 'FBA objective value:',objectiveValue(mindistsol2,C)

    print '\nx0 = 0, optreq = 1'
    mindistsol3 =  mindist(cobramodel,None,CobylaFBAobjectivevalue,optreq = 1)
    print 'Distance:',compdistcomplete(mindistsol3,cobramodel)
    print 'FBA objective value:',objectiveValue(mindistsol3,C)

    print '\nx0 = FBAres, optreq = 1'
    mindistsol4 = mindist(cobramodel,CobylaFBAres,CobylaFBAobjectivevalue,optreq = 1)
    print 'Distance:',compdistcomplete(mindistsol4,cobramodel)
    print 'FBA objective value:',objectiveValue(mindistsol4,C)

    print 'Using CobraPy-calculated objective value as target:'
    print '\nx0 = 0, optreq = 0'
    mindistsol5 = mindist(cobramodel,None,CobraPyobjectivevalue)
    print 'Distance:',compdistcomplete(mindistsol5,cobramodel)
    print 'FBA objective value:',objectiveValue(mindistsol5,C)

    print '\nx0 = FBAres, optreq = 0'
    mindistsol6 = mindist(cobramodel,CobylaFBAres,CobraPyobjectivevalue)
    print 'Distance:',compdistcomplete(mindistsol6,cobramodel)
    print 'FBA objective value:',objectiveValue(mindistsol6,C)

    print '\nx0 = 0, optreq = 1'
    mindistsol7 =  mindist(cobramodel,None,CobraPyobjectivevalue,optreq = 1)
    print 'Distance:',compdistcomplete(mindistsol7,cobramodel)
    print 'FBA objective value:',objectiveValue(mindistsol7,C)

    print '\nx0 = FBAres, optreq = 1'
    mindistsol8 = mindist(cobramodel,CobylaFBAres,CobraPyobjectivevalue,optreq = 1)
    print 'Distance:',compdistcomplete(mindistsol8,cobramodel)
    print 'FBA objective value:',objectiveValue(mindistsol8,C)

