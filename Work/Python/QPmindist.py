#QPmindist.py
#Reimplementation of QPmindist.m

import numpy as np
import scipy.io
import gurobipy as gurobi


def getgurobisolutionvector(model):
        sol = [v.x for v in model.getVars()]
        return sol

def getgurobisolutiondict(model):
        sol = {v.Varname:v.x for v in model.getVars()}
        return sol

def computeFBAobjval(fluxsolution,model):
        #Assumes that reactions are in the same order in the fluxsolution vector and the model reactions list!
        objval = 0
        for i in range(len(fluxsolution)):
            objval += fluxsolution[i]*model.reactions[i].objective_coefficient
            #print i #DEBUG
        return objval



def makeFBAmodel(cobramodel):
        '''
        Creates a gurobi model from a CobraPy model
        '''
        gurobimodel = gurobi.Model('QP')
        gurobimodel.setAttr ('ModelSense',-1) 
        for reaction in cobramodel.reactions:
                newvar = gurobimodel.addVar(lb = reaction.lower_bound, ub = reaction.upper_bound, name = reaction.id, obj = reaction.objective_coefficient)
                gurobimodel.update()
        for metabolite in cobramodel.metabolites:
                reactions = metabolite.get_reaction()
                newconstr = gurobi.LinExpr([rxn.get_coefficient(metabolite.id) for rxn in reactions],[gurobimodel.getVarByName(rxn.id) for rxn in reactions])
                gurobimodel.addConstr(newconstr, gurobi.GRB.EQUAL, 0, metabolite.id)
        return gurobimodel


def setFBAobjective(gurobimodel,cobramodel):
        FBAobjective = gurobi.LinExpr()
        for variable in gurobimodel.getVars():
                rxid = variable.getAttr('name')
                coef = cobramodel.reactions.get_by_id(rxid).objective_coefficient
                FBAobjective.add(variable * coef)
        gurobimodel.setObjective(FBAobjective, gurobi.GRB.MAXIMIZE)
        gurobimodel.update()
        return gurobimodel

def getFBAobjective(gurobimodel):
        FBAobjective = gurobi.LinExpr()
        for variable in gurobimodel.getVars():
                coef = variable.Obj
                FBAobjective.add(variable * coef)
        return FBAobjective

def gurobiFBA(cobramodel):
        gurobimodel = makeFBAmodel(cobramodel)
        #gurobimodel = setFBAobjective(gurobimodel)
        gurobimodel.update()
        gurobimodel.optimize()
        return gurobimodel

def QPmindist(cobramodel,fluxvalues,reactionmap,optreq,useoptreq = True,debug = False):
   
    #Create a dictionary to hold the experimental flux values:
    Ydict = {}

    #For every entry in the reaction map:
    for linkdict in reactionmap:
            if ((len(linkdict['modrxns']) > 1) or (len(linkdict['exprxns']) >1)): 
                    raise Exception('Model to experimental reaction mapping must be one to one.')
            modrxnid = linkdict['modrxns'][0]['rxid']
            exprxnid = linkdict['exprxns'][0]['rxid']
            modcoef = linkdict['modrxns'][0]['coef']
            #print(linkdict)
            if (exprxnid in fluxvalues):
                    Ydict[modrxnid] = fluxvalues[exprxnid]*modcoef


    
    gurobimodel = gurobiFBA(cobramodel)
    FBAsolution = [v.x for v in gurobimodel.getVars()]
    FBAobjval = gurobimodel.Objval
    FBAobjective = getFBAobjective(gurobimodel)
    if debug:
        print("FBA objective value:",FBAobjval)

    #Add optimality requirement constraint to original model:
    if useoptreq:
        gurobimodel.addConstr(FBAobjective, gurobi.GRB.GREATER_EQUAL,FBAobjval*optreq)
        gurobimodel.update()    

    #Create a new QuadExpr object to store the objective function.
    QPobjective = gurobi.QuadExpr()

    reactlist = []
    terms = []

    for key in Ydict:
            x = gurobimodel.getVarByName(key)
            fluxval = float(Ydict[key])
            newterm = x * x - 2*fluxval*x
            newterm.addConstant((fluxval**2))
            QPobjective.add(newterm)
            gurobimodel.setObjective(QPobjective, gurobi.GRB.MINIMIZE)
            gurobimodel.update()
            terms.append(newterm)
            
    #Perform QP optimization:

    gurobimodel.optimize()
    gurobimodel.update()
    return gurobimodel


if __name__ == "__main__": #If the module is executed as a program, run a test.
    from cobra.io.sbml import read_sbml_model
    cobramodel = read_sbml_model('../SBML/SCHUETZR.xml')

    import loadData as load
    rmap = load.ReactionMapfromXML('reactionmaps.xml','Perrenoud','SCHUETZR')
    
    optreq = 0.0

    fluxvalues = load.ExpFluxesfromXML('expdata.xml','Perrenoud','Batch','aerobe')

    gurobimodel = QPmindist(cobramodel,fluxvalues,rmap,optreq)
    QPsolutionvector = getgurobisolutionvector(gurobimodel)
    QPsolutiondict = getgurobisolutiondict(gurobimodel)
    QPFBAobjval = computeFBAobjval(QPsolutionvector,cobramodel)

    import extractflux
    extractfluxdict = extractflux.extractfluxdict(QPsolutiondict,rmap)

    import fluxreport
    #fluxreport.fluxreport(extractedfluxes,fluxvalarray)

    import compdist

    dist = compdist.compdistdict(extractfluxdict)
    #print 'Optimality requirement:',optreq
    print('Optreq = 0:')
    print('compdist distance:',dist)

    #print 'gurobimomdel gurobi objective value:',gurobimodel.ObjVal
    import math
    #print 'square root of gurobi objective value:',math.sqrt(gurobimodel.ObjVal)
    print("QP solution FBA objective value:",QPFBAobjval)

    #print 'Performing FBA only:'
    #QPFBAres = gurobiFBA(cobramodel)

    optreq = 1
    print('Optreq = 1:')
    gurobimodel2 = QPmindist(cobramodel,fluxvalues,rmap,optreq)
    QPsolutionvector2 = getgurobisolutionvector(gurobimodel2)
    QPFBAobjval2 = computeFBAobjval(QPsolutionvector2,cobramodel)
    QPsolutiondict2 = getgurobisolutiondict(gurobimodel2)
    extractfluxdict2 = extractflux.extractfluxdict(QPsolutiondict2,rmap)
    dist2 = compdist.compdistdict(extractfluxdict2)
    print('compdist distance:',dist2)
    print("QP solution FBA objective value:",QPFBAobjval2)
    
