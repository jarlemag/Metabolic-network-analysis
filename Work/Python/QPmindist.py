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
        return objval

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


    #Create a new Gurobi model
    gurobimodel =gurobi.Model("QP")
    #Create a new QuadExpr object to store the objective function.
    QPobjective = gurobi.QuadExpr()

    #Add decision variables (reaction fluxes) to the Gurobi model, and construct the FBA objective function:

    FBAobjective = gurobi.LinExpr() #Create a new Gurobi linear expression for the FBA objective function
    #For every reaction in the Cobra model:
    for reaction in cobramodel.reactions:
        #Create a new decision variable in the Gurobi model, with upper and lower bounds as specified in the Cobra model:
        newvar = gurobimodel.addVar(lb = reaction.lower_bound, ub = reaction.upper_bound, name = reaction.id)
        gurobimodel.update()
        FBAobjective.add(reaction.objective_coefficient * newvar) #Construct the FBA objective
        gurobimodel.update()
        if debug:
            print FBAobjective
            z = raw_input('Press enter to continue.')
        #Update the Gurobi model
        

    #Add steady-state constraints for every metabolite:
    for metabolite in cobramodel.metabolites:
        #Get the list of reactions in which the metabolite partakes:
        reactions = metabolite.get_reaction()
        newconstr = gurobi.LinExpr([rxn.get_coefficient(metabolite.id) for rxn in reactions],[gurobimodel.getVarByName(rxn.id) for rxn in reactions])
        gurobimodel.addConstr(newconstr, gurobi.GRB.EQUAL, 0, metabolite.id)
    gurobimodel.update()

    #Perform FBA to determine the target objective function value to inform the optimality requirement constraint:

    if debug:
        print 'FBAobjective variable:',FBAobjective
        print 'gurobimodel objective before update:',str(gurobimodel.getObjective())
    gurobimodel.setObjective(FBAobjective, gurobi.GRB.MAXIMIZE)
    gurobimodel.update()
    if debug:
        print 'gurobimodel objective after update:',str(gurobimodel.getObjective())
    gurobimodel.optimize()

    FBAsolution = [v.x for v in gurobimodel.getVars()]
    FBAobjval = gurobimodel.Objval
    if debug:
        print "FBA objective value:",FBAobjval

    #Add optimality requirement constraint to original model:
    if useoptreq:
        gurobimodel.addConstr(FBAobjective, gurobi.GRB.GREATER_EQUAL,FBAobjval*optreq)
        gurobimodel.update()    

    #Construct the QP objective:

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
    from cobra.io.sbml import create_cobra_model_from_sbml_file
    cobramodel = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')

    import loadData as load
    rmap = load.ReactionMapfromXML('reactionmaps.xml','Perrenoud','SCHUETZR')
    
    optreq = 1.0

    fluxvalues = load.ExpFluxesfromXML('expdata.xml','Perrenoud','Batch','aerobe')

    gurobimodel = QPmindist(cobramodel,fluxvalues,rmap,optreq)
    QPsolutionvector = getgurobisolutionvector(gurobimodel)
    QPsolutiondict = getgurobisolutiondict(gurobimodel)
    QPFBAobjval = computeFBAobjval(QPsolutionvector,cobramodel)

    import extractflux2
    extractfluxdict = extractflux2.extractfluxdict(QPsolutiondict,rmap)

    import fluxreport
    #fluxreport.fluxreport(extractedfluxes,fluxvalarray)

    import compdist2

    #dist = compdist2.compdist2(extractedfluxes)
    #print 'Optimality requirement:',optreq
    #print 'compdist2 distance:',dist

    #Sanity check:
    # terms[1] = -10.8 pgi -> experimental pgi flux should be 5.4
    #cobramodel.reactions[12] = pgi #pgi reaction is reaction # 13 in model
    #Look in reactionmap, the reaction corresponding to model reaction 13 is experimental reaction 4.
    #Check experimental fluxvalues. fluxvalues[3] = 5.4 --> OK!

    #terms[2] = -12.4 fbaAb, experimental fbaAb flux should be 6.2
    #How to find index of a certain element:
    #f = [i for i in range(len(cobramodel.reactions)) if cobramodel.reactions[i].id == "fbaAB"]
    #cobramodel.reaction[14]= fbaAB -> fbaAB is reaction #15 in models
    #Look in reactionmap, the reaction corresponding to modelreaction  is experimental reaction 6
    #Check experimental fluxvalues. fluxvalues[5] = 6.2 -> OK!

    print 'gurobimomdel gurobi objective value:',gurobimodel.ObjVal
    import math
    print 'square root of gurobi objective value:',math.sqrt(gurobimodel.ObjVal)
    print "QP solution FBA objective value:",QPFBAobjval
