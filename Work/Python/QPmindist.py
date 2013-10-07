#QPmindist.py
#Reimplementation of QPmindist.m

import numpy as np
import scipy.io
import gurobipy as gurobi


def getgurobisolution(model):
        sol = [v.x for v in model.getVars()]
        return sol

def computeFBAobjval(fluxsolution,model):
        #Assumes that reactions are in the same order in the fluxsolution vector and the model reactions list!
        objval = 0
        for i in range(len(fluxsolution)):
            objval += fluxsolution[i]*model.reactions[i].objective_coefficient
        return objval

def QPmindist(cobramodel,fluxvalues,reactionmap,optreq,useoptreq = True,debug = False):
   
    n = len(cobramodel.reactions)
    Y = np.zeros(n)
    #Create a vector to hold the experimental flux values:
    Y[:] = np.NAN

    #For every entry in the reaction map:
    for entry in reactionmap:
        #print entry[0]
        #Update the experimental flux values vector based on the entries in the reaction map and the flux values vector
        Y[entry[1]-1] = fluxvalues[entry[0]-1]

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
    #For every element in the vector of experimental flux values:
    for i in range(len(Y)):
        #If an experimental flux value is given (the value is not "not a number")
        if ~np.isnan(Y[i]):
            reactlist.append(cobramodel.reactions[i]) #For debugging/verbosity
            #Set the current decision variable to the corresponding model reaction
            x = gurobimodel.getVarByName(cobramodel.reactions[i].id)
            #Create an objective function term for the reaction:
            newterm = x * x - 2*Y[i]*x
            newterm.addConstant((Y[i]**2))
            #Note: Adding constant terms to the objective changes the computed solution slightly in at least some cases.
            #Add the term to the objective function:
            QPobjective.add(newterm)
            #Apply the objective function to the Gurobi model:
            gurobimodel.setObjective(QPobjective, gurobi.GRB.MINIMIZE)
            #Update the Gurobi model
            gurobimodel.update()
            terms.append(newterm)
    if debug:
        print reactlist 
        print (len(reactlist)) 
    #Perform QP optimization:

    gurobimodel.optimize()
    gurobimodel.update()

    return gurobimodel


if __name__ == "__main__": #If the module is executed as a program, run a test.
    from cobra.io.sbml import create_cobra_model_from_sbml_file
    cobramodel = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')

    mat = scipy.io.loadmat('reactionmaps.mat')
    rmaps = mat['reactionmaps']

    Fmap = rmaps[0][0][0]
    Fmap2 = rmaps[0][0][3]

    reactionmap = Fmap2
    '''
    TO DO: Change reaction maps to work solely with reaction IDs (Mapping between experimental reaction IDs and model reaction IDs.
    Stop messing around with numerical indexes.
    '''

    optreq = 0.8

    expdata = scipy.io.loadmat('expdata.mat')
    perrenoud = expdata['expdata']['perrenoud']
    fluxvalarray = perrenoud[0][0][0][0][0][0][0][0][0][0][0][0][0][0]
    fluxvalues = [row[0] for row in fluxvalarray] #expdata.perrenoud.abs.batch.aerobe.fluxvalues

    gurobimodel = QPmindist(cobramodel,fluxvalues,reactionmap,optreq)
    QPsolution = getgurobisolution(gurobimodel)
    QPFBAobjval = computeFBAobjval(QPsolution,cobramodel)

    import extractflux2
    extractedfluxes = extractflux2.extractflux(QPsolution,reactionmap)

    import fluxreport
    fluxreport.fluxreport(extractedfluxes,fluxvalarray)

    import compdist2

    dist = compdist2.compdist2(extractedfluxes)
    print 'Optimality requirement:',optreq
    print 'compdist2 distance:',dist

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
