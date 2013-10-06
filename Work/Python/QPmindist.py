#QPmindist.py
#Reimplementation of QPmindist.m

import numpy as np
import scipy.io

#test code START. Move to separate test script later.
#def QPmindist(model):
from cobra.io.sbml import create_cobra_model_from_sbml_file
cobramodel = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')

import gurobipy as gurobi

mat = scipy.io.loadmat('reactionmaps.mat')
rmaps = mat['reactionmaps']

Fmap = rmaps[0][0][0]
Fmap2 = rmaps[0][0][3]

reactionmap = Fmap2
'''
TO DO: Change reaction maps to work solely with reaction IDs (Mapping between experimental reaction IDs and model reaction IDs.
Stop messing around with numerical indexes.
'''

#test code END


expdata = scipy.io.loadmat('expdata.mat')
perrenoud = expdata['expdata']['perrenoud']
fluxvalarray = perrenoud[0][0][0][0][0][0][0][0][0][0][0][0][0][0]
fluxvalues = [row[0] for row in fluxvalarray] #expdata.perrenoud.abs.batch.aerobe.fluxvalues

n = len(cobramodel.reactions)
Y = np.zeros(n)
#Create a vector to hold the experimental flux values:
Y[:] = np.NAN


#For every entry in the reaction map:
for entry in reactionmap:
    #print entry[0]
    #Update the experimental flux values vector based on the entries in the reaction map and the flux values vector
    Y[entry[1]-1] = fluxvalues[entry[0]-1]

print Y

#Create a new Gurobi model
gurobimodel =gurobi.Model("QP")
#Create a new QuadExpr object to store the objective function.
objective = gurobi.QuadExpr()

#For every reaction in the Cobra model:
for reaction in cobramodel.reactions:
    #Create a new decision variable in the Gurobi model, with upper and lower bounds as specified in the Cobra model:
    newvar = gurobimodel.addVar(lb = reaction.lower_bound, ub = reaction.upper_bound, name = reaction.id)
    #Update the Gurobi model
    gurobimodel.update()
reactlist = []
terms = []

#Construct the objective:

#For every element in the vector of experimental flux values:
for i in range(len(Y)):
    #If an experimental flux value is given (the value is not "not a number")
    if ~np.isnan(Y[i]):
        reactlist.append(cobramodel.reactions[i]) #For debugging/verbosity
        #Set the current decision variable to the corresponding model reaction
        x = gurobimodel.getVarByName(cobramodel.reactions[i].id)
        #Create an objective function term for the reaction:
        newterm = x * x - 2*Y[i]*x
        #linterm = gurobi.LinExpr(float(Y[i]**2)) #Y[i]**2 is numpy float. Must convert to regular float before passing to LinExpr.
        newterm.addConstant((Y[i]**2)) #More succint than the above line
        #Note: Adding constant terms to the objective changes the computed solution slightly in at least some cases.
        #Add the term to the objective function:
        objective.add(newterm)
        #Apply the objective function to the Gurobi model:
        gurobimodel.setObjective(objective)
        #Update the Gurobi model
        gurobimodel.update()
        terms.append(newterm)
print reactlist #For debugging/verbosity
print (len(reactlist)) #For debugging/verbosity

#Check the objective:
theobjective = gurobimodel.getObjective()
print str(theobjective)

#For every metabolite:
for metabolite in cobramodel.metabolites:
    #Get the list of reactions in which the metabolite partakes:
    reactions = metabolite.get_reaction()
    #Make a new constraint in the Gurobi model, such that the fluxes of the metabolite are balanced:
    newconstr = gurobi.LinExpr([rxn.get_coefficient(metabolite.id) for rxn in reactions],[gurobimodel.getVarByName(rxn.id) for rxn in reactions])
    #Add the new constraint to the Gurobi model
    gurobimodel.addConstr(newconstr, gurobi.GRB.EQUAL, 0, metabolite.id)
gurobimodel.modelsense = 1 #Set model sense to minimize
gurobimodel.update()

gurobimodel.optimize()
gurobimodel.update()

'''
for v in gurobimodel.getVars():
    print v.varName, v.x 
'''

solvec =[v.getAttr("x") for v in gurobimodel.getVars()]

import extractflux2

extractedfluxes = extractflux2.extractflux(solvec,Fmap)

import fluxreport

fluxreport.fluxreport(extractedfluxes,fluxvalarray)

import compdist2

dist = compdist2.compdist2(extractedfluxes)
print 'dist:',dist

#Sanity check:
# terms[1] = -10.8 pgi -> experimental pgi flux should be 5.4
#cobramodel.reactions[12] = pgi #pgi reaction is reaction # 13 in model
#Look in reactionmap, the reaction corresponding to model reaction 13 is experimental reaction 4.
#Check experimental fluxvalues. fluxvalues[3] = 5.4 --> OK!

#terms[2] = -12.4 fbaAb, experimental fbaAb flux should be 6.2
#How to find index of a certain element:
f = [i for i in range(len(cobramodel.reactions)) if cobramodel.reactions[i].id == "fbaAB"]
#cobramodel.reaction[14]= fbaAB -> fbaAB is reaction #15 in models
#Look in reactionmap, the reaction corresponding to modelreaction  is experimental reaction 6
#Check experimental fluxvalues. fluxvalues[5] = 6.2 -> OK!


print gurobimodel.ObjVal
