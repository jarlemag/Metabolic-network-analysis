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
Y[:] = np.NAN


for entry in reactionmap:
    #print entry[0]
    Y[entry[1]-1] = fluxvalues[entry[0]-1]

print Y


gurobimodel =gurobi.Model("QP")
objective = gurobi.QuadExpr()
for reaction in cobramodel.reactions:
    newvar = gurobimodel.addVar(lb = reaction.lower_bound, ub = reaction.upper_bound, name = reaction.id)
    gurobimodel.update()
reactlist = []
terms = []
for i in range(len(Y)):
    if ~np.isnan(Y[i]):
        reactlist.append(cobramodel.reactions[i])
        x = gurobimodel.getVarByName(cobramodel.reactions[i].id)
        newterm = x * x - 2*Y[i]*x
        objective.add(newterm)
        gurobimodel.update()
        terms.append(newterm)
print reactlist
print (len(reactlist))

for metabolite in cobramodel.metabolites:
    reactions = metabolite.get_reaction()
    newconstr = gurobi.LinExpr([rxn.get_coefficient(metabolite.id) for rxn in reactions],[gurobimodel.getVarByName(rxn.id) for rxn in reactions])
    gurobimodel.addConstr(newconstr, gurobi.GRB.EQUAL, 0, metabolite.id)
gurobimodel.modelsense = 1 #Set model sense to minimize
gurobimodel.update()

gurobimodel.optimize()
gurobimodel.update()

for v in gurobimodel.getVars():
    print v.varName, v.x 



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
