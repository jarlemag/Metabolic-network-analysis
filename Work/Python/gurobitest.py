#gurobitest.py
#Testing of gurobi
import gurobipy as gurobi


m = gurobi.Model("mip1")

x = m.addVar(vtype=gurobi.GRB.BINARY, name = "x")
y = m.addVar(vtype=gurobi.GRB.BINARY, name = "y")
z = m.addVar(vtype=gurobi.GRB.BINARY, name = "z")
m.update()
m.setObjective(x + y + 2 * z, gurobi.GRB.MAXIMIZE)
m.addConstr(x + 2 * y + 3 * z <=4, "c0")
m.addConstr(x + y >= 1,"c1")
m.update()
m.optimize()
for v in m.getVars():
    print v.varName, v.x 

print 'Obj:', m.objVal



#Use gurobi to solve an FBA problem:

from cobra.io.sbml import create_cobra_model_from_sbml_file

cobramodel = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')


gurobimodel =gurobi.Model("Schuetz") #Initialize a new Gurobi model


#add one decision variable for each reaction in the cobrapy model:

objective = gurobi.LinExpr() 
variables = []
for reaction in cobramodel.reactions:
    newvar = gurobimodel.addVar(lb = reaction.lower_bound, ub = reaction.upper_bound, name = reaction.id) #Add reactions as decision variables, with bounds.
    objective.add(reaction.objective_coefficient * newvar) #Add the reactions to the objective with their respective objective coefficients.
    #print reaction.id
    #variables += newvar
gurobimodel.update()
gurobimodel.setObjective(objective) #Apply the objective defined above.
gurobimodel.update()

gurobimodel.getObjective() #Show the objective

for metabolite in cobramodel.metabolites:
    reactions = metabolite.get_reaction()
    '''
    #newconstr = gurobi.LinExpr()
    for rxn in reactions:
        coef = rxn.get_coefficient(metabolite.id)
        #newconstr.add(coef* gurobimodel.getVarByName(rxn.id))
    '''
    newconstr = gurobi.LinExpr([rxn.get_coefficient(metabolite.id) for rxn in reactions],[gurobimodel.getVarByName(rxn.id) for rxn in reactions])
    gurobimodel.addConstr(newconstr, gurobi.GRB.EQUAL, 0, metabolite.id)    
        #print coef
gurobimodel.update()

#Set sense to maximize:
gurobimodel.modelsense = -1
gurobimodel.update()

gurobimodel.optimize()


gurobimodel.NumVars
gurobimodel.NumConstrs

f = gurobimodel.optimize()
