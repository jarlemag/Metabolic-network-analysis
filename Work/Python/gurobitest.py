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
