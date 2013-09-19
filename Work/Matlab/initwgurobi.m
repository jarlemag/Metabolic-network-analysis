%Initialize COBRA toolbox and set solver to Gurobi5

initCobraToolbox

setLPsolvertogurobi5 = changeCobraSolver('gurobi5','LP')
setMILPsolvertogurobi5 = changeCobraSolver('gurobi5','MILP')
setMIQPsolvertogurobi5 = changeCobraSolver('gurobi5','MIQP')
setQPsolvertogurobi5 = changeCobraSolver('gurobi5','QP')