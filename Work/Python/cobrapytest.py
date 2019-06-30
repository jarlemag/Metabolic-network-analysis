#cobrapytest.py

from cobra.io.sbml import read_sbml_model
from cobra.io.sbml import write_sbml_model

from cobra.test import test_all

#test_all()

#Load models
ECME = read_sbml_model('../SBML/ECME.xml')

SCHUETZR = read_sbml_model('../SBML/SCHUETZR.xml')

iJO1366b = read_sbml_model('../SBML/iJO1366b.xml')

#Perform FBA:

#ECME.optimize(solver='gurobi')

#print 'ECME:',ECME.solution.f


SCHUETZR.optimize()

#print 'SCHUETZR:',SCHUETZR.solution.f


#iJO1366b.optimize(solver='gurobi')

#print 'iJO1366b:',iJO1366b.solution.f

