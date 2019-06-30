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

ECMEsolution = ECME.optimize()

print ('ECME status:',ECMEsolution.status)
print ('ECME:',ECMEsolution.objective_value)



SCHUETZRsolution = SCHUETZR.optimize()

print ('SCHUETZR status:',SCHUETZRsolution.status)
print ('SCHUETZR objective value:',SCHUETZRsolution.objective_value)


iJO1366bsolution = iJO1366b.optimize()

print ('iJO1366b status:',iJO1366bsolution.status)
print ('iJO1366b objective value:',iJO1366bsolution.objective_value)
