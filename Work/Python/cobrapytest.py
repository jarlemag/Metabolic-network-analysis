#cobrapytest.py

from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.io.sbml import write_cobra_model_to_sbml_file

from cobra.test import test_all

#test_all()

#Load models
ECME = create_cobra_model_from_sbml_file('../SBML/ECME.xml')

SCHUETZR = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')

iJO1366b = create_cobra_model_from_sbml_file('../SBML/iJO1366b.xml')

#Perform FBA:

ECME.optimize(solver='gurobi')

print ECME.solution.f


SCHUETZR.optimize(solver='gurobi')

print SCHUETZR.solution.f


iJO1366b.optimize(solver='gurobi')

print iJO1366b.solution.f
