#testcompdist2.py
#test script for compdist2

from compdist2 import *

from extractflux2 import *

from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.io.sbml import write_cobra_model_to_sbml_file

model = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')

model.optimize(solver='gurobi')

import scipy.io

mat = scipy.io.loadmat('reactionmaps.mat')
rmaps = mat['reactionmaps']

Fmap = rmaps[0][0][0]
Cmap = rmaps[0][0][1]
Gmap = rmaps[0][0][2]

Fmap2 = rmaps[0][0][3]
Cmap2 = rmaps[0][0][4]
Gmap2 = rmaps[0][0][5]


extractedfluxes = extractflux2(model.solution.x,Fmap2)

distance = compdist2(extractedfluxes)
print 'distance:',distance


#Add code here to verify that output is correct:
