#fluxreport.py

import numpy as np

def fluxreport(fluxresult,fluxvalarray):

    expfluxes = fluxvalarray[:,0]
    experrors = fluxvalarray[:,1]

    difference = expfluxes - fluxresult

    uncertfraction = difference / experrors

    number = range(len(fluxresult))

    textheader = ["Ex rxn #", "Com. flux","Exp flux", "Diff.","Exp. uncert.", "Diff/uncert."]
    
    outsidefluxes = 0
    for fraction in uncertfraction:
        if (abs(fraction) > 1) and (~np.isinf(abs(fraction))):
            outsidefluxes +=1

    import datetime
    now = datetime.datetime.now()
    print "Flux report generated at:",now
    print "Solver:"
    print "Model:"
    print "Objective:"
    print "Constraints:"
    print "Number of computed fluxes outside experimental uncertainty bounds:",outsidefluxes

    print "%15s %15s %15s %15s %15s %15s" % tuple(textheader)
    for i in range(len(expfluxes)):
        print "%15.1f %15.1f %15.1f %15.1f %15.1f %15.1f" % (number[i],fluxresult[i],expfluxes[i],difference[i],experrors[i],uncertfraction[i])
    
    #for 

    return difference
    




#TEST CODE START:
if __name__ == "__main__":
    from cobra.io.sbml import create_cobra_model_from_sbml_file

    SCHUETZR = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')
    SCHUETZR.optimize(solver='gurobi')

    FBAresult = SCHUETZR.solution.x

    import extractflux2

    import scipy.io

    mat = scipy.io.loadmat('reactionmaps.mat')
    rmaps = mat['reactionmaps']

    Fmap = rmaps[0][0][0]

    fluxresult = extractflux2.extractflux(FBAresult,Fmap)


    expdata = scipy.io.loadmat('expdata.mat')
    perrenoud = expdata['expdata']['perrenoud']
    fluxvalarray = perrenoud[0][0][0][0][0][0][0][0][0][0][0][0][0][0]


    t = fluxreport(fluxresult,fluxvalarray)

#TEST CODE END.


