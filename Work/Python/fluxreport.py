#fluxreport.py

import numpy as np
import loadData as load
import copy


def fluxreport(fluxdict,expfluxdict = load.ExpFluxesfromXML('expdata.xml','Perrenoud','Batch','aerobe'), exp_errordict = load.ExpErrorsfromXML('expdata.xml','Perrenoud','Batch','aerobe')):

    fluxdictcopy = copy.deepcopy(fluxdict)
    for key in fluxdictcopy:
        if key not in expfluxdict:
            del fluxdict[key]

    diff_dict = {key:(expfluxdict[key] - fluxdict[key]) for key in fluxdict}
    #bug here: KeyError: '26'
    

    #uncertfraction = difference / experrors
    uncertfrac_dict = {key:(np.divide(diff_dict[key],exp_errordict[key])) for key in fluxdict}
    
    textheader = ["Ex rxn #", "Com. flux","Exp flux", "Diff.","Exp. uncert.", "Diff/uncert."]
    
    outsidefluxes = 0
    for key in uncertfrac_dict:
        if (abs(uncertfrac_dict[key]) > 1) and (~np.isinf(abs(uncertfrac_dict[key]))):
            outsidefluxes +=1

    import datetime
    now = datetime.datetime.now()
    print "Flux report generated at:",now
    print "Solver:"
    print "Model:"
    print "Objective:"
    print "Constraints:"
    print "Number of computed fluxes outside experimental uncertainty bounds:",outsidefluxes


    #print 'fluxdict:',fluxdict #DEBUG
    #print 'expfluxdict:',expfluxdict #DEBUG

    print "%15s %15s %15s %15s %15s %15s" % tuple(textheader)

    
    for key in fluxdict:
        #print key #DEBUG
        #print fluxdict[key]
        #print expfluxdict[key]
        #print diff_dict[key]
        #print exp_errordict[key]
        #print uncertfrac_dict[key]
        #print uncertfrac_dict
        print "%s %15.1f %15.1f %15.1f %15.1f %15.1f" % (key,float(fluxdict[key]),float(expfluxdict[key]),float(diff_dict[key]),float(exp_errordict[key]),float(uncertfrac_dict[key]))

    return diff_dict
    

#TEST CODE START:
if __name__ == "__main__":
    from cobra.io.sbml import create_cobra_model_from_sbml_file

    SCHUETZR = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')
    SCHUETZR.optimize(solver='gurobi')

    FBAresult = SCHUETZR.solution.x_dict
    
    import extractflux2

    rmap = load.ReactionMapfromXML('reactionmaps.xml','Perrenoud','SCHUETZR')
    expfluxes = load.ExpFluxesfromXML('expdata.xml','Perrenoud','Batch','aerobe')

    fluxresult = extractflux2.extractfluxdict(FBAresult,rmap)


    t = fluxreport(fluxresult)

#TEST CODE END.


