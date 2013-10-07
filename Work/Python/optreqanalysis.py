#optreqanalysis.py
#Analyse how the minimal achievable distance varies with the relative optimality requirement:

from QPmindist import *
import scitools.easyviz as ev
import scipy.io

def optreqanalysis(models,objectives,experiments,steps = 10,makeplots = True, fluxvalues = "Default",reactionmap = "Default"):

    for objective in objectives:
        pass #Placeholder

    for experiment in experiments:
        pass #Placeholder

    if fluxvalues == "Default":
            expdata = scipy.io.loadmat('expdata.mat')
            perrenoud = expdata['expdata']['perrenoud']
            fluxvalarray = perrenoud[0][0][0][0][0][0][0][0][0][0][0][0][0][0]
            fluxvalues = [row[0] for row in fluxvalarray] #expdata.perrenoud.abs.batch.aerobe.fluxvalues


    if reactionmap == "Default":
            mat = scipy.io.loadmat('reactionmaps.mat')
            rmaps = mat['reactionmaps']

            Fmap = rmaps[0][0][0]
            Fmap2 = rmaps[0][0][3]

            reactionmap = Fmap2
        
    if ~(type(models) == type([])):
        models = [models]
    for model in models:
            cobramodel = model
            print 'cobramodel:',model
            print type(cobramodel)

            optreqs = []
            for i in range(steps+1):
                optreqs.append(i/steps)
            results = []
            for optreq in optreqs:
                result = QPmindist(cobramodel,fluxvalues,reactionmap,optreq)
                results.append(result)
            ev.plot(optreqs,results)
            ev.show()


if __name__ == "__main__":
    from cobra.io.sbml import create_cobra_model_from_sbml_file
    cobramodel = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')
    #print cobramodel
    #print type(cobramodel)
    optreqanalysis(cobramodel,[],[])
