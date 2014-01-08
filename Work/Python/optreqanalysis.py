#optreqanalysis.py
#Analyse how the minimal achievable distance varies with the relative optimality requirement:

from QPmindist import *
from compdist2 import *
from extractflux2 import *
import scitools.easyviz as ev
import scipy.io
import math
import prettyplotlib as ppl
# This is "import matplotlib.pyplot as plt" from the prettyplotlib library
from prettyplotlib import plt

def optreqanalysis(modeldicts,objectives,experiments,steps = 10,makeplots = True, expfluxdict = "Default",plotter = "easyviz"):
    '''Function for evaluating how the minimally achievable distance between a flux solution and a set of experimental fluxes
    changes when varying the required relative FBA objective-optimality of the flux solution.

    syntax: optreqanalysis(modeldicts,objectives,experiments,steps = 10,makeplots = True,expfluxdict)

    arguments:
    modeldicts: A python dictionary with the following fields:
        modelobject: A CobraPy model object for the model in question.
        reactionmap: A Cogs-FBA reaction map array mapping between model reactions and experimental reactions for the model and experiment in question.
    objectives: A list of objectives to be used in the analysis.
    experiments: A list of experiments to be used in the analysis
    steps: The granularity of the analysis and produced plots.
    makeplots: Whether plots should be produced or not
    expfluxdict: A python dictionary with experiment-ids as keys and the corresponding flux value arrays as values.
    '''
    for objective in objectives:
        pass #Placeholder

    for experiment in experiments:
        pass #Placeholder

    if expfluxdict == "Default": #This section should be rewritten to use an "expflux" dict.
            expdata = scipy.io.loadmat('expdata.mat') #load experimentaldata
            perrenoud = expdata['expdata']['perrenoud']
            fluxvalarray = perrenoud[0][0][0][0][0][0][0][0][0][0][0][0][0][0]
            fluxvalues = [row[0] for row in fluxvalarray] #expdata.perrenoud.abs.batch.aerobe.fluxvalues

    #if ~(type(models) == type([])):
     #   models = [models]
      #  print "TEST!"
    resultlist = []
    for modeldict in modeldicts:
            cobramodel = modeldict['modelobject']
            #print 'modeldict:',modeldict
            reactionmap = modeldict['reactionmap']
            #print 'cobramodel:',cobramodel
            #print type(cobramodel)

            optreqs = []
            for i in range(steps+1):
                optreqs.append(float(i)/steps)
            results = []
            #print 'opreqs:',optreqs
            #q = raw_input('Continue?')
            for optreq in optreqs:
                optimizedmodel = QPmindist(cobramodel,fluxvalues,reactionmap,optreq)
                results.append(math.sqrt(optimizedmodel.ObjVal))
            resultlist.append(results)
            #Using easyviz
            #ev.plot(optreqs,results)
            #ev.hold()
            
            #Using prettyplotlib
            #fig, ax = plt.subplots(1)
            #ppl.scatter(ax, x, y)
            #ppl.plot(optreqs,results)
    return resultlist
            
            
    #ev.show()
    #fig.savefig('optreqplot.png')

if __name__ == "__main__":
    from cobra.io.sbml import create_cobra_model_from_sbml_file
    SCHUETZR = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')
    ECME = create_cobra_model_from_sbml_file('../SBML/ECME.xml')
    iJO1366b = create_cobra_model_from_sbml_file('../SBML/iJO1366b.xml')
 
    mat = scipy.io.loadmat('reactionmaps.mat')
    rmaps = mat['reactionmaps']

    Fmap = rmaps[0][0][0]
    Cmap = rmaps[0][0][1]
    Gmap = rmaps[0][0][2]

    Fmap2 = rmaps[0][0][3]
    Cmap2 = rmaps[0][0][4]
    Gmap2 = rmaps[0][0][5]


    SCHUETZRdict = {}
    SCHUETZRdict['modelobject'] = SCHUETZR
    SCHUETZRdict['reactionmap'] = Fmap2
    ECMEdict = {}
    ECMEdict['modelobject'] = ECME
    ECMEdict['reactionmap'] =  Cmap2

    iJOdict = {}
    iJOdict['modelobject'] = iJO1366b
    iJOdict['reactionmap'] = Gmap2

    models = [SCHUETZRdict,iJOdict] #ECME model optimization fails.
    results = optreqanalysis(models,[],[])
