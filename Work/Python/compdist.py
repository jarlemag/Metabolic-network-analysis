#compdist2.py
#Rewrite of compdist


import numpy as np
import scipy.io
import loadData as load
import extractflux

def compdistdict(fluxdict,expdata = None,options = None,sense = None):
    if expdata == None:
        expdata = load.ExpFluxesfromXML('expdata.xml','Perrenoud','Batch','aerobe')

    fluxvector = []
    fluxvalues = []

    
    for key in fluxdict:
        if (key in expdata):
            fluxvector.append(fluxdict[key])
            fluxvalues.append(float(expdata[key]))
 
    dist = np.linalg.norm(np.array(fluxvector)-np.array(fluxvalues))
    return dist

def compdist(fluxvector,expdata = None, options = None, sense = None, extract = True, rmap = 'Default'):
    rmap = load.ReactionMapfromXML('reactionmaps.xml','Perrenoud','SCHUETZR')
    if extract:
        fluxvector = extractflux.extractfluxvector(fluxvector,rmap)

    if expdata == None:
        expdata = load.ExpFluxesfromXML('expdata.xml','Perrenoud','Batch','aerobe', vector = True)
    dist = np.linalg.norm(np.array(fluxvector)-np.array(expdata))
    return dist




if __name__ == "__main__": #If the module is executed as a program, run a test.
    from cobra.io.sbml import read_sbml_model
    cobramodel = read_sbml_model('../SBML/SCHUETZR.xml')
    import loadData as load
    fluxvalues = load.ExpFluxesfromXML('expdata.xml','Perrenoud','Batch','aerobe')
    rmap = load.ReactionMapfromXML('reactionmaps.xml','Perrenoud','SCHUETZR')

    optreq = 1.0

    import QPmindist as qp

    gurobimodel = qp.QPmindist(cobramodel,fluxvalues,rmap,optreq)

    QPsolutiondict = qp.getgurobisolutiondict(gurobimodel)

    import extractflux2 as extract

    extractfluxdict = extract.extractfluxdict(QPsolutiondict,rmap)


    dictdist = compdistdict(extractfluxdict)

    print('dist:',dictdist)

    simple1 = {"A":1,"B":2,"C":3}
    simple2 = {"A":1,"B":2,"C":4}
    simpledist = compdistdict(simple1, expdata = simple2)
