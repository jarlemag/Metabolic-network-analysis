#computesplits.py

def computesplits(fluxdict,splitmap):
    '''Computes split ratios from a set of flux values and a definition of split ratios.
    '''
    output = {}
    for splitkey in splitmap:
        denomfluxes = 0
        numfluxes = 0
        for reaction in splitmap[splitkey]['denominator']:
            rxid = reaction['id']
            rawflux = fluxdict[rxid]
            if 'coef' in reaction:
                rawflux = rawflux*int(reaction['coef'])
            if 'modification' in reaction:
                if reaction['modification'] == 'subplus':
                    rawflux = max(0,rawflux)
            denomfluxes += rawflux
        if denomfluxes < 1e-4:
            denomfluxes = 0
        for reaction in splitmap[splitkey]['numerator']:
            rxid = reaction['id']
            rawflux = fluxdict[rxid]
            if 'coef' in reaction:
                rawflux = rawflux*int(reaction['coef'])
            if 'modification' in reaction:
                if reaction['modification'] == 'subplus':
                    rawflux = max(0,rawflux)
            numfluxes += rawflux
        if numfluxes < 1e-4:
            numfluxes = 0
        if denomfluxes == 0:
            output[splitkey] = 0 #To avoid divide by zero
        else:
            output[splitkey] = numfluxes/denomfluxes
    return output
            


def printsplits(splits):

    keys = splits.keys()
    keys.sort(key=int)
    for key in keys:
        print key,'%.2f' %splits[key]

if __name__ == "__main__":

    import loadData as load
    from cobra.io.sbml import create_cobra_model_from_sbml_file
    import extractflux as extract

    expfluxes = load.ExpFluxesfromXML('expdata.xml','Perrenoud','Batch','aerobe')

    splitmap = load.SplitMapfromXML('reactionmaps.xml','SCHUETZR','EXPDATA')

    splits = computesplits(expfluxes,splitmap)
    print 'Experimental splits:'
    printsplits(splits)
    splitmap = load.SplitMapfromXML('reactionmaps.xml','SCHUETZR','SCHUETZR')

    cobramodel = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')
    cobramodel.optimize(solver='gurobi')

    compfluxes = extract.vectorToDict(cobramodel.solution.x,cobramodel)
    
    compsplits = computesplits(compfluxes,splitmap)
    print 'SCHUETZR FBA result splits:'
    printsplits(compsplits)
