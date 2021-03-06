import numpy as np

def extractfluxvector(fluxvector,reactionmap):
    
    fluxes = []

    for row in reactionmap:
        rxns = row[1:]
        totflux = 0
        for rxn in rxns:
            totflux += np.sign(rxn)*fluxvector[rxn-1]
        fluxes.append(totflux)
    return fluxes


def vectorToDict(fluxvector,model):
    tup = zip([reaction.id for reaction in model.reactions],fluxvector,)
    fluxdict = {reactionid:fluxvalue for reactionid,fluxvalue in tup}
    return fluxdict


def dictToVector(fluxdict,model):
    fluxvector = [fluxdict[reaction.id] for reaction in model.reactions]
    return fluxvector

def dictsToVectors(dict1,dict2):
    return [(dict1[key],dict2[key]) for key in dict1]

def extractfluxdict(fluxdict_in,reactionmap):
    fluxdict_out = {}
    for linkdict in reactionmap:
        if len(linkdict['exprxns']) > 1:
            raise Exception("Reaction map links one or more model reactions to more than one experimental reactions. A single model reaction or groups of model reactions may only map to one experimental reaction.")
        #print('linkdict:',linkdict)
        expid = linkdict['exprxns'][0]['rxid']
        totflux = 0
        for rxndict in linkdict['modrxns']:
            rxid = rxndict['rxid']
            coef = rxndict['coef']
            #print('rxid:',rxid)
            #print(fluxdict_in[rxid])
            #print('coef:',coef)
            #print(rxndict)
            totflux += fluxdict_in[rxid] * int(coef)
        fluxdict_out[expid] = totflux
    return fluxdict_out

def computesplits(fluxdict,splitsmap):
    pass
