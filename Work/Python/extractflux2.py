#extractflux2.py
#Rewrite from scratch of extractflux

import numpy as np

def extractflux(fluxvector,reactionmap):
    
    fluxes = []

    for row in reactionmap:
        rxns = row[1:]
        totflux = 0
        for rxn in rxns:
            totflux += np.sign(rxn)*fluxvector[rxn-1]
        fluxes.append(totflux)
    return fluxes


