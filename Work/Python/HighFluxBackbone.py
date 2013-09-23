#HighFluxBackbone.py
# -*- coding: cp1252 -*-
import numpy as np
import cobra as cb
def HFB(model,flux_vector)
''' Based on the HighFluxBackbone.m MATLAB script by E. Almaas. Converted to Python by J. Pahr.
    Calculates the High-Flux Backbone.
    Citation: E. Almaas, B. Kovas, T Vicsek, Z.N. Oltvai and A.-L. Barabási.
    "Global organization of metabolic fluxes in the bacterium Escherichia coli."
    Nature 427, 839 (2004).

    Variables:
    model:  contains the CobraPy-format constraint-based model
    flux_vector contains the CobraPy-format flux vector from which the HFB will
    be extracted.
    numb: Returns the number of reactions participating in the HFB.
    res: Is a vector of all reactions with 0 (no) or 1 (yes) denoting
    if that reaction is part of the HFB.
    names: Is a vector of reaction names of the HFB participant reactions.

    EXAMPLE use:
    model.optimize()
    [numb,res,names] = HFB(model,model.solution.x)
'''

#Find all largest consuming and producing reactions:

    model = cb.core.ArrayBasedModel(model)
    mets = len(model.metabolites)
    rxns = len(model.reactions)


    prod_can = zeros(rxns,1)
    pp = 0
    cons_can = np.zeros(rxns)
    pc = 0

    for i in range(len(mets)): #For every metabolite
        j =  np.nonzero(model.S.getrow(i))[1] # Find the indices of all non-zero elements in the stoichoimetric matrix row corresponding to the current metabolite.
        print 'j:',j
        n = len(j)
        print 'n:',n
        bigp = 1e-7
        bign = 1e-7
        indp = np.zeros([1,200])
        indnx = np.zeros([1,200])
        fp = 0
        fn = 0
        if n > 0:
            #Loop over all the reactions metabolite i participates in.
            for k in range(n):
                value = model.S[i,j[k]]*flux_vector[j[k]]
                if (value >= bigp):
                    #Given that the metabolite is a product
                    if value == bigp:
                        fp +=1
                        indp[fp] = j[k]
                    else:
                        fp = 1
                        bigp = value
                        indp[fp] = j[k]
                elif value < bign:
                    #Given that the metabolite is a substrate
                    if value == bign:
                        fn = fn +1
                        indn{fn] = j[k]
                    else:
                        fn = 1
                        bign = value
                        indn[fn] = j[k]
            #Take care of any producing-reaction candidates identified:
            if fp > 0:
                for l in range(fp):
                    if #matlab code: ismember(indp(l),prod_cand)):
                        pass
                    else:
                        pp +=1
                        prop_cand[pp] = indp[l]
            #Take care of any consuming-reaction candidates identified:
            if fn > 0:
                for l in range(fn):
                    if #matlab code: ismember(indn(l),cons_cand))
                        pass
                    else:
                        pc +=1
                        cons_cand[pc] = indn[l]
    res = np.zeros([rxns,1])
    '''
    Consolidate the two candidate lists:For a reaction to be member of the HFB,
    it must be a maximal production and consumption reaction, not necessarily
    for the same metabolite.
    '''

    for i in range(pp):
        if #matlab code: ismember(prod_can(i),cons_cand)
            res[prod_cand[i] = 1
    #matlab code: list = find(res)
    #matlab code: numb = size(list,1)
    #names = #get the names of the reactions.

    return numb,res,names
