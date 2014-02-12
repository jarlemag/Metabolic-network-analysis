#FCF.py
#Implementation of Flux Coupling Finder
import cobra
import numpy as np
import itertools
from cobra.io.sbml import create_cobra_model_from_sbml_file

#Perform flux variability analysis


class FluxCouplings(object):
    def __init__(self,cobramodel):
        self.couplings,self.couplingsets = findcouplings(cobramodel)
        

class FluxCoupling(object):
    def __init__(self,reaction1,reaction2,couplingtype)
        self.reaction1 = reaction1
        self.reaction2 = reaction2
        self.couplingtype = couplingtype


def findcouplings(cobramodel):

    FVAres = cobra.flux_analysis.variability.flux_variability_analysis(cobramodel,fraction_of_optimum = 0)

    metabolites = cobramodel.metabolites
    nonblocked = [cobramodel.reactions.get_by_id(key) for key in FVAres if FVAres[key]['maximum'] != 0]
    for reaction in nonblocked:
        reaction.alreadyCoupled = 0

    #alreadyCoupled = np.zeros(len(metabolites))

    reactionpairs = [tuitertools.combinations(nonblocked, r = 2)
    couplings = []
    couplingsets = []
    #Iterate over couples of reactions:
    for reactionpair in reactionpairs:
        reaction1 = reactionpair[0]
        reaction2 = reactionpair[1]
        if reaction1.alreadycoupled == 1:
            continue
        #Solve for Rmin and Rmax (minimum maximum ratios of v_j/v_'j)

        #A. If R_min = 0 and Rmax is unbounded, the reactions are uncoupled.
        if R_min < tol and R_max == 'unbounded':
                     continue
        #B. If R-min = 0 and R_max = c >0 then  v:j -> v_j (The fluxes are directionally coupled)
        if R_min < tol and np.isinf(R_max) == False:
            couplings.append(FluxCoupling(reaction1.id,reaction2.id,couplingtype='directional')
        #C. if Rmin = 0 and Rmax = c > 0 then:
        if R_min < tol and R_max > tol and (np.isinf(R_max) == False):
            #If (c2-c2) > 0 then (v_j <-> v_'j)  (partially coupled)
                if (c2-c1) > 0:
                    couplings.append(FluxCoupling(reaction1.id,reaction2.id,couplingtype = 'partial')
            #if (c2-c1) = 0 then (v_j <=> v_'j) (fully coupled)
                if abs(c2-c1) < tol:
                    couplings.append(FluxCoupling(reaction1.id,reaction2.id,couplingtype='full'))
            reaction1.alradycoupled = 1
            reactions_found_in_set = False
            #If a coupled reaction set already containing one of the reactions exists,
            #Add the reactions to the set.
            for reactionset in couplingsets:
                if any(reaction.id in reactionset for reaction in reactionpair):
                    reactionset.append(reaction1.id)
                    reactionset.append(reaction2.id)
                    isinset = True
            #If not, make a new set containing the two reactions:
            if not reactions_found_in_set:
                couplingsets.append([reaction1.id,reaction2.id])
        #D. if Rmin = c > 0 and R_max is unbounded, then v'j -> vj
        if R_min > tol and np.isinf(R_max):
            couplings.append(FluxCoupling(reaction2.id,reaction1.id,couplingtype ='directional')

        return couplings,couplingsets


def maxminfluxratio(cobramodel,reaction1,reaction2):

    #Wipe the default objective funcito of the model:
    for reaction in cobramodel.reactions:
        reaction.objective_coefficient = 0
    #Set the new objective function:
    cobramodel.reactions.get_by_id(reaction1).objective_coefficient = 1

    t = 2 #Any number above 0?
    #Multiply all other constraints by t:
    for reaction in cobramodel.reactions:
        reaction.upper_bound = reaction.upper_bound *t
        reaction.lower_bound = reaction.lower_bound *t

    #Multiply al objective coefficients by t:
    for reaction in cobramodel.reactions:
        reaction._metabolites = [rx._metabolites[key] * t for key in rx._metabolites]

    #Set the constraint for the denominator reaction:
    cobramodel.reactions.get_by_id(reaction2).upper_bound = 1
    cobraodel.reactions.get_by_id(reaction2).lower_bound = 1

                         
    pass

def normalizeFluxesByReaction(fluxdict,reaction):
    pass

if __name__ == '__main__':

    SCHUETZR = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')
                 
    
