import numpy as np
import cobra
import sympy
#See Famili & Palsson 2003. "Systemic metabolic reactions are obtained by singular value decomposition (...)"



def findEigenReaction(cobramodel):
    cobramodel.to_array_based_model()
    u,s, vh = np.linalg.svd(cobramodel.S)
    pass

def plotSingularValues():
    pass





if __name__ == '__main__':
    testmodel = cobra.core.Model('Testmodel')
    metabolite_a = cobra.core.Metabolite('A')
    metabolite_b= cobra.core.Metabolite('B')
    metabolite_c= cobra.core.Metabolite('C')

    reaction_1 = cobra.core.Reaction('v1')
    reaction_2  = cobra.core.Reaction('v2')

    reaction_1.add_metabolites({metabolite_a:-1,metabolite_b:-1,metabolite_c:1})
    reaction_2.add_metabolites({metabolite_a:-1,metabolite_b:1})

    testmodel.add_reactions([reaction_1,reaction_2])
    
