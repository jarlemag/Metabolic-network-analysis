#gx-FBA.p
import cobra
from cobra.io.sbml import create_cobra_model_from_sbml_file
import numpy as np
from collections import defaultdict
from math import log
from copy import copy

def gxFBA(cobramodel_1,cobramodel_2,gene_expressions,GPRlist,maxflux = 500, exprType = 2):

    if exprType == 1:
        logvals = {gene:log(gene_expressions[gene],2) for gene in gene_expressions}
        ratios = gene_expressions
    else:
        logvals = gene_expressions
    

    FVAcondition = 2 #Specify which condition FVA should be performed for.
    #Local variables:
    eps2 = 0.00001
    
    cobramodel_1.optimize(solver='gurobi')

    wildtype = cobramodel_1.solution.x_dict
    wildtype_averageflux = np.mean([fluxvalue for fluxvalue in wildtype.values() if fluxvalue > eps2]) #average wildtype flux

    FVA_result = cobra.flux_analysis.variability.flux_variability_analysis(cobramodel_1,fraction_of_optimum = 0)

    
    wt_avg = {key:np.mean([FVA_result[key]['minimum'],FVA_result[key]['maximum']])  for key in FVA_result}

    '''
    T = defaultdict(list)
    for gene_id in expressionvalues: #For every gene expression value
        if gene_id in GPRlist: #If the gene identifier exists in the gene/reaction-connecting dictionary
            for reaction_id in GPRlist[gene_id]: #For every reaction associated with the gene
                print gene_id
                print reaction_id
                T[reaction_id].append(gene_id) #Add that gene to the reaction's list of available gene expression values
            #T[reaction_id].append([expressionentry for reaction_id in GPRlist[expressionentry] ])
    '''
    genemap = defaultdict(list)
    for gene in genes_expressions:
        if gene in GPRlist:
            for rx in GPRlist[gene]:
                genemap[rx].append(gene)

    T = []
    #Make a list of reactions for which data  T:
    for reaction in cobramodel_2.reactions:
        #Skip reactions for which gene expression data is unavailable
        if reaction.id not in genemap:
            continue
        #Do not include the reaction if the average flux is too high:
        if wt_avg[reaction.id] >= maxflux:
            continue
        #Do not include the reaction if it cannot carry any flux (as per FVA):
        if abs(wt_avg[reaction.id]) < eps2:
            continue
        #Do not include the reaction if it is not irreversible under the given constraints:
        bounds = FVA_result[reaction.id].values()
        irreversible = not min(bounds) < 0 < max(bounds)
        if not irrerversible:
            continue
        #Do not include the reaction if gene expression data is inconsistent (both up- and down-regulation):
        exprvalues = [logvals[gene] for gene in genemap[reaction.id] ]
        signs = [cmp(logvalue,0) for logvalue in exprvalues]
        if not (signs[1:] == lst[:-1]): #If all the signs are not the same value. Alternative: if not (signs.count(signs[0]) == signs(x))
            continue
        #If none of the above checks failed, add the reaction to the list:
        T.append(reaction)

    
    for reaction in T:
        logvalues = [logvals[gene] for gene in genemap[reaction.id] ]
        sign = cmp(logvalues[0],0)

        '''
        if wildtype[reaction.id] == 0: #If the wildtype flux is zero
            #If the reaction is up-regulated
            if sign > 0:
                Ci_mRNA = 2**max(logvalues) 
                reaction.upper_bound = Ci_mRNA * wildtype_averageflux
            else:
                Ci_mRNA = 2**min(logvalues)
                #Don't do change the lower reaction limit if the reaction is down-regulated.
        else:
            if sign > 0:
                 Ci_mRNA = 2**max(logvalues)
                 reaction.upper_bound = Ci_mRNA * wildtype[reaction.id]
            elif sign <0:
                Ci_mRNA = 2**min(logvalues)
                reaction.lower_bound = Ci_mRNA * wildtype[reaction.id]
        #Set the objective function coefficient:
                
        reaction.objective_coefficient = Ci_
        '''
        if sign > 0:
            Ci_mRNA = 2**max(logvalues)
            if wildtype[reaction.id] == 0:
                reaction.upper_bound = Ci_mRNA * wildtype_averageflux
            else:
                reaction.upper_bound = Ci_mRNA * wildtype[reaction.id]
            reaction.objective_coefficient = max(logvalues)*wildtype[reaction.id] / wt_avg[reaction.id]
        elif sign < 0:
            Ci_mRNA = 2**min(logvalues)
            if wildtype[reaction.id] == 0:
                pass
            else:
                reaction.lower_bound = Ci_mRNA * wildtype[reaction.id]
            reaction.objective_coefficient = min(logvalues)*wildtype[reaction.id] / wt_avg[reaction.id]
    return cobramodel_2.solution
##
##'''
##1) Generate the wild-type flux distribution v_i_wt for the starting condition using an Interior Point algorithm.
##'''
##
##
###'''
###2) Perform FBA with nutritional constraints per condition 2 and minimum biomass flux = 0 to calculat the upper and lower fluxes each
###    reaction can carry.
###    Calculate mean possible flux value for each reaction, and the flux average for all active reactions.
##
##    #cobramodel2.optimize(solver='gurobi')
##
##
##'''
##3) Identify the set of reactions (called T) for which an mRNA expression value can be associated.
##    Eliminate reactions that carry "unreasonably high" flux values.
##    For reactions dependent on several genes (several mRNA values), use the maximal up/down-regulation value.
##    If expression values are inconsistent (both up- and down-regulation of different genes associated with a single reaction),
##    exclude the reaction.
##
##4)  For each reaction in T; define a new upper bound UB = Ci_mRNA * v_i_wt if up-regulated or a new lower bound LB =
##    Ci_mRNA * v_i_wt if downregulated, where Ci_mRNA is the mRNA expression ratio.
##    Exclude reversible reactions.
##
##5) Construct the new objective function: Z = sum over reactions in T:  log(Ci_mRNA) *v_i / (v_i average)
##
##    For reactions that are inactive in the wild_type, v_i is set equal to the average of all active reactions (v_all average).
##    The upper bound is set equal to Ci_mRNA * (v_all average), and the lower bound equal to 0.
##
##'''
##    


if __name__ == "__main__":

    model = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')
    model.optimize(solver='gurobi')

    #fva = gxFBA(model,[],[])
