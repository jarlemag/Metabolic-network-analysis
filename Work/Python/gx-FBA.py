#gx-FBA.p
import cobra
from cobra.io.sbml import read_sbml_model
from cobra.io.sbml import write_sbml_model
import numpy as np
from collections import defaultdict
from math import log
from copy import deepcopy
import QPmindist
import random

def hold():
    x = raw_input('Press any key to continue.')
    print('\n')


def vectorToText(fluxvector,filename,decimal = ',', append = False, newline = True, linestart = None):
    if append:
        target = open(filename, 'a')
    else:
        target = open(filename, 'w')
    if linestart is not None:
        target.write(linestart)
        target.write('\t')
    for flux in fluxvector:
        target.write(str(flux).replace('.',decimal))
        target.write('\t')
    if newline:
        target.write('\n')
    target.close()
    return


def FAMEtoTXT(filein,fileout, decimal = ','):
    with open(filein, 'r') as sourcefile:
        target = open(fileout, 'w')
        for line in sourcefile:
            flux = line.split()[1].replace('.',decimal)
            #print 'line:',line
            #print 'flux:',flux
            target.write(flux)
            target.write('\t')
        target.close()
    return
        
    


def metabolitefluxsum(metabolite_id,cobramodel, abs = False):
    if abs:
        return sum([abs(cobramodel.solution.x_dict[reaction.id]*reaction.get_coefficient(metabolite_id)) for reaction in cobramodel.metabolites.get_by_id(metabolite_id).get_reaction()])
    else:
        return sum([cobramodel.solution.x_dict[reaction.id]*reaction.get_coefficient(metabolite_id) for reaction in cobramodel.metabolites.get_by_id(metabolite_id).get_reaction()])
    
def metabolitefluxes(metabolite_id,cobramodel):
    for reaction in cobramodel.metabolites.get_by_id(metabolite_id).get_reaction():
        print(reaction.id,reaction.get_coefficient(metabolite_id)*cobramodel.solution.x_dict[reaction.id])

def gxFBA(cobramodel_1,cobramodel_2,gene_expressions,GPRlist,maxflux = 500, exprType = 2, verbose = True,FVAcondition = 2, wait = False,treshold = 0.5,objectiveweights = 'log',dumpmodel = False):
    if verbose:
        print('Performing gxFBA. Verbose = True.')
        print('Provided gene expression ratios:')
        print(gene_expressions)
        if exprType == 1:
            print('Input gene expression ratios handled as absolute ratios.')
        else:
            print('Input gene expression ratios handled as log ratios.')

    if exprType == 1:
        logvals = {gene:log(gene_expressions[gene],2) for gene in gene_expressions}
        ratios = gene_expressions
        if verbose:
            print('Calculated log ratios:')
            print(logvals)
    else:
        logvals = gene_expressions
        ratios = {gene:2**gene_expressions[gene] for gene in gene_expressions}
    
    #Local variables:
    eps2 = 0.00001

    ##1) Generate the wild-type flux distribution v_i_wt for the starting condition:
    cobramodel_1.optimize(solver='gurobi')
    wildtype = cobramodel_1.solution.x_dict
    if verbose:
        print('generated wild-type flux distribution.')
        #print cobramodel_1.solution.x_dict
        for key,value in sorted(cobramodel_1.solution.x_dict.items()):
            print(key,value)
    if wait:
        hold()
    ###2) Perform FVA with nutritional constraints for condition 1 or 2 (Default: condition 2) and minimum biomass flux = 0 to calculate the upper and lower fluxes each
    ###    reaction can carry.

    if FVAcondition == 1:
        FVA_result = cobra.flux_analysis.variability.flux_variability_analysis(cobramodel_1,fraction_of_optimum = 0)
    else:
        FVA_result = cobra.flux_analysis.variability.flux_variability_analysis(cobramodel_2,fraction_of_optimum = 0)
    if verbose:
        print('Performed FVA for condition',FVAcondition)
        for key,value in sorted(FVA_result.items()):
            print(key,value)
        #print FVA_result

    #Reset the objective function for condition 2:
    for reaction in cobramodel_2.reactions:
        reaction.objective_coefficient = 0

    ###    Calculate mean possible flux value for each reaction under condition 2, and the flux average for all active reactions under condition 1
    wt_avg = {key:np.mean([FVA_result[key]['minimum'],FVA_result[key]['maximum']])  for key in FVA_result}
    active_fluxes = [fluxvalue for fluxvalue in wildtype.values() if fluxvalue > eps2]
    wildtype_averageflux = np.mean(active_fluxes) #average wildtype flux

    if verbose:
        print('FVA average fluxes:')
        #print wt_avg
        for key,value in wt_avg.items():
            print(key,value)
        print('Total number of reactions in wildtype model:',len(wildtype.values()))
        print('Number of active reactions in wildtype FBA solution:',len(active_fluxes))
        print('Average wildtype flux (Average of flux value for all fluxes in wildtype FBA solution):',wildtype_averageflux)

    genemap = defaultdict(list)
    for gene in gene_expressions:
        if gene in GPRlist:
            for rx in GPRlist[gene]:
                genemap[rx].append(gene)

    if verbose:
        print('Gene map:')
        print(genemap)


    if wait:
        hold()
    T = []
    ##3) Identify the set of reactions (called T) for which an mRNA expression value can be associated.
    if verbose:
        print('Identifying reactions eligible for gx-regulation.')
    for reaction in cobramodel_2.reactions:
        #Skip reactions for which gene expression data is unavailable
        if reaction.id not in genemap:
            if verbose:
                print(reaction.id,'not in genemap. Skipping reaction.')
            continue
        #Skip reactions for which the maximal up/down-regulation is below the set treshold:
        expressionratios = [ratios[gene] for gene in genemap[reaction.id] ]
        if max([abs(ratio-1) for ratio in expressionratios]) < treshold:
            if verbose:
                print(reaction.id,': Level of regulation below treshold. Skipping reaction.')
            continue
        #Do not include the reaction if the average flux is too high:
        if wt_avg[reaction.id] >= maxflux:
            if verbose:
                print(reaction.id,': FVA average flux too high. Skipping reaction.')
            continue
        #Do not include the reaction if it cannot carry any flux (as per FVA):
        if abs(wt_avg[reaction.id]) < eps2:
            if verbose:
                print(reaction.id,': Reaction always inactive under specified conditions. Skipping reaction.')
            continue
        #Do not include the reaction if it is not irreversible under the given constraints:
        bounds = FVA_result[reaction.id].values()
        irreversible = not min(bounds) < 0 < max(bounds)
        if not irreversible:
            if verbose:
                print(reaction.id,': Reaction may carry flux in both directions under specified conditions. Skipping reaction.')
            continue
        #Do not include the reaction if gene expression data is inconsistent (both up- and down-regulation):
        exprvalues = [logvals[gene] for gene in genemap[reaction.id] ]
        signs = [cmp(logvalue,0) for logvalue in exprvalues]
        if not (signs[1:] == signs[:-1]): #If all the signs are not the same value. Alternative: if not (signs.count(signs[0]) == signs(x))
            if verbose:
                print('Gene expression information for reaction indicates mixed regulation (up/down). Skipping reaction.')
            continue
        #If none of the above checks failed, add the reaction to the list:
        if verbose:
            print('Adding reaction',reaction.id,'to T.')
        T.append(reaction)


    ##4)  For each reaction in T; define a new upper bound UB = Ci_mRNA * v_i_wt if up-regulated or a new lower bound LB =
    ##    Ci_mRNA * v_i_wt if downregulated, where Ci_mRNA is the mRNA expression ratio.
    if verbose:
        print('Calculating new bounds and objective function coefficients:')


    def reactionInfo(reaction_id):
        print('Ci_mRNA:',Ci_mRNA)
        print('Wild-type flux:',wildtype[reaction.id])
        print('Original reaction bounds:')
        print('Lower:',cobramodel_1.reactions.get_by_id(reaction.id).lower_bound,'Upper:',cobramodel_1.reactions.get_by_id(reaction.id).upper_bound)
        print('New reaction bounds:')
        print('Lower:',reaction.lower_bound,'Upper:',reaction.upper_bound)
        print('Objective coefficient:',reaction.objective_coefficient)
    for reaction in T:
        if wait:
            hold()
        logvalues = [logvals[gene] for gene in genemap[reaction.id] ]
        expressionratios = [ratios[gene] for gene in genemap[reaction.id] ] 
        sign = cmp(logvalues[0],0)
        if verbose:
            print('Processing reaction:',reaction.id)
            print('Logvalues:',logvalues)
            print('Ratios:',expressionratios)
            
        if sign > 0:
            Ci_mRNA = 2**max(logvalues) ##    For reactions dependent on several genes (several mRNA values), use the maximal up/down-regulation value.
            if wildtype[reaction.id] == 0:
                reaction.upper_bound = Ci_mRNA * wildtype_averageflux
            else:
                reaction.upper_bound = Ci_mRNA * wildtype[reaction.id]
            #Set the objective function coefficient: Ci = log(Ci_mRNA) / (v_i average):
            if objectiveweights == 'log':
                reaction.objective_coefficient = max(logvalues) / wt_avg[reaction.id]
            elif objectiveweights =='logabs':
                reaction.objective_coefficient = max(logvalues) / wt_avg[reaction.id]
            elif objectiveweights == 'ratios':
                reaction.objective_coefficient = Ci_mRNA / wt_avg[reaction.id]
            elif objectiveweights == 'pureratios':
                reaction.objective_coefficient = Ci_mRNA
            
        elif sign < 0:
            Ci_mRNA = 2**min(logvalues) ##    For reactions dependent on several genes (several mRNA values), use the maximal up/down-regulation value.
            if wildtype[reaction.id] == 0:
                pass
            else:
                reaction.lower_bound = Ci_mRNA * wildtype[reaction.id]
            #Set the new objective function coefficient: Ci = log(Ci_mRNA) / (v_i average):
            if objectiveweights == 'log':
                reaction.objective_coefficient = min(logvalues) / wt_avg[reaction.id]
            elif objectiveweights =='logabs':
                reaction.objective_coefficient = abs(min(logvalues)) / wt_avg[reaction.id]
            elif objectiveweights =='ratios':
                reaction.objective_coefficient = Ci_mRNA / wt_avg[reaction.id]
            elif objectiveweights == 'pureratios':
                reaction.objective_coefficient = Ci_mRNA
        if verbose:
            reactionInfo(reaction.id)

    cobramodel_2.reactions.get_by_id('R5_akgdh').lower_bound = 0.65
    cobramodel_2.reactions.get_by_id('R5_akgdh').upper_bound = 1.3
    cobramodel_2.reactions.get_by_id('R2_lowergly').lower_bound = 0.8
    cobramodel_2.reactions.get_by_id('R2_lowergly').upper_bound = 3.42857

    cobramodel_2.reactions.get_by_id('R2_lowergly').objective_coefficient = 1.5
    cobramodel_2.reactions.get_by_id('R5_akgdh').objective_coefficient = 0.5
    cobramodel_2.reactions.get_by_id('R1_uppergly').objective_coefficient = 1.5
    cobramodel_2.reactions.get_by_id('R3_citsyn').objective_coefficient = 0.5
    if verbose:
        print('New objective function:')
        #print {reaction.id:reaction.objective_coefficient for reaction in T}
        print({reaction.id:reaction.objective_coefficient for reaction in cobramodel_2.reactions})
        print('New lower bounds:',[reaction.lower_bound for reaction in cobramodel_2.reactions])
        print('New upper bounds:',[reaction.upper_bound for reaction in cobramodel_2.reactions])
        print('Performing gx-FBA.')
        

    #Solve the gx-FBA problem:
    cobramodel_2.optimize()
    if verbose:
        print('Flux solution:')
        print(cobramodel_2.solution.fluxes)
        print('Objective value:',cobramodel_2.solution)

    #Using Gurobi manually:
    #Commenting out. Unable to use gurobi at the moment due to lack of license.

    #gurobisolution = QPmindist.gurobiFBA(cobramodel_2)
    #gurobisolutiondict = QPmindist.getgurobisolutiondict(gurobisolution)
    #print('gurobi solution:')
    #print(gurobisolutiondict)
    #print('Objective value:',gurobisolution.Objval)

    #Check that mass conservation is preserved:
    print('Total metabolite fluxes:')
    for metabolite in cobramodel_2.metabolites:
        print(metabolite.id,metabolitefluxsum(metabolite.id,cobramodel_2),metabolitefluxsum(metabolite.id,cobramodel_2, abs = True))

    #print 'ATP fluxes:',metabolitefluxes('atp_c',cobramodel_2)
    #print ''

    print('Solution:')
    print(cobramodel_2.solution.x)
    vectorToText(cobramodel_2.solution.x,'flux.txt')
    print('Flux solution dumped to flux.txt')
    #vectorToText([reaction.objective_coefficient for reaction in cobramodel_2.reactions],'objective.txt')
    vectorToText([reaction.objective_coefficient for reaction in cobramodel_2.reactions],'objective.txt', decimal = ',')
    vectorToText([reaction.lower_bound for reaction in cobramodel_1.reactions],'FBA-lb.txt')
    vectorToText([reaction.upper_bound for reaction in cobramodel_1.reactions],'FBA-ub.txt')
    vectorToText([reaction.lower_bound for reaction in cobramodel_2.reactions],'gxFBA-lb.txt')
    vectorToText([reaction.upper_bound for reaction in cobramodel_2.reactions],'gxFBA-ub.txt')

    vectorToText([reaction.objective_coefficient for reaction in cobramodel_2.reactions],'A.txt', decimal = ',', linestart = 'C:')
    #vectorToText([reaction.lower_bound for reaction in cobramodel_1.reactions],'A.txt',append = True, linestart = 'LB:')
    #vectorToText([reaction.upper_bound for reaction in cobramodel_1.reactions],'A.txt', append = True,linestart = 'UB:')
    vectorToText([reaction.lower_bound for reaction in cobramodel_2.reactions],'A.txt', append = True, linestart = 'LB:')
    vectorToText([reaction.upper_bound for reaction in cobramodel_2.reactions],'A.txt.', append = True, linestart = 'UB:')
    vectorToText(cobramodel_2.solution.x,'A.txt.', append = True, linestart = 'V:')

    if dumpmodel:
        write_sbml_model(cobramodel_2,'dumpedmodel.xml')
    
    return cobramodel_2



def sampleOptimalSolutions(cobramodel,biomassreaction = None):
    tol = 1e-3
    iterations = 100
    #Identify the optimal value of the FBA objective function:
    cobramodel.optimize()
    z = cobramodel.solution.f
    #For each model reaction, use FVA to identify the upper and lower flux bound, subject to FBA objective optimality:
    FVAres = cobra.flux_analysis.variability.flux_variability_analysis(cobramodel)
    #Identify the reaction set R for which the upper and lower bounds are unequal (the flux is allowed to vary):
    R = [reaction for reaction in FVAres.keys() if abs(FVAres[reaction]['minimum']-FVAres[reaction]['maximum']) > tol]

    itercount = 0
    Z_gxFBA = []
    BMmax = []
    BMmin = []
    while itercount < iterations:
        itercount +=1

        #Select a random reaction:
        randomreaction = random.choice(R)
        randomflux = random.uniform(FVAres[randomreaction]['minimum'],FVAres[randomreaction]['maximum'])
        cobramodel.reactions.get_by_id(randomreaction).upper_bound = randomflux
        cobramodel.reactions.get_by_id(randomreaction).lower_bound = randomflux
        cobramodel.optimize()

        gxFBAsolution = gxFBA(cobramodel) #Does this ,make sense?

        Z_gxFBA.append(gxFBAsolution.solution.f)
        if biomassreaction is not None:
            bmFVA = cobra.flux_analysis.variability.flux_variability_analysis(cobramodel,the_reactions = [biomassreaction])
            BMmax.append(bmFVA[biomassreaction]['maximum'])
            BMmin.append(bmFVA[biomassreaction]['minimum'])
   
    return R,FVAres

if __name__ == "__main__":

    model = read_sbml_model('../SBML/SCHUETZR.xml')
    model.optimize()

    #fva = gxFBA(model,[],[])

    example_a = read_sbml_model('../SBML/gx-fba.xml')
##    gx-FBA example model 1:
##    Parameters:
##    Objective: Biomass
##    ATP maintenance: 12
##    Max glucose import (glc_e -> glc_c): 10
##    Max puruvate export (pyr_c -> pyr_e): 5
    example_a.optimize()
    #print 'Example A FBA result:'
    #print example_a.solution
    #print example_a.solution.x_dict

    example_b = read_sbml_model('../SBML/gxfba_example.xml')


   # example_c = create_cobra_model_from_sbml_file('../SBML/gx-fba_corrected.xml') #The model file has disappeared. WTH?
    #example_c.optimize(solver='gurobi')
##    gx-FBA example model 2:
##    Parameters:
##    Objective: Biomass (R8_Biomass)
##    ATP maintenancee: 12 (R9_atpm)
##    Max glucose exchange: -10 (R_glcb_ex)
##    Max pyruvate exchange: 5 (R_pyr_ex)
    example_b.optimize()
    #print 'Example B FBA result:'
    #print example_b.solution
    #print example_b.solution.x_dict

    gene_expressions = {'V1':1.5,'V2':1.5,'V3':0.5,'V5':0.5}
    #GPRlist = {'R1_uppergly':'V1','R2_lowergly':'V2','R3_citsyn':'V3','R5_akgdh':'V5'}
    #GPRlist_b = {y:[x] for x,y in GPRlist.iteritems()}
    GPRlist = {'V1':['R1_uppergly'],'V2':['R2_lowergly'],'V3':['R3_citsyn'],'V5':['R5_akgdh']}

    cobramodel_1 = example_b
    cobramodel_2 = deepcopy(cobramodel_1)

    #gx = gxFBA(cobramodel_1,cobramodel_2,gene_expressions,GPRlist,exprType = 1, wait = True)

    #gx2 = gxFBA(cobramodel_1,cobramodel_2,gene_expressions,GPRlist,exprType = 1, wait = True, objectiveweights = 'logabs')

    #gx3 = gxFBA(cobramodel_1,cobramodel_2,gene_expressions,GPRlist,exprType = 1, wait = True, objectiveweights = 'ratios')

    #cobramodel_1 = example_b
    cobramodel_2 = deepcopy(cobramodel_1)
    
    #gx4 = gxFBA(cobramodel_1,cobramodel_2,gene_expressions,GPRlist,exprType = 1, wait = True, objectiveweights = 'pureratios')

    #gx5 = gxFBA(cobramodel_1,cobramodel_2,gene_expressions,GPRlist,exprType = 1, wait = True, dumpmodel = True)
    #Confirm FVA result for upper_gly by optimizing with upper_gly as objective:
    #example_b.reactions.get_by_id('R8_Biomass').objective_coefficient = 0
    #example_b.reactions.get_by_id('R1_uppergly').objective_coefficient = 1
    #example_b.reactions.optimize()
    #print example_b.solution
    atp_reactions = ['R1_uppergly','R2_lowergly','R6_fumar','R8_biomass','R9_atpm']

    #metabolitefluxes('ATP_c',model)

    gprlist2 = {'V2':['R2_lowergly'],'V5':['R5_akgdh']}

    #gx6 = gxFBA(cobramodel_1,cobramodel_2,gene_expressions,GPRlist,exprType = 1, wait = True, dumpmodel = True,objectiveweights = 'pureratios')
    
   
    TEST = example_b
    
    TEST.reactions.get_by_id('R2_lowergly').lower_bound = 0.8
    TEST.reactions.get_by_id('R2_lowergly').upper_bound = 3.42857


    TEST.reactions.get_by_id('R5_akgdh').lower_bound = 0.65
    TEST.reactions.get_by_id('R5_akgdh').upper_bound = 1.3

    TEST.reactions.get_by_id('R8_Biomass').objective_coefficient = 0
    
    TEST.reactions.get_by_id('R2_lowergly').objective_coefficient = 1.5
    TEST.reactions.get_by_id('R5_akgdh').objective_coefficient = 0.5
    TEST.reactions.get_by_id('R1_uppergly').objective_coefficient = 1.5
    TEST.reactions.get_by_id('R3_citsyn').objective_coefficient = 0.5

    #TEST.reactions.get_by_id('R8_Biomass').objective_coefficient = -1
    TEST.optimize(solver='gurobi')
    
    print('solution:')
    
    print('R1:',TEST.solution.x_dict['R1_uppergly'])
    print('R2:',TEST.solution.x_dict['R2_lowergly'])
    print('R3:',TEST.solution.x_dict['R3_citsyn'])
    print('R4:',TEST.solution.x_dict['R4_acon'])
    print('R5:',TEST.solution.x_dict['R5_akgdh'])
    print('R6:',TEST.solution.x_dict['R6_fumar'])
    print('R7:',TEST.solution.x_dict['R7_oxysh'])
    print('R8:',TEST.solution.x_dict['R8_Biomass'])
    print('R9:',TEST.solution.x_dict['R9_atpm'])
    print('R13:',TEST.solution.x_dict['R13_pyrexport'])
    print('R15:',TEST.solution.x_dict['R15_co2export'])

    print('Objective:',[reaction.objective_coefficient for reaction in TEST.reactions])
    print('Lower bounds:',[reaction.lower_bound for reaction in TEST.reactions])
    print('Upper bounds:',[reaction.upper_bound for reaction in TEST.reactions])

    vectorToText([reaction.objective_coefficient for reaction in TEST.reactions],'B.txt', decimal = ',', linestart = 'C:')
    vectorToText([reaction.lower_bound for reaction in TEST.reactions],'B.txt', append = True, linestart = 'LB:')
    vectorToText([reaction.upper_bound for reaction in TEST.reactions],'B.txt.', append = True, linestart = 'UB:')
    vectorToText(TEST.solution.x,'B.txt.', append = True, linestart = 'V:')
    

    SCHUETZR = read_sbml_model('../SBML/SCHUETZR.xml')
    R,FVAres = sampleOptimalSolutions(SCHUETZR)    
    
    
    
