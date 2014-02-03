#modelmanipulation.py
'''
Functions for manipulating COBRA models.
'''



def couplereactions(rx1,rx2,factor):
    pass
'''
Couple reactions together so that they always occur in a specific proportion to each other.
'''


def ratiolimit(model,reactionID,fraction):
    limit = -model.reactions.get_by_id(reactionID).lower_bound*fraction
    for reaction in model.reactions:
        if 'EX_' in reaction.id == False:
            pass
    return


def maximizeMetaboliteProduction(model,metabolite_id):
    metabolite = model.metabolites.get_by_id(metabolite_id)
    metabolitereactions = metabolite._reaction
    for reaction in model.reactions:
            reaction.objective_coefficient = 0
    for reaction in metabolitereactions:
        reaction.objective_coefficient = reaction.get_coefficient(metabolite)
    return model


def objective_coefficient_list(model):
    for reaction in model.reactions:
        print reaction.id, reaction.objective_coefficient

def getObjectiveVector(cobramodel):
    C = [reaction.objective_coefficient for reaction in cobramodel.reactions]
    return C

def objectiveValue(x,C):
    return  np.dot(np.array(C)*np.array(x))

def getNonzeroIndexes(C):
    return [i for i, e in enumerate(C) if e != 0]

def fluxlimitlist(model):
    for reaction in model.reactions:
        print reaction.id,reaction.upper_bound,reaction.lower_bound

def getUpperBounds(model):
    return [reaction.upper_bound for reaction in model.reactions]

def getLowerBounds(model):
    return [reaction.lower_bound for reaction in model.reactions]


'''
upperrbounds = []
for upperbound in ub:
    def upperlimit(x, upperbound = upperbound): #Beware of closures/wrong scoping
        return upperbound - x
    upperbounds.append(upperlimit)


lowerbounds = []
for lowerbound in lb:
    def lowerlimit(x, lowerbound = lowerbound):
        return x - lowerbound
    lowerbounds.append(lowerlimit)
'''

def constrainFunctionsUB(ub):
    upperboundfuncs = []
    for index, upperbound in enumerate(ub):
        def upperlimit(x, upperbound = upperbound):
            return upperbound - x[index]
        upperboundfuncs.append(upperlimit)
    return upperboundfuncs


def constrainFunctionsLB(lb):
    lowerboundfuncs = []
    for index, ilowerbound in enumerate(lb):
        def lowerlimit(x, lowerbound = lowerbound):
            return x[index] - lowerbound
        lowerboundfuncs.append(lowerlimit)
    return lowerboundfuncs


def constrainFunctions(ub,lb):
    ub_funcs = []
    lb_funcs = []
    for index, upperbound in enumerate(ub):
        def upperlimit(x, index = index, upperbound = upperbound):
            return upperbound - x[index]
        ub_funcs.append(upperlimit)
    for index, lowerbound in enumerate(lb):
        def lowerlimit(x, index = index, lowerbound = lowerbound):
            return x[index] - lowerbound
        lb_funcs.append(lowerlimit)
    return ub_funcs, lb_funcs

def ssConsFuncspos(S):
    ss_funcs = []
    for row in S:
        def steadystatconstr(x, row = row ):
            return sum(np.multiply(row,x))
        ss_funcs.append(steadystatconstr)
    return ss_funcs

def ssConsFuncsposneg(S):
    ss_funcs = []
    for row in S:
        def steadystatconstr(x, row = row ):
            return -sum(np.multiply(row,x))
        ss_funcs.append(steadystatconstr)
    return ss_funcs


def createTokens(model,reactionmap):
    from cobra import Reaction
    from cobra import Metabolite
    '''
    Compares a CobraPy model and a reaction map. If several model reactions are
    associated with the same experimental reaction, a token reaction is created
    whose flux will equal a linear combination of those reactions (the token flux). This token
    reaction can be used as a decision variable during quadratic optimization attempting
    to match the token flux to that of the experimental reaction. 
    Noe that limiting the token flux to zero will prevent any flux through the separate reactions.
    '''
    for linkdict in reactionmap:
        if len(linkdict['modrxns']) > 1:
            experimental_r_id = linkdict['exprxns'][0]['rxid']
            token_id = 't_{exprxnid}'.format(exprxnid = experimental_r_id)
            new_token_reaction = Reaction(token_id)
            new_token_reaction.name = 'Token exchange reaction for experimental reaction {exprxnid}'.format(exprxnid = experimental_r_id)
            new_token_reaction.lower_bound = -1000
            new_token_reaction.upper_bound = 1000
            new_token_metabolite = Metabolite(token_id,name = 'Token metabolite for experimental reaction {exprxnid}'.format(exprxnid = experimental_r_id),compartment='c')
            new_token_reaction.add_metabolites({new_token_metabolite: -1.0})
            model.add_reactions(new_token_reaction)
            for dictionary in linkdict['modrxns']:
                modelreaction = model.reactions.get_by_id(dictionary['rxid'])
                modelreaction.add_metabolites({new_token_metabolite: float(dictionary['coef'])})
    return model


if __name__ == "__main__":

    from cobra.io.sbml import create_cobra_model_from_sbml_file
    from cobra.io.sbml import write_cobra_model_to_sbml_file
    import loadData as load

    notokenmodel = create_cobra_model_from_sbml_file('../SBML/SCHUETZR_notokens.xml')
    notokenmap = load.ReactionMapfromXML('reactionmaps_notokens.xml','Perrenoud','SCHUETZR')
                                            
    tokenmodel = createTokens(notokenmodel,notokenmap)
    sbml_out_file = 'SCHUETZR_withtokens.xml'
    sbml_level = 2
    sbml_version = 1
    write_cobra_model_to_sbml_file(tokenmodel, sbml_out_file, sbml_level,sbml_version, print_time=False)
    
