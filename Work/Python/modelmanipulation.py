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
    
