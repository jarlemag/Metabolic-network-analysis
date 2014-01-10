#modelmanipulation.py
'''
Functions for manipulating COBRA models.
'''



def couplereactions(rx1,rx2,factor):
'''
Couple reactions together so that they always occur in a specific proportion to each other.
'''


def ratiolimit(model,reactionID,fraction):
    limit = -model.reactions.get_by_id(reactionID).lower_bound*fraction
    for reaction in model.reactions:
        if 'EX_' in reaction.id == False:



def createTokens(model,reactionmap):
    '''
    Compares a CobraPy model and a reaction map. If several model reactions are
    associated with the same experimental reaction, a token reaction is created
    whose flux will equal a linear combination of those reactions (the token flux). This token
    reaction can be used as a decision variable during quadratic optimization attempting
    to match the token flux to that of the experimental reaction. 
    Noe that limiting the token flux to zero will prevent any flux through the separate reactions.
    '''

