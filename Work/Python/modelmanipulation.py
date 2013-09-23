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
