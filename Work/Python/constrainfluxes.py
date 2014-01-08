#constrainfluxes.py

def constrainfluxes(model,*args):
'''
Takes as input a COBRA model and constrains one or more reactions with a
defined corresponding experimental reaction to limits specified by the
experimental value and by a tolerance relative to the experimental uncertainty.
'''
    if len(args) < 1:
        options = {}
    else:
        options = args[0]
    if ('model_id'in options) == False:
        options['model_id'] =1

    if ('exp_id' in options) == False:
        options['exp_id'] = 1
    if ('tolerance' in options) == False:
        options['tolerance'] = 1
    if ('fluxconstrainset' in options) == False:
        options['fluxconstrainset'] = 0
    model_id = options['model_id']
    exp_id = options['exp_id']
    tolerance = options['tolerance']
    fluxconstrainset = options['fluxconstrainset']

    #matlab code: load('reactionmaps.mat')
    #matlab code: load('expdata')

    if ('description' in model) == True:
        modelname = model['description']
    else:
        modelname = '?'

    print 'constrainfluxes.py: Usng reaction map for model with model ID ',model_id
    if model_id == 1:
        #matlab code: reactionmap = reactionmaps.Fmap2
        pass
    elif model_id == 2:
        #matlab code: reactionmap = reactionmaps.Cmap2
        pass
    elif model_id == 3:
        #matlab code: reactionmap = reactionmaps.Gmap2
        pass
    print 'constrainfluxes.py: Using experimental fluxes from experiment',exp_id
    if exp_id == 1:
        #matlab code: expflux = expdata.perrenoud.abs.batch.aerobe.fluxvalues
        pass
    else:
        raise ValueError('Data not available.')

    if fluxconstrainset == 0:
        print 'Constraining fluxes. Tolerance factor:',tolerance
    '''
    for i in range(shape(reactionmap),1):
        #matlab/python mashup:if np.count_nonzero(reactionmap(i,:)) < 3
        pass
    '''      
