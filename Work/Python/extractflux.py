#extractflux.py
'''
Extracts/calculates chosen experimental fluxes from a flux vector.
#fluxes = extractflux(fluxvector,options)
'''
import scipy.io

def extractflux(fluxvector,options):

    
    reactionmaps = scipy.io.loadmat('reactionmaps.mat')

    if ('extractmode' in options) == False:
        options['extractmode'] = 'default'
    if ('usealtmap' in options) == False:
        options['usealtmap'] = 0

    if options['extractmode'] == 'default'
        options['extractmode'] = 3
    elif options['extractmode'] == 'hardcoded'
        options['extractmode'] =1
    elif options['extractmode'] == 'simplemap'
        options['extractmode'] = 2
    elif options['extractmode'] == 'tokens'
        options['extractmode'] =3

    if ('model_id' in options) == False:
        options['model_id'] = 1
        print ' WARNING: extractflux.py: Model not specified! Model ID set to 1 by default.'

    if ('verbflag' in options) == False:
        options['verbflag'] =1
    if ('debugmode' in options) == False:
        options['debugmode'] = 1

    debugmode = options['debugmode']
    verbflag = options['verbflag']
    model_id = options['model_id']

    if debugmode >0:
        print 'extractflux.py: Debugmode:',debugmode
        print 'extractflux.py: extractmode',options['extractmode']

        
