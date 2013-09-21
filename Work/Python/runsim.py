#runsim.py
import numpy as np

#def runsim(*arg):
    
'''
%runsim(options)
%Perform a simulation to find minimum (and in some cases, maximum) distance
%between computed and experimental fluxes given selected model constraints and objective-optimality requirement:

%syntax: result = runsim(options)

%"options" is a structure with the following possible fields.
%model_id
%constraints
%objective
%exp_id
%optreq
%usesplits
%FBAonly
%verbflag
%iterations
%usegurobi
%fluxreport
%normchoice
%usescaledfluxes
%example:
%result = runsim()
'''

args = [] #Placeholder. Fix the function!
if len(args)<1:
    options = {}
    options['verbflag'] =1
    result = {}

if 'verbflag' in options == False:
    options['verbflag'] = 1 #By default, show output about the process as the program runs.

verbflag = options['verbflag']

if 'debugmode' in options == False:
    options['debugmode'] = 0

'''
if options['debugmode'] == 1: #Keyerror. Fix this.
    options['verbflag'] = 1
'''

if 'model_id' in options == False:
    options['model_id'] = 1
if verbflag == 1:
    print 'Runsim.py: ModelID not specified in input. Set to 1 (SCHUETZR) by default'

if 'objective' in options == False:
    options['objective'] = 1
    if verbflag == 1:
        print 'Objective not specified in input. Set to biomass by default.'

if 'exp_id' in options == False:
    options['exp_id'] = 1
    if verbflag == 1:
        print 'Experimental condition to compare with not specified in input. Set to aerobic batch culture by default.'

if 'optreq' in options == False:
    options['optreq'] = 1
    if verbflag == 1:
        print 'Optimality requirement not specified in input. Set to 1 by default.'

if 'usesplits' in options == False:
    options['usesplits'] = 0
    if verbflag == 1:
        print 'Use of raw fluxes or split ratios not specified in input. Set to raw fluxes by default.'

'''
if options['usesplits'] == 2: #Fix this.
    options['substractBM'] = 1
'''

if 'FBAonly' in options == False:
    options['FBAonly'] = 0

if 'iterations' in options == False:
    options['iterations'] = 0
    if verbflag == 1:
              print 'Number of random starting points for fmincon not specified. Set to zero by default (Use only FBA result).'

if 'usegurobi' in options == False:
    options['usegurobi'] = 1

if 'maxwithgurobi' in options == False:
    options['maxwithgurobi'] = 0 #%Change this if sucessful in getting maximization with Gurobi to work.

''' <<<This code doesn't do anything. Evaluate for correctness/necessity.>>>
if ('usegurobi' in options == 0) & ('usefmincon' in options == False):
    if verbflag == 1:
              print 'Solver not specified. Set to fmincon by default.'
'''

if 'usefmincon' in options == False:
    if options['model_id'] == 3:
        options['usefmincon'] = 0
    else:
        options['usefmincon'] = 1

if 'preloadmodels' in options == False:
    options['preloadmodels'] = 1
'''%Use to specify if model structure should be loaded from a .mat file, or be constructed from a SBML file using the COBRA function readCbModel.
    %By default, use pre-loaded model data. Remember to update model-data
    %if changes are made to the models!'''

'''
%Allow the posssibility of excluding some reactions from optimization
%procedures'''

if 'excludereactions' in options == False:
    options['excludereactions'] = 0 #By default, no reactions are excluded.
'''  %To exclude experimental reactions, set options.excludereactions = a
    %vector with the experimental reaction #s.
    %example: options.excludereactions = [1 12 24].'''

if 'useweights' in options == False:
    options['uncertweights'] = 0

if 'makeplots' in options == False:
    options['makeplots'] = 0

if 'usescaledfluxes' in options == False:
    options['usescaledfluxes'] = 0

if 'fluxreport' in options == False:
    options['fluxreport'] = 1

if 'computercorrcoef' in options == 0:
    options['computecorrcoef'] = 1

if 'normchoice' in options == False:
    options['normchoice'] = 'euclidean'

if 'usealloptions' in options == False:
    options['usealloptions'] = 0

if 'printfile' in options == False:
    options['writefile'] = 0

if 'usetokens' in options == False:
    options['usetokens'] = 1

if 'plotFBA' in options == False:
    options['plotFBA'] = 0

if 'extractmode' in options == False:
    options['extractmode'] = 'default'

if 'fluxconstrainset' in options == False:
    options['fluxconstrainset'] = 0

if 'reportsplits' in options == False:
    options['reportsplits'] = 0

if 'constrainglucose' in options == False:
    options['constrainglucose'] = 2

if 'constrainBM' in options == False:
    options['constrainBM'] = 1

if 'usealtmap' in options == False:
    options[usealtmap] = 0

''' #Fix this.
if options['usealloptions'] == 1:
    options['fluxreport'] = 1
    options['computecorrcoef'] = 1
    options['makeplots'] = 1
    options['usegurobi'] = 1
    options['usefmincon'] = 1
    options['iterations'] = 10
    options['verbflag'] = 1
'''

''' #Fix this.
model_id = options['model_id'] Fix this
constraints = options['constraints'] 
objective = options['objective']
exp_id = options['exp_id']
optreq = options['optreq']
usesplits = options['usesplits']
FBAonly = options['FBAonly']
usegurobi = options['usegurobi']
iterations = options['iterations']
makeplots = options['makeplots']
debugmode = options['debugmode']
usefmincon = options['usefmincon']
computecorrcoef = options['computecorrcoef']
fluxconstrainset = options['fluxconstrainset']
constrainglucose = options['constrainglucose']
constrainBM = options['constrainBM']
'''

''' #Fix this.
if options['preloadmodels'] == 1:
    if verbflag == 1:
        print 'Loading COBRA model data...'
              
    #Load preprocessed model data.
'''

''' #Fix this.
if options['excludereactions'] == 'default':
    if model_id == 1:
              options['excludereactions'] = [1,12,18,19,22,23,24,25]
    elif model_id == 2:
              options['excludereactions'] = [1,11,12,19,22,23,25,25]
    elif model_id == 3:
              options['excludereactions'] = [1,11,12,18,19,22,23,24,25]
'''

''' #Fix this:
if verbflag == 1:
    print 'Runsim.m: Excluding the folloing experimental reactions:'

    for reaction in options['excludereactions']:
              print reaction
'''
if verbflag == 1:
    print 'Loading experimental data...'

#Load experimental data


if verbflag == 1:
    print 'Setting reference experimental fluxes...'

#Load the experimental flux data

#Originally, .mat files were loaded here. Have to find out what kind of filetype to use with Python.
'''
if exp_id == 1: 
    #refexpfluxes = expdata.perrenoud.abs.batch.aerobe.fluxvalues
    if verbflag == 1:
        print 'Referenc experimental fluxes set: Perrenoud 2005 batch aerobe.'

elif exp_id == 2:
    print 'Experimental flux data not available. #Not really true. I got the data in the end..

elif exp_id == 3:
    print 'Experimental flux data not available. #Not really true. I got the data in th
'''

''' Fix this.             
if debugmode != 1:
              #Turn off warnings from COBRA regarding model read-in.
    pass
'''
if verbflag == 1:
    print 'Options used are:'
    for i in zip(options.keys(),options.values()):
        print i[0],':',i[1]

''' Fix this
if model_id == 1:
    if verbflag == 1:
              print 'Using Schuetz revised (SCHUETZR) model.'
    if options['preloadmodels'] != 1:
        model = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')
    else:
              #model = modeldata.SCHUETZR
        pass #Remove pass when above line is fixed.
elif model_id == 2:
    if verbflag == 1:
        print 'Using expanded E. coli core (ECME) model.'
    if options['preloadmodels'] != 1:
        model = create_cobra_model_from_sbml_file('../SBML/ECME.xml')    
    else:
              #model = modeldata.ECME
        pass #Remove when above line is fixed

elif model_id ==3:
    if verbflag == 1:
        print 'Using JO1366 genome-scale model.'
    if options['preloadmodels'] !=1:
        model = create_cobra_model_from_sbml_file('../SBML/iJO1366b.xml')
    else:
        #model = modeldata.iJO1366b.
        pass #Remove when above line is fixed.
'''

#options['modelname'] = model.description  #Check this.

''' #Fix this.
if verbflag == 1:
    print model

'''

''' #Fix this
rawmodel = model
'''
#Set constraints:

#Options (Refer to Schuetz et al. 2007):
