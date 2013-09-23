#runsim.py
import numpy as np
import matplotlib.pylab as pl
import scipy.stats as st
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.io.sbml import write_cobra_model_to_sbml_file
from constrainfluxes import *

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

if ('verbflag' in options) == False:
    options['verbflag'] = 1 #By default, show output about the process as the program runs.

verbflag = options['verbflag']

if ('debugmode' in options) == False:
    options['debugmode'] = 0

if options['debugmode'] == 1: #Keyerror. Fix this.
    options['verbflag'] = 1


if ('model_id' in options) == False:
    options['model_id'] = 1
if verbflag == 1:
    print 'Runsim.py: ModelID not specified in input. Set to 1 (SCHUETZR) by default'

if ('constraints' in options) == False:
    options['constraints'] = 0
    if verbflag == 1:
        print 'Constraints not specified in input. Set to default model constraints.'

if ('objective' in options) == False:
    options['objective'] = 1
    if verbflag == 1:
        print 'Objective not specified in input. Set to biomass by default.'

if ('exp_id' in options) == False:
    options['exp_id'] = 1
    if verbflag == 1:
        print 'Experimental condition to compare with not specified in input. Set to aerobic batch culture by default.'

if ('optreq' in options) == False:
    options['optreq'] = 1
    if verbflag == 1:
        print 'Optimality requirement not specified in input. Set to 1 by default.'

if ('usesplits' in options) == False:
    options['usesplits'] = 0
    if verbflag == 1:
        print 'Use of raw fluxes or split ratios not specified in input. Set to raw fluxes by default.'

'''
if options['usesplits'] == 2: #Fix this.
    options['substractBM'] = 1
'''

if ('FBAonly' in options) == False:
    options['FBAonly'] = 0

if ('iterations' in options) == False:
    options['iterations'] = 0
    if verbflag == 1:
              print 'Number of random starting points for fmincon not specified. Set to zero by default (Use only FBA result).'

if ('usegurobi' in options) == False:
    options['usegurobi'] = 1

if ('maxwithgurobi' in options) == False:
    options['maxwithgurobi'] = 0 #%Change this if sucessful in getting maximization with Gurobi to work.

''' <<<This code doesn't do anything. Evaluate for correctness/necessity.>>>
if ('usegurobi' in options == 0) & ('usefmincon' in options == False):
    if verbflag == 1:
              print 'Solver not specified. Set to fmincon by default.'
'''

if ('usefmincon' in options) == False:
    if options['model_id'] == 3:
        options['usefmincon'] = 0
    else:
        options['usefmincon'] = 1

if ('preloadmodels' in options) == False:
    options['preloadmodels'] = 0 #Change to 1 when preloading is implemented.
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

if ('useweights' in options) == False:
    options['uncertweights'] = 0

if ('makeplots' in options) == False:
    options['makeplots'] = 0

if ('usescaledfluxes' in options) == False:
    options['usescaledfluxes'] = 0

if ('fluxreport' in options) == False:
    options['fluxreport'] = 1

if ('computercorrcoef' in options) == 0:
    options['computecorrcoef'] = 1

if ('normchoice' in options) == False:
    options['normchoice'] = 'euclidean'

if ('usealloptions' in options) == False:
    options['usealloptions'] = 0

if ('printfile' in options) == False:
    options['writefile'] = 0

if ('usetokens' in options) == False:
    options['usetokens'] = 1

if ('plotFBA' in options) == False:
    options['plotFBA'] = 0

if ('extractmode' in options) == False:
    options['extractmode'] = 'default'

if ('fluxconstrainset' in options) == False:
    options['fluxconstrainset'] = 0

if ('reportsplits' in options) == False:
    options['reportsplits'] = 0

if ('constrainglucose' in options) == False:
    options['constrainglucose'] = 2

if ('constrainBM' in options) == False:
    options['constrainBM'] = 1

if ('usealtmap' in options) == False:
    options['usealtmap'] = 0

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

#Fix this.
model_id = options['model_id'] #Fix this
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

#Fix this
if model_id == 1:
    if verbflag == 1:
              print 'Using Schuetz revised (SCHUETZR) model.'
    if options['preloadmodels'] != 1:
        model = create_cobra_model_from_sbml_file('../SBML/SCHUETZR.xml')
        print 'model:',model
        n = raw_input('Press any key to proceed')
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


#options['modelname'] = model.description  #Check this.

''' #Fix this.
if verbflag == 1:
    print model

'''

 #Fix this
rawmodel = model

#Set constraints:

#Options (Refer to Schuetz et al. 2007):


if constraints == 1:
    options['constraintsname'] = 'P/O = 1'
    if model_id == 1:
        #Invoke P/O ratio = 1 in fixed model:
        POobjrxns = ['nuo','cydB']
        if verbflag == 1:
            print 'Constraint set: P/O ratio = 1'
    else:
        raise ValueError('Model/Objective/Constraints combination not supported')
elif constraints == 2:
    #Invoke glucose/O2 uptake ratio of 2/3 in fixed model:
    options['constraintsname'] = 'Glucose/O2 uptake ratio = 2/3'
    if model_id == 1:
        #glclimit = model.lb(2) #Replace by cobrapy equivalent.
        glxO2rxns = ['EX_O2(e)','EX_GLC(e)']
        for rxn in glxO2rxns:
            model.reactions.get_by_id(rxn).lower_bound = 0
        model.reactions.get_by_id('EX_GLX_O2(e)').lower_bound = glclimit/2 #Activate coupled uptake reaction
        if verbflag == 1:
            print 'Constrain set: Glucose/O2 uptake ratio = 2/3'
    elif model_id == 2:
        raise ValueError('Model/Objective/Constraints combination not supported')
    elif momdel_id == 3:
        raise ValueError('Model/Objective/Constraints combination not supported')

elif (constraints == 3) or (constraints == 4):
    if constraints == 3:
        O2limit = -11.5
        options['constraintsname'] = 'O2<11.5'
    elif constraints == 4:
        O2limit = -14.75
        options['constraintsname'] = 'O2<14.75'

    if model_id == 1:
        #Invoke qO2_max < 11.5 in SCHUETZR:
        R1 = model.reactions.get_by_id("EX_O2_e")  
    elif model_id == 2:
        R1 =  model.reactions.get_by_id("EX_O2_e")
    elif model_id == 3:
        R1 =  model.reactions.get_by_id("EX_O2_e")
    R1.lower_bound = O2limit
    if verbflag == 1:
        print 'Constraint set: Max O2 uptake rate = %.2f mmol/g*h' % (-O2limit)

elif constraints == 5:
    options['constraintsname'] = 'ATP maint.'
    if model_id == 1:
        R1 = model.reactions.get_by_id('maint')
    elif model_id == 2 or model_id == 3:
        R1 = model.reactions.get_by_id('ATPM')
    R1.lower_bound = 7.6
    R1.upper_bound = 7.6
    if verbflag == 1:
        print 'Constraint set: ATP maintenance requirement = 7.6 mmol/g*h'
        

elif constraints == 5.1:
    if model_id == 1:
        R1 = model.reactions.get_by_id('maint')
    elif model_id == 2 or model_id == 3:
        R1 = model.reactions.get_by_id('ATPM')
    R1.lower_bound = 7.6

elif constraints == 6:
    options['constraintsname'] = 'Rx bounds'
    #Set max values for all fluxes in fixed model to 200 % of maximum glucose uptake rate.
    #(Slightly different than described by Schuetz et al.)
    if model_id == 1:
        limit = -2*model.reactions.get_by_id('EX_GLC_e').lower_bound

        for reaction in model.reactions:
            if ('EX_' in reaction.id == False): #If the reaction is not an exchange reaction. (Assuming that all exchange reactions and no others contain 'EX' in the ID)
                reaction.upper_bound = limit
                if reaction.lower_bound != 0:
                    reaction.lower_bound = -limit
    else:
        raise ValueError('Model/Objective/Constraints combination not supported')
    if verbflag == 1:
        print 'Constraint set: Bounds on all reactions set to 200 % of max glucose uptake rate'

elif constraints == 7: #Combination of previous constraints
    options['constraintsname'] = 'All'
    if model_id == 1:
        #Case1
        POobjrxns = ['nuo','cydAB']
        for rxn in POobjrxns:
            model.reactions.get_by_id(rxn).lower_bound = 0
        #Case 2:
        glclimit = model.reactions.get_by_id('EX_GLC_e').lower_bound
        glcO2rxns = ['EX_O2(e)','EX_GLC(e)']
        for rxn in glcO2rxns:
            model.reactions.get_by_id(rxn).lower_bound = 0 #Deactivate normal glucose and O2 uptake
        model.reactions.get_by_id('EX_GLC_O2(e)').lower_bound = glclimit/2

        #Case 3:
            #Overriden by case 2.

        #Case 5:
        model.reactions.get_by_id('maint').upper_bound = 7.6
        model.reactions.get_by_id('maint').lower_bound = 7.6
        #Case 6:
        limit = -2*model.reactions.get_by_id('EX_GLC_O2(e)').lower_bound
        for rxn in model.reactions:
            if ('EX_' in reaction.id == False):
                reaction.upper_bound = limit
                if reaction.lower_bound != 0:
                    reaction.lower_bound = -limit
    else:
        raise ValueError('Model/Objective/Constraints combination not supported')
    if verbflag == 1:
        print('Constraint set: All constraints.')
elif constraints == 8: #Combination of constraints 3,5 and 6:
    #Constraint 3:
    if model_id == 1:
        #Invoke qO2_max < 11.5 in SCHUETZR
        model.reactions.get_by_id('EX_O2(e)').lower_bound = -11.5
        model.reactions.get_by_id('maint').upper_bound = 7.6
        model.reactions.get_by_id('maint').lower_bound = 7.6
    if model_id == 2:
        model.reactions.get_by_id('EX_O2(e)').lower_bound = -11.5
        model.reactions.get_by_id('ATPM').upper_bound = 7.6
        model.reactions.get_by_id('ATPM').lower_bound = 7.6
elif constraints == 0:
    options['constraintsname'] = 'Default'
    model = rawmodel
    print 'model:',model
    print 'rawmodel:',rawmodel
    if verbflag == 1:
        print 'Using default model constraints'

#Define growth-media:

if exp_id == 1:
    if verbflag == 1:
        print 'Simulating aerobic culture. Maximal uptake rates at default values.'
elif exp_id == 2:
    model.reactions.get_by_id('EX_o2(e)').upper_bound = 0
    model.reactions.get_by_id('EX_o2(e)').lower_bound = 0
    if verbflag == 1:
        print 'Simulating anaerobic culture. Maximal oxygen uptake set to zero.'
elif exp_id == 3:
    raise ValueError('Model/Objective/Constraints combination not supported')

#Set growth rate to experimentally determined value

BMreactnum = 25 # #of the biomass reactionin in the experimental dataset. This line should be changed to read in the value. Hardcoded for now

if constrainBM == 1:
    print 'runsim.py: Constraining growth rate to experimentally determined value.'
    tempopt = options
    tempopt['fluxconstrainset'] = BMreactnum
    #model = constrainfluxes(model,tempopt) #call constrainfluxes.py #Uncomment when constrainfluxes() is ready.

#Set objective:


FBArunflag = 0 #Change to 1 when FBA simulation has been run.

#Parameters for optimizeCbModel, and other relevant flags

options['objective_sense'] = 'maximize'
#minnorm = 0 # By default, minimization of fluxes is turned off. #Haven't found an optinon for this isn CobrapY.

options['objectivetype'] = 'linear' #Describe what kind of objective is optimized. for use in QPmindist.py

if objective == 1:
    objectivename = 'BM'
    if model_id == 1:
        model.reactions.get_by_id('biomass').lower_bound = 0
        model.reactions.get_by_id('biomass').upper_bound = 1000
    elif model_id == 2:
        model.reactions.get_by_id('biomass_Ecoli_core_w_GAM').lower_bound = 0
        model.reactions.get_by_id('biomass_Ecoli_core_w_GAM').upper_bound = 1000
    elif model_id == 3:
        model.reactions.get_by_id('Ec_biomass_iJO1366_core_53p95M').lower_bound = 0
        model.reactions.get_by_id('Ec_biomass_iJO1366_core_53p95M').upper_bound = 1000
    if verbflag == 1:
        print 'Objective set to max biomass production. Constraints on biomass have been lifted.'
    if constrainglucose > 0:
        if verbflag == 1:
            print 'runsim.m constraining glucose flux.\n'
        tempopt = options
        tempopt['fluxconstrainset'] = 1
        #model = constrainfluxes(model,tempopt) #Uncomment when constrainfluxes is ready.
        if debugmode < 0:
            header = ['reaction #','Reaction name','LB','UB']
            print '{0:10}\t{1:10}\t{2:10}\t{3:10.}'.format(*header)
            for i in len(model.reactions):
                print '{0:10}\t{1:10}\t{2:10}\t{3:10]'.format(i,model.reactions[i].name,model.reactions[i].lower_bound,model.reactions[i].upper_bound)
    if constrainglucose > 1:
        #load(expdata.mat)
        #glucoseconsumption = expdata.perrenoud.abs.batch.aerobe.glucoseconsumption; #Uncomment when data import is ready.
        if verbflag == 1:
            print 'runsim.m: Constraining glucose uptake flux.'
        if model_id == 1:
            rn = 'EX_GLC(e)'
        elif (model_id == 2) or (model_id == 3):
            rn = 'EX_glc(e)'
                
        #matlab code: model = changeRxnBounds(model,rn,-glucoseconsumption(1)-glucoseconsumption(2),'l');
        #matlab code: model = changeRxnBounds(model,rn,-glucoseconsumption(1)+glucoseconsumption(2),'u');

        #model.reactions.get_by_id(rn).lower_bound = -glucoseconsumption[0]-glucoseconsumption[1] #Uncomment when glucoseconsumption is defined.
        #model.reactions.get_by_id(rn).upper_bound = -glucoseconsumption[0]+glucoseconsumption[1] #                                   
elif objective == 2:
    objectivename = 'max ATP'
    if model_id == 1:
        atpOBJrxns =['pgk','pykAF','sucCD','atp','ackAB_tdcD_purT']
    elif model_id ==2:
         atpOBJrxns =['ATPS4r','Pyk']
    elif model_id == 3:
        atpOBJrxns =['AP5AH','ATPS4rpp','PYK']
    cobra.flux_analysis.objective.update_objective(model, atpOBJrxns)
    if verbflag == 1:
        print 'Objective reactions set: ATP producing fluxes'
elif objective == 3.1:
    objectivename='Flux minimization (euclidean norm)';
    if verbflag == 1:
        print 'Objective is minimization of overall intracellular flux (euclidean norm)'
        print 'Performing custom FBA...'
        minnorm = 1e-6
        options['objectivetype'] = 'euclidean'
        import warnings as wn
        wn.warn('RUNSIM:UnsupportedObjective','Untested objective. Proceed with caution.')
elif objective == 3.2:
    objectivename='Flux minimization (manhattan norm)'
    if verbflag == 1:
            print 'Objective is minimization of overall intracellular flux (taxicab norm)'
            print 'Performing custom FBA...'
            minnorm = 'one'
            options['objectivetype'] = 'taxicab'
elif objective == 4:
    objectivename='min glucose'
    if model_id == 1:
        glcRXNs=['EX_GLC(e)','EX_GLC_O2(e)']
    elif model_id == 2 or model_id == 3:
        glcRXNs=['EX_glc(e)']
    cobra.flux_analysis.objective.update_objective(model,glcRXNs)
    options['objective_sense'] = 'maximize'

    #lift constraints on glucose flux:
    print 'Lifting constraints on glucose flux.'
    if model_id == 1:
        glucoseID = 'EX_GLC(e)'
        
    elif model_id ==2 or model_id == 3:
        glucoseID = 'EX_glc(e)'

    model.reactions.get_by_id(glucoseID).upper_bound = rawmodel.reactions.get_by_id('EX_GLC(e)').upper_bound
    model.reactions.get_by_id(glucoseID).lower_bound = rawmodel.reactions.get_by_id('EX_GLC(e)').lower_bound
    if verbflag == 1:
        print 'Objective set: Minimization of glucose fluxes.'
elif objective == 5:
    objectivename= 'Min redox'
    if model_id == 1:
        NADHrxns=['aceEF','maeA','sucAB','udhA','fdhF_fdoGHI','ldhA','maeB','zwf','gnd','pntAB','frdABCD_sdhAB','dld']
        cobra.flux_analysis.objective.update_objective(model,NADHrxns)
        options['objective_sense'] = 'minimize'
    else:
        raise ValueError('Objective not supported.')
    if verbflag == 1:
        print 'Objective reactions set. Minimization of redox reactions.'

elif objective == 6:
    objectivename = 'Min ATP'
    if model_id == 1:
        atpOBJrxns = ['pgk','pykAF','sucCD','atp','ackAB_tdcD_purT']
    elif model_id == 2:
        atpOBJrxns =['ATPS4r','Pyk']
    elif model_id == 3:
        atpOBJrxns =['AP5AH','ATPS4rpp','PYK']
    cobra.flux_analysis.objective.update_objective(model,atpOBJrxns)
    options['objective_sense'] = 'minimization'
    if verbflag == 1:
        print 'Objective reactions set. Minimization of ATP fluxes.'

#Optional: Constrain one or some fluxes:

if fluxconstrainset != 0:
    #model = constrainfluxes(model,options) #Uncomment when constrainfluxes is ready
    pass

print 'MODEL:',model
print 'RAWMODEL:',rawmodel

options['objectivename'] = objectivename
#Perform FBA:

if ('objective_sense' in options) == False:
    options['objective_sense'] = 'maximize' 
objective_sense = options['objective_sense']

if FBArunflag == 0:
    if verbflag == 1:
        print 'Performing FBA...'
        #print '%s %s \n' %('Model:',options['modelname']) #Uncomment when modelname is defined.
        print '%s %s \n' %('Constraints:',options['constraintsname'])
        print '%s %s \n%s %s \n%s %s \n' %('Objective type:',
                                           options['objectivetype'],
                                           'Objective:',options['objectivename'],
                                           'Sense:',options['objective_sense'])
        model.optimize(objective_sense = objective_sense)
if debugmode == 1:
    print 'RBArunflag:'
    print FBArunflag
    print 'Solution status:',model.solution.status
    print 'Objective value:',model.solution.f

#Send distance from FBA result to output

tempoptions = options
tempoptions['excludereactions'] = 0 #In calculating the distance of the optimal reaction, don't exclude any reactions.

model.solution.constraintsdescription = options['constraintsname']
#model.solution.FBAdistance = compdist(model.solution.x,tempoptions,1) #Uncomment when compdist is finished
#model.solution.FBAsplits = computesplits(model.solution.x,model_id) #Uncomment when compdist is finished

# model.solution.expvalues = refexpfluxes #Add reference experimental fluxes to output. #Uncomment when refexpfluxes is ready.

def getobjvec(model):
    c =[]
    for i in range(len(model.reactions)):
        c.append(model.reactions[i].objective_coefficient)
    return c

model.solution.c = getobjvec(model)

model.solution.modelname = model.description

#Search for distance-optimal solution, unless input specified that only FBA should be performed:

if FBAonly != 1:
    if verbflag == 1:
        print 'Starting minimizing/maximization of distance procedure(s)'

    #Call fmindist to minimize/maximize the distance:
    if usefmincon == 1:
        if verbflag == 1:
            print 'Minimizing/maximizing distance of raw fluxes using fmincon...Starting point: FBA result.'
        if usesplits == 1:
            print 'Minimizing/maximizing distance of split ratios using fmincon...Starting point: FBA result.'
            splitresult ={}
            tempopt = options
            tempopt['optsplits'] = 1
            tempopt['compsplits'] = 1
            splitresult['splitminsolution'] = fmindist(model.solution,model,tempopt,1)
            splitresult['splitmaxsolution'] = fmindist(model.solution,model,tempopt,-1)
            splitresult['expsplitvalues']= computesplits(model.solution.expvalues,0)
            splitresult['splitminvalues'] = computesplits(splitresult['splitminsolution'],model_id)
            splitresult['splitmaxvalues'] = computesplits(splitresult['splitmaxsolution'],model_id)
            tempopt['compslits'] = 1
            splitresult['splitmindistance'] = compdist(splitresult['splitminsolution'],tempopt,1)
            splitresult['splitmaxdistance'] = compdist(splitresult['splitmaxsolution'],tempopt,1)
            splitresult['splitrawmindistance'] = compdist(splitresult['splitminsolution'],options,1)
            splitresult['splitrawmaxdistance'] = compdist(splitresult['splitmaxsolution'],options,1)
            #%Not currently working
        options['optsplits'] = 0
        model.solution.Fmin_minsol = fmindist(model.solution,model,options,1)
        model.solution.Fmin_maxsol = fmindist(model.solution,model,options,-1)

        tempoptions = options
        tempoptions['excludereactions'] = 0

        model.solution.Fmin_mindistance = compdist(model.solution,Fmin_minsol,tempoptions,1)
        model.solution.Fmin_maxdistance = compdist(model.solution,Fmin_maxsol,tempoptions,1)
        model.solution.Fmin_minsol_exp = extractflux(model.solution.Fmin_minsol,options)
        model.solution.Fmin_maxsol_exp = extractflux(model.solution.Fmin_maxsol,options)

        #%Use fmincon iteratively with a specified number of random starting points:
        if iterations != 0:
            if verbflag == 1:
                print 'Running fmincon with multiple random starting points.'
                print 'Number of iterations:',iterations
            k = iterations
            randresultmatrix = np.zeros(len(model.reactions),k)
            randminsolutions = np.zeros(len(model.reactions),k)
            randmaxsolutions = np.zeros(len(model.reactions),k)
            randmindistance = np.zeros(1,k)
            randmaxdistance = np.zeros(1,k)
            #Make a copy of the model so the objective can be changed without affecting other operations:
            randmodel = model

            def applyobjective(c,model):
                if len(c) != len(model.reactions):
                    raise ValueError('Number of reactions do not match.')
                else:
                    for i in range(len(c)):
                        model.reactions[i].objective_coefficient = c[i]
            for i in range(k):
                randmodel.c = np.random.random((len(model.reactions),1))
                randmodel.optimize(solver='gurobi') #Perform FBA with the random objective.
                #Matlab code: randresultmatrix(:,1) = randresult.x # Record the resulting flux vector
                #Matlab code:randminsolutions(:,1) = fmindist(randresult,model,options,1)
                #Matlab code: randmaxsolutions(:,i) = fmindist(randresult,model,options,-1)
                #matlab code: model.solution.randmindistance(i) = compdist(randminsolutions(:,1),options,1)
                #matlab code: model.solution.randmaxdistance(i) = compdist(randmaxsolutions(:,i),options,1)
                
            #Compute further outputs:

            model.solution.meanrandmindistance = mean(model.solution.randmindistance)
            model.solution.meanrandmaxdistance = mean(model.solution.randmaxdistance)
            model.solution.stddevrandmindistance = std(model.solution.randmindistance)
            model.solution.stddevrandmaxdistance = std(model.solution.randmaxdistance)
            #matlab code: result.bestrandminsolution = randminsolutions(:,bestminindex);
            #matlab code: result.bestrandmaxsolution = randmaxsolutions(:,bestmaxindex);
    
        if makeplots == 1:
            if verbflag == 1:
                print 'Generating plot: Minimal and maximal distance found with fmincon using starting points generated by FBA using random objective functions.'
        xlabel = 'Objective #'
        ylabel = 'Distance'
        fig = pl.figure()
        left = [i for i in range(len(model.solution.randmindistance))]
        ax = fig.add_subplot(1,1,1)
        ax.bar(left,[i for i in model.solution.randmindistance],align = 'center')
    #END of fmincon iterative routine.

    if usegurobi == 1:
        if verbflag == 1:
            print 'Minimizing distance of raw fluxes using Gurobi.'
        gurobi_output  = QPmindist(model,result,options,1)
        #Operations on gurobi output goes here.
#%end of "if FBAonly~=1..."


if verbflag == 1:
    print 'Summary:'
    print 'Relative optimality required is:',optreq
    if usefmincon == 1:
        print 'Minimal distance found using fmincon (starting point: FBA result)',result.Fmin_mindistance
        print 'Maximal distance found using fmincon (starting point: FBA result):',result.Fmin_maxdistance
    if usegurobi == 1:
        print 'Minimal achievable distance found using Gurobi is:',result.gurobi_mindistance
        #print 'Maximal achievable distance found using Gurobi is:',result.gurobi_maxdistance
    if iterations != 0:
        print 'Minimzal achievable distance found using fmincon (Random starting points) is:',result.bestandmindistance
        print 'Maximal achievable distance found using fmincon (Random starting points) is:',result.bestrandmaxdistance
        print 'Number of random objectives used to generate starting points:',iterations

#End of "if verbflag == 1..:"

#Compute correlation coefficient between computed and experimental fluxes:

if computecorrcoef == 1:
    if hasattr(result,'Fmin_minsolution') == True:
        #matlab/python mashup: result.fmincorrcoef = st.pearson(result.expvalues(:,1),extractflux(result.Fmin_minsol,options))
        #matlab code: result.mincorrcoef = result.fmincorrcoef(2,1)
        pass
    if hasattr(result,'gurobi_minsolution') == True:
        #matlab code: result.gurobimincorrcoef = corrcoef(result.expvalues(:,1),extractflux(result.gurobi_minsol,options));
        #matlab code: result.gurobimincorrcoef = result.gurobimincorrcoef(2,1);
        pass
    if hasattr(result,'bestrandminsolution') == 1:
        #matlab code: result.bestrandmincorrcoef = corrcoef(result.expvalues(:,1),extractflux(result.bestrandminsolution,options));
        #matlab code: result.bestrandmincorrcoef = result.bestrandmincorrcoef(2,1);
        pass
    if verbflag == 1:
        print 'Correlation coefficients:'
        if hasattr(result,'fmincorrcoef') == True:
            print 'Fmincon solution (starting point FBA result):',result.fmincorrcoef
        if hasattr(result,'bestrandmincorrcoef') == True:
            print 'Fmincon solution (best of solutions generated from random starting points):',result.bestrandmincorrcoef
        if hasattr(result,'gurobimincorrcoef') == True:
            print 'Gurobi solution:',result.gurobimincorrcoef
#End of "if computcorrcoef..."

if exp_id == 1:
    #matlab code: expfluxvalues = expdata.perrenoud.abs.batch.aerobe.fluxvalues
    #matlab code: expsolution = expfluxvalues(:,1)
    #matlab code:experrors = expfluxvalues(:,2)
    pass

#Generate plot:
#Plot the best fit solution(s) together with experimental values and uncertainties.
#Also include FBA solution?

if makeplots == 1:
    if verbflag == 1:
        print 'Generating result plot...'
    fig = pl.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_title('Experimental and computed reactionrates')
    legendstring = []
    if hasattr(result,'Fmin_minsol') == True:
        fluxresult = extractflux(result.Fmin_minsol,options)
        pl.scatter([i for i in range(len(fluxresult))],fluxresult)
        legendstring.append('Fmincon solution')
    if hasattr(result,'gurobi_minsol') == 1:
        fluxresult = extractflux(result.gurobi_minsol,options)
        pl.scatter([i for i in range(len(fluxresult))],fluxresult)
        legendstring.append('Quadratic minimization')
    if options['plotFBA'] == 1:
        fluxresult = extractflux(result.x,options)
        pl.scatter([i for i in range(len(fluxresult))],fluxresult)
        legendstring.append('FBA solution')
    #xlabel('Experimental reaction #')
    #ylabel('Flux value (mmol/g*h)
    #legend(legendstring)
#End of "if makeplots == 1..."

#Compute objective function values of generated solutions:

if hasattr(result,'Fmmin_minsol') == True:
    #matlab code:result.Fmin_minsol_objval = model.c'*result.Fmin_minsol
    pass

if hasattr(result,'gurobi_minsol') ==True:
    #matlab code:result.gurobi_minsol_objval = model.c'*result.gurobi_minsol
    pass

#Generate flux report:

if options['fluxreport'] == 1:
    if hasattr(result,'Fmin_minsol') == True:
        tempopt = options
        tempopt['solvername'] = 'fmincon'
        fluxresult = extractflux(result.Fmin_minsol,tempopt)
        result.fluxresult = fluxresult
        result.fluxreport.fmincon = fluxreport(fluxresult,expfluxvalues,options)
    if hasattr(result,'gurobi_minsol') == 1:
        tempopt = options
        tempopt.solvername = 'gurobi'
        gurobifluxresult = extractflux(result.gurobi_minsol,options)
        result.fluxreport.gurobi = fluxreport(gurobifluxresult,expfluxvalues,tempopt)

#Compare fmincon and Gurobi results:

if (hasattr(result,'Fmin_minsol') == 1) and (hasattr(result,'gurobi_minsol') == True):
    Fminexp = extractflux(result.Fmin_minsol,options)
    Gurobiexp = extractflux(result.gurobi_minsol,options)
    diff = Fminexp-Gurobiexp
    totdiff = sum(diff)
    result.compare = [Finexp,Gurobiexp,diff]
    if options['verbflag'] == 1:
        print 'Comparison of Fmincon and Gurobi solutions:\n'
        print 'Excluded reactions (reaction IDs):',options['excludedreactions']
        textheader = ['Reaction #','#Fmincon','Gurobi','Difference']
        print '{0:10}\t{1:10}\t{2:10}\t{3:10.}'.format(*textheader)
        #Code for printing out table
    if makeplots > 0:
        fig = pl.figure()
        ax = fig.add_subplot(1,1,1)
        ax.set_title('Differences in predicted reaction rates')
        #matlab code: title('Differences in predicted reaction rates')
        #matlab code: xlabel('Experimental reaction #')
        #matlab code: ylabel('Difference (mmol/g*h)')
        #pl.show()

#Select the best (most distance-optimal) solution of those calculated

'''
MATLAB CODE:

k = 1;
if isfield(result,'Fmin_minsol') == 1
    resultarray(k).name = 'Fmincon';
    resultarray(k).solution = result.Fmin_minsol;
    resultarray(k).mindist = result.Fmin_mindistance;
    k = k+1;
end

if isfield(result,'gurobi_minsol') == 1
   resultarray(k).name = 'Gurobi';
   resultarray(k).solution =  result.gurobi_minsol;
   resultarray(k).mindist = result.gurobi_mindist;
   k = k+1;
end

if isfield(result,'FBAdistance') == 1
     resultarray(k).name = 'FBA';
     resultarray(k).solution = result.x;
     resultarray(k).mindist = result.FBAdistance;
   
end
    
mindists =  {resultarray(:).mindist};
mindists =  cell2mat(mindists);
[C,I] = min(mindists); %Find the smallest distance
equiv = zeros(length(mindists)); %Vector to store indexes of equally optimal solutions

if isfield(options,'devtol') == 0
    options.devtol = 1e-2;
end

devtol = options.devtol;
maxdev = 0;
for i = 1:length(mindists)
    if i ~=I
        dev = C - mindists(i);
        if abs(dev) < devtol
            equiv(i) = i;
            if abs(dev) > maxdev
                maxdev = abs(dev);
            end
        end
    end
end

    if verbflag > 0
 
        if nnz(equiv) == 0
            fprintf(1,'%s %s%s \r\n','The most distance-optimal solution was generated by',resultarray(I).name,'.');
            
        elseif nnz(equiv) > 0
            fprintf(1,'%s %s','Equivalent solutions were found by',resultarray(I).name);
            for i = 1:length(equiv)
                if equiv(i) > 0
                    fprintf(1,'%s %s',' and',resultarray(equiv(i)).name);
                   
                       
                end
            end
             fprintf(1,'%s \r\n','.');
             fprintf(1,'%s %.4f \r\n','Maximal distance deviation between equivalent solutions:',maxdev);
             fprintf(1,'%s %.4f \r\n','Maximal distance deviation to be considered equivalent solutions:',devtol);
            
        end
            
            fprintf(1,'%s %.2f %s \n','Minimal distance achieved:',C,'mmol/g*h');
            fprintf(1,'%s %.1f \r\n','Optimality requirement:',options.optreq);
    end
   

   result.options = options;
   result.model = model;

(end of runsim.m)
'''
