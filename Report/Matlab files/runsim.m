%runsim.m
%Perform a simulation to find minimum (and in some cases, maximum) distance
%between computed and experimental fluxes given selected model constraints and objective-optimality requirement:
%
%syntax: result = runsim(options)
%
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
%Executes the program with all selections set to the default.
%To set options, use the syntax options.field = fieldvalue 
%Example: options.objective = 2;
%
%For details on constraints, objectives and experimental data, see Schuetz
%et al. 2007: Systematic evaluation of objective functions.
%
%Select model by setting the model_id field. Options:
%SHUETZR: 1
%ECME: 2
%iJO1366: 3 
%
%Select constraints by setting the constraints field. Options
%0: Default model constraints (default)
%1: P/O ratio = 1
%2: Glucose/oxygen uptake ratio = 2/3
%3: qO2_max < 11.5
%4: qO2_max < 14.75:
%5:  ATP Maintenance
%5.1: ATP maintenance, additional waste allowed
%6: Bounds (UB set to 200 % of maximal glucose uptake for all fluxes)
%7: Combination of above (1, 2, 3, 5 & 6)
%
%Select experimental fluxes to compare with by setting the exp_id field. Options: 
%1: Aerobic batch (default)
%2: Anaerobic batch, NO3 limited
%3: Anaerobic batch
%
%Set objective function by setting the objective parameter.
%options:
%
%1: biomass (default)
%2: Max ATP production
%3.1: Minimize overall intracellular flux (euclidean norm)
%3.2: Minimize overall intracellular flux (taxicab norm)
%4: Min glc/biomass
%5: Min redox reactions 
%6: Min ATP production
%
%Set requirement for optimality of solution by setting the optreq
%field. Valid range is 0 (no requirement for optimality) to 1 (consider
%only optimal solutions).
%
%Determine if raw fluxes or split ratios should be used by setting the
%usesplits field. %Options:
%0: Raw fluxes (default)
%1: Split ratios
%2: Split ratios, BM substracted

function result = runsim(options)

if nargin <1
    options = struct();
    options.verbflag = 1;
end

result = struct();

if isfield(options,'verbflag') == 0;
    options.verbflag = 1; %By default, show output about the process as the program runs.
end

verbflag = options.verbflag;

if isfield(options,'debugmode') == 0
    options.debugmode = 0;
end

if options.debugmode == 1 %If debug mode is on, activate verbose mode.
    options.verbflag = 1;
end


if isfield(options,'model_id') ==0
    options.model_id = 1; %Default model is SCHUETZR
    if verbflag == 1
        disp('Runsim.m: Model ID not specified in input. Set to 1 (SCHUETZR) by default.')
    end
end

if isfield(options,'constraints') == 0
    options.constraints = 0; %Default constraints are default model constraints
    if verbflag == 1
        disp('Constraints not specified in input. Set to default model constraints by default.')
    end
end

if isfield(options,'objective') == 0
    options.objective = 1; %Default objective is maximization of biomass reaction
    if verbflag == 1
        disp('Objective not specified in input. Set to biomass by default')
    end
end

if isfield(options,'exp_id') == 0
    options.exp_id = 1; %Default experimental condition is batch aerobe.
    if verbflag ==1
        disp('Experimental condition to compare with not specified in input. Set to batch aerobe by default.')
    end
end

if isfield(options,'optreq') == 0;
    options.optreq = 1; %By default, consider only optimal solutions.
    if verbflag == 1
       disp('Optimality requirement not specified in input. Set to 1 (consider only optimal solutions) by default.')
    end
end

if isfield(options,'usesplits') == 0
    options.usesplits = 0; %By default, use raw fluxes, not split ratios.
    if verbflag == 1
     disp('Use of raw fluxes or split ratios not specified in input. Set to raw fluxes by default.')
    end
     
end

if options.usesplits == 2
    options.substractBM = 1;
end

if isfield(options,'FBAonly') ==0
    options.FBAonly = 0; %Perform full analysis by default. Elsewise, perform only standard FBA.
end

if isfield(options,'iterations') == 0
    options.iterations = 0; %by default, fmincon is only run once.
    if verbflag == 1
        disp('Number of random starting points for fmincon not specified. Set to zero by default (use only FBA result).')
    end
end

if isfield(options,'usegurobi') == 0
    options.usegurobi = 1; 
end

if isfield(options,'maxwithgurobi') == 0
    options.maxwithgurobi =0;
end
%Change this if sucessful in getting maximization with Gurobi to work.

if isfield(options,'usegurobi') == 0 && isfield(options,'usefmincon')
    if verbflag == 1
        disp('Solver not specified. Set to fmincon by default.')
    end
end

if isfield(options,'usefmincon') == 0
    if options.model_id == 3
       options.usefmincon = 0;
    else
        options.usefmincon =1;
    end
end

if isfield(options,'preloadmodels') ==0
    options.preloadmodels =1;
    %Use to specify if model structure should be loaded from a .mat file, or be constructed from a SBML file using the COBRA function readCbModel.
    %By default, use pre-loaded model data. Remember to update model-data
    %if changes are made to the models!
end


%Allow the posssibility of excluding some reactions from optimization
%procedures
if isfield(options,'excludereactions') ==0 
    options.excludereactions = 0; %By default, no reactions are excluded from the analysis.
    %To exclude experimental reactions, set options.excludereactions = a
    %vector with the experimental reaction #s.
    %example: options.excludereactions = [1 12 24].
end


if isfield(options,'useweights') == 0
    options.uncertweights = 0;    
    %Weigh flux distances using the experimental uncertainties. This way, differences are more tolerated if the experimental uncertainty is high.
end

if isfield(options,'makeplots') == 0
    options.makeplots =0;
end

if isfield(options,'usescaledfluxes') ==0 
    %To implement later: Use glucose-scaled fluxes (relative fluxes) used in Schuetz 2007 as input for optimization with fmincon.
    options.usescaledfluxes=0;
end

if isfield(options,'fluxreport') ==0
    options.fluxreport = 1;  %Produce a report comparing experimental and computed fluxes
end


if isfield(options,'computecorrcoef') ==0
    options.computecorrcoef =1; %Compute correlation coefficients. 
end

if isfield(options,'normchoice') ==0
    options.normchoice='euclidean'; %Use euclidean distance measure as default. Second possibility: 'taxicab'
end


if isfield(options,'usealloptions') == 0
    options.usealloptions =0;
end

if isfield(options,'printfile') == 0
    options.writefile = 0;
end

if isfield(options,'usetokens') ==0
    options.usetokens = 1;
end

if isfield(options,'plotFBA') ==0
    options.plotFBA = 0; %Option to include FBA result in generated plots.
end

if isfield(options,'extractmode') == 0
    options.extractmode = 'default';
    %Set which data source to use when extracting flux data. Options: Hardcoded, simple reaction map ('simplemap') and reaction map with token reactions ('tokens'). 
end

if isfield(options,'fluxconstrainset') == 0
    options.fluxconstrainset = 0;
end

if isfield(options,'reportsplits') == 0
    options.reportsplits = 0;
end

if isfield(options,'constrainglucose') == 0
    options.constrainglucose = 2;
end

if isfield(options,'constrainBM') == 0
    options.constrainBM = 1;
end

if isfield(options,'usealtmap') == 0
    options.usealtmap = 0;
end


if options.usealloptions == 1 %Activate all optional features.
    options.fluxreport = 1;
    options.computecorrcoef=1;
    options.makeplots = 1;
   options.usegurobi =1;
   options.usefmincon =1;
   options.iterations = 10;
   options.verbflag =1;
end


model_id = options.model_id;
constraints = options.constraints;
objective = options.objective;
exp_id = options.exp_id;
optreq = options.optreq;
usesplits = options.usesplits;
FBAonly = options.FBAonly;
usegurobi = options.usegurobi;
iterations = options.iterations;
makeplots = options.makeplots;
debugmode = options.debugmode;
usefmincon = options.usefmincon;
computecorrcoef = options.computecorrcoef;
fluxconstrainset = options.fluxconstrainset;
constrainglucose = options.constrainglucose;
constrainBM = options.constrainBM;
%Set options to default values if not specified:

if options.preloadmodels == 1
    if verbflag == 1
        disp('Loading COBRA model data...')    
    end
      load('modeldata.mat')
end


%Specify 'default' set of excluded reactions:
if strcmp(options.excludereactions,'default') == 1
    switch model_id
        case 1
            options.excludereactions = [1 12 18 19 22 23 24 25];
        case 2
            options.excludereactions = [1 11 12 19 22 23 24 25];
        case 3
            options.excludereactions = [1 11 12 18 19 22 23 24 25];
    end
end

 if verbflag == 1
    fprintf(1,'%s','runsim.m: Excluding following experimental reactions: ');
    for i = 1:length(options.excludereactions)
        if i < length(options.excludereactions)
        fprintf(1,'%d%s',options.excludereactions(i),',');
        else
        fprintf(1,'%d%s \n',options.excludereactions(i),'.');    
        end
    end
  
    %disp('runsim.m: Excluding following experimental reactions:')
    %disp(options.excludereactions)
 end


if verbflag == 1
        disp('Loading experimental data...')
end
load('expdata');

if verbflag == 1
    disp('Setting reference experimental fluxes...')
end


%Load the experimental flux data:
switch exp_id
    case 1
        refexpfluxes = expdata.perrenoud.abs.batch.aerobe.fluxvalues;  %Cannot add to "result" structure yet, because this gets overwritten by FBA analysis result below.
        
        if verbflag ==1
            disp('Reference experimental fluxes set: Perrenoud 2005 batch aerobe')
        end
        
    case 2
                disp('Experimental flux data (absolute values) not available')
    case 3
                disp('Experimental flux data (absolute values) not available')
end


if debugmode ~=1
    warning('OFF', 'MATLAB:strrep:InvalidInputType') %Turn off warnings from COBRA regarding model read-in.
end

if verbflag == 1
    disp('Options used are:')
    disp(options)
end

switch model_id

    case 1
        if verbflag == 1
        disp('Using Schuetz revised (SCHUETZR) model')
        end
        
        if options.preloadmodels ~= 1
        model = readCbModel('SCHUETZR');
        else model = modeldata.SCHUETZR;
        end
    case 2
        if verbflag == 1
        disp('Using expanded E coli core (ECME) model')
        end
        if options.preloadmodels ~= 1
        model = readCbModel('ECME');
        else model = modeldata.ECME;
        end
    case 3 
        if verbflag == 1
        disp('Using iJO1366 genome-scale model')
        end
        if options.preloadmodels ~= 1
        model = readCbModel('iJO1366b');
        else model = modeldata.iJO1366b;
        end
end

%Save model description in result and options structures:
options.modelname = model.description;


if verbflag == 1
    disp(model)
end


rawmodel = model; %Save a copy of the raw model without additional constraints

%Set constraints:

%Options (Refer to Schuetz et al. 2007):


%0: Reset constraints to default
switch constraints
        case 1
        options.constraintsname='P/O = 1';   
            switch model_id 
    
            case 1
            %Invoke P/O ratio = 1 in fixed model:
            POobjrxns = {'nuo','cydAB'};
            model = changeRxnBounds(model,POobjrxns,0,'b');
                if verbflag == 1
            disp('Constraint set: P/O ratio = 1')
                end
    
            case 2
            error('Model/Objective/Constraints combination not supported')     
            case 3
            error('Model/Objective/Constraints combination not supported')  
            end

    case 2
%Invoke glucose/O2 uptake ratio of 2/3 in fixed model:
    options.constraintsname='Glucose/O2 uptake ratio = 2/3';
    switch model_id
        
        case 1

    glclimit = model.lb(2);
    glcO2rxns = {'EX_O2(e)','EX_GLC(e)'};
    model = changeRxnBounds(model,glcO2rxns,0,'l'); %Deactivate normal glucose and O2 uptake
    model = changeRxnBounds(model,'EX_GLC_O2(e)',glclimit/2,'l'); %Activate coupled uptake reaction
    if verbflag == 1
    disp('Constraint set: Glucose/O2 uptake ratio = 2/3')
    end
    
        case 2
            error('Model/Objective/Constraints combination not supported')  
        case 3
            error('Model/Objective/Constraints combination not supported')  
    end


    case 3
         options.constraintsname='O2<11.5';
       switch model_id 
           
           case 1
          %Invoke qO2_max < 11.5 in SCHUETZR:
            model = changeRxnBounds(model,'EX_O2(e)',-11.5,'l');
                     
               
           case 2
               model = changeRxnBounds(model,'EX_O2(e)',-11.5,'l');
           case 3
               model = changeRxnBounds(model,'EX_o2(e)',-11.5,'l');
       
       end
   
        if verbflag == 1
            disp('Constraint set: Max O2 uptake rate = 11.5 mmol/g*h')
        end  
       
    case 4
       options.constraintsname='O2<14.75';
       switch model_id 
           
           case 1
            %Invoke qO2_max < 14.75 in SCHUETZR:
            model = changeRxnBounds(model,'EX_O2(e)',-14.75,'l'); 
           
                 
           case 2
               model = changeRxnBounds(model,'EX_O2(e)',-14.75,'l'); 
           case 3
       end
       
        if verbflag == 1
            disp('Constraint set: Max O2 uptake rate = 14.75 mmol/g*h')
        end
    
    case 5  %Invoke ATP Maintenance:
        options.constraintsname='ATP maint.';
        switch model_id 
            
            case 1
       
            model = changeRxnBounds(model,'maint',7.6,'b');         
            case 2
            model = changeRxnBounds(model,'ATPM',7.6,'b');    
            case 3
             model = changeRxnBounds(model,'ATPM',7.6,'b');
        
        end
       
    if verbflag == 1
        disp('Constraint set: ATP maintenance requirement = 7.6 mmol/g*h')
    end
    case 5.1
        switch model_id 
            
            case 1
       
            model = changeRxnBounds(model,'maint',7.6,'l');         
            case 2
            model = changeRxnBounds(model,'ATPM',7.6,'l');    
            case 3
             model = changeRxnBounds(model,'ATPM',7.6,'l');
        
        end
              
    case 6
        options.constraintsname='Rx. bounds';
        %Set max values for all fluxes in fixed model to 200% of maximum glucose uptake rate
        %(Slightly different than described by Schuetz et al. 

    switch model_id
        
        case 1
        limit = -2*model.lb(2);
        for i = 12:85 %Exclude exchange reactions
        model.ub(i) = limit; 
        if model.lb(i) ~=0
        model.lb(i) = -limit;
        end
        end
            
            
        case 2
            error('Model/Objective/Constraints combination not supported') 
            
        case 3
            error('Model/Objective/Constraints combination not supported') 
            
       
    end
    
     if verbflag == 1
        disp('Constraint set: Bounds on all reactions set to 200 % of max glucose uptake rate')
     end

    
    case 7 %Combination of previous constraints:
   options.constraintsname='All';
        switch model_id
           
           case 1
            %Case 1
            POobjrxns = {'nuo','cydAB'};
            model = changeRxnBounds(model,POobjrxns,0,'b');
        
            %Case 2
            glclimit = model.lb(2);
            glcO2rxns = {'EX_O2(e)','EX_GLC(e)'};
            model = changeRxnBounds(model,glcO2rxns,0,'l'); %Deactivate normal glucose and O2 uptake
            model = changeRxnBounds(model,'EX_GLC_O2(e)',glclimit/2,'l'); %Activate coupled uptake reaction
            %Case 3:
            model = changeRxnBounds(model,'EX_O2(e)',-11.5,'l');
        
            %Case 5:
            model = changeRxnBounds(model,'maint',7.6,'b');
        
            %Case 6:
            limit = -2*model.lb(86);
            for i = 12:85 %Exclude exchange reactions
            model.ub(i) = limit; 
            if model.lb(i) ~=0
            model.lb(i) = -limit;
            end
            end        
               
            case 2
               error('Model/Objective/Constraints combination not supported') 
               
               
           case 3
               error('Model/Objective/Constraints combination not supported') 
               
               
        end
         if verbflag == 1
            disp('Constraint set: All constraints')
         end   
               
            case 8 %Combination of constraints 3,5 and 6. 
                %Constraint 3:
                switch model_id 
           
                case 1
                %Invoke qO2_max < 11.5 in SCHUETZR:
                model = changeRxnBounds(model,'EX_O2(e)',-11.5,'l');
                     
               
                model = changeRxnBounds(model,'maint',7.6,'b'); 
                
                case 2
                model = changeRxnBounds(model,'EX_O2(e)',-11.5,'l');
                model = changeRxnBounds(model,'ATPM',7.6,'b'); 
                    
                case 3
                model = changeRxnBounds(model,'EX_o2(e)',-11.5,'l');
                model = changeRxnBounds(model,'ATPM',7.6,'b');         
                
                
                end
        
        case 0
            options.constraintsname='Default';
            model = rawmodel;
            if verbflag == 1
                disp('Using default model constraints')
            end
end
       
%Define growth media:

switch exp_id
    case 1
        if verbflag == 1
        disp('Simulating aerobic culture. Maximal uptake rates at default values.')
        end
    case 2
           model = changeRxnbounds(model,'EX_o2(e)',0,'b');
          if verbflag == 1
            disp('Simulating anaerobic culture. Maximal oxygen uptake set to zero.')
          end
                
    case 3
         error('Model/Objective/Constraints combination not supported') 
end

           
%Set growth rate to experimentally determined value:


BMreactnum = 25; %This should be changed as to be read in together with experimental data. Hardcoded for now.
if constrainBM == 1
    disp('runsim.m: Constraining growth rate to experimental range.')
    tempopt = options;
    tempopt.fluxconstrainset = BMreactnum;
    model = constrainfluxes(model,tempopt);
end


  
%Set objective

FBArunflag =0;
%Change to 1 when FBA simulation has been run.

%Parameters for optimizeCbModel, and other relevant flags:
options.osensestr='max'; %By default, maximize the objective function.
minnorm=0; %By default, minimization of fluxes is turned off.
options.objectivetype = 'linear'; %Describe what kind of objective is optimized. For use in QPmindist.

switch objective
    case 1
        objectivename ='BM';
        switch model_id
           %Default objective is biomass. Need only ensure that biomass
           %reaction is unconstrained:
            case 1
                model=changeRxnBounds(model,'biomass',0,'l');
                model=changeRxnBounds(model,'biomass',1000,'u');
            case 2
               model=changeRxnBounds(model,'Biomass_Ecoli_core_w_GAM',0,'l');
               model=changeRxnBounds(model,'Biomass_Ecoli_core_w_GAM',1000,'u');
            case 3
                model =changeRxnBounds(model,'Ec_biomass_iJO1366_core_53p95M',0,'l');
                model =changeRxnBounds(model,'Ec_biomass_iJO1366_core_53p95M',1000,'u');
        end
        
        
        if verbflag == 1
            disp('Objective set to max biomass production. Constraints on biomass have been lifted.')
        end
             
        if constrainglucose > 0
            if verbflag == 1
                fprintf(1,'%s \n','runsim.m. Constraining glk flux.');
            end
            tempopt = options;
            tempopt.fluxconstrainset = 1;
            
            model = constrainfluxes(model,tempopt);  
            
            
            if debugmode < 0
                header = {'Reaction #','Reaction name','LB','UB'};
                fprintf(1,'%10s \t %10s \t %10s \t %10s\n',header{1,:});
                for i = 1:length(model.rxns)
                    fprintf(1,'%10d \t %10s \t %10.1f \t %10.1f \n',i,model.rxns{i},model.lb(i),model.ub(i));
                end
            end
            
        end
        
         if constrainglucose > 1
            load('expdata.mat')
            glucoseconsumption = expdata.perrenoud.abs.batch.aerobe.glucoseconsumption;
             if verbflag == 1
                fprintf(1,'%s \n','runsim.m. Constraining glucose uptake flux.');
            end
          
            switch model_id
                case 1
                    rn = 'EX_GLC(e)';
                   
                case 2
                    rn = 'EX_glc(e)';
                   
                case 3
                    rn = 'EX_glc(e)';
            end
            
          model = changeRxnBounds(model,rn,-glucoseconsumption(1)-glucoseconsumption(2),'l');
          model = changeRxnBounds(model,rn,-glucoseconsumption(1)+glucoseconsumption(2),'u');
         end
        
        
        
    case 2
        objectivename='Max ATP';
        switch model_id
            
            case 1
        atpOBJrxns ={'pgk','pykAF','sucCD','atp','ackAB_tdcD_purT'};
        
            case 2
        atpOBJrxns ={'ATPS4r','Pyk'};               
            case 3
               atpOBJrxns ={'AP5AH','ATPS4rpp','PYK'};
        end
            model = changeObjective(model,atpOBJrxns);
            
            if verbflag == 1
                disp('Objective reactions set: ATP producing fluxes')
            end
    case 3.1
        objectivename='Flux minimization (euclidean norm)';
        if verbflag == 1
           disp('Objective is minimization of overall intracellular flux (euclidean norm)')
           disp('Performing custom FBA...')
        end
        
        %setMIQPsolvertogurobi5 = changeCobraSolver('gurobi5','MIQP');
        %setQPsolvertogurobi5 = changeCobraSolver('gurobi5','QP');
        
        minnorm = 1e-6;
        options.objectivetype ='euclidean';
        warning('RUNSIM:UnsupportedObjective','Untested objective. Proceed with caution.')
    case 3.2
        objectivename='Flux minimization (manhattan norm)';
        if verbflag == 1
            disp('Objective is minimization of overall intracellular flux (taxicab norm)')
            disp('Performing custom FBA...')
        end
        minnorm = 'one';
        options.objectivetype ='taxicab';
    case 4
        objectivename='min glucose';
        switch model_id
            case 1
                glcRXNs={'EX_GLC(e)','EX_GLC_O2(e)'};
            case 2
                glcRXNs={'EX_glc(e)'};
            case 3
                glcRXNs={'EX_glc(e)'};
        end
        model = changeObjective(model,glcRXNs);
        options.osensestr='max';
        
        
        %Lift constraints on glucose flux: 
        disp('Lifting constraints on glucose flux.')
        switch model_id
            case 1
                model.lb(2) = rawmodel.lb(2)
                model.ub(2) = rawmodel.ub(2)
            case 2
                model.lb(28) = rawmodel.lb(28)
                model.ub(28) = rawmodel.ub(28)
            case 3
                model.lb(164) = rawmodel.lb(164)
                model.ub(164) = rawmodel.ub(164)
        end
        if verbflag == 1
            disp('Objective set: Minimization of glucose fluxes.')
        end
        
     case 5
         objectivename='Min redox';
         switch model_id
             case 1
            NADHrxns={'aceEF','maeA','sucAB','udhA','fdhF_fdoGHI','ldhA','maeB','zwf','gnd','pntAB','frdABCD_sdhAB','dld'};
            model = changeObjective(model,NADHrxns);
        %'gapA','mdh','icd',adhE_mhpF, adhPC_adhE_r2 and sdhABCD excluded due to reversibility
            options.osensestr = 'min'; %Objective is MINIMIZATION of redox fluxes
        
            %model.c(
        
         if verbflag == 1
             disp('Objective reactions set: Minimization of redox reactions.')
         end
         
         end
    case 6
        objectivename='Min ATP';
         switch model_id
            
            case 1
        atpOBJrxns ={'pgk','pykAF','sucCD','atp','ackAB_tdcD_purT'};
       
            case 2
                atpOBJrxns ={'ATPS4r','Pyk'};
                
            case 3
                 atpOBJrxns ={'AP5AH','ATPS4rpp','PYK'};
        end
          model = changeObjective(model,atpOBJrxns);
        
           options.osensestr = 'min'; 
         if verbflag == 1
             disp('Objective reactions set: Minimization of ATP producing reactions.')
         end
         
end

%Optional, constrain one or some fluxes:

if fluxconstrainset ~=0
   model = constrainfluxes(model,options);
end





options.objectivename=objectivename;
%perform FBA:
if FBArunflag ==0;

    if verbflag == 1
    disp('Performing FBA...') 
    fprintf(1,'%s %s \n','Model:',options.modelname);
    fprintf(1,'%s %s \n','Constraints:',options.constraintsname);
    fprintf('%s %s \n%s %s \n%s %s \n','Objectivetype:',options.objectivetype,'Objective:',options.objectivename,'Sense:',options.osensestr);
    end
    result = optimizeCbModel(model,options.osensestr,minnorm); %Do not attempt to add anything to "result" structure above this point. It will be overwritten here!
    FBArunflag =1;
    result.objectivename = objectivename;
    

if debugmode == 1
    disp('FBArunflag:')
    disp(num2str(FBArunflag))
    disp(result)
end
    
%Send distance from FBA result to output:

tempoptions = options; %Make a copy of the options structure
tempoptions.excludereactions =0; %In calculatig the distance of the optimal solution, do not exclude any reactions


result.constraintsdescription = options.constraintsname;
result.FBAdistance = compdist(result.x,tempoptions,1);
result.FBAsplits = computesplits(result.x,model_id);

end

result.expvalues = refexpfluxes; %Add reference experimental fluxes to output.
    
result.c = model.c; %add model objective vector to output
result.modelname =model.description; %Add model name to the output


%Search for a distance-optimal solution, unless input specified that only FBA should be performed. 
if FBAonly ~=1


    if verbflag == 1
    disp('Starting minimizing/max procedure(s)...')
    end
    
    %Call fmindist to minimize/maximize the distance, using the FBA solution as starting point.:
    if usefmincon == 1
        if verbflag == 1
            disp('Minimizing/maximizing distance of raw fluxes using fmincon... Starting point: FBA result.')
        end
        
        if  usesplits == 1
                if verbflag == 1
                disp('Minimizing/maximizing distance of split ratios using fmincon...Starting point: FBA result.')
                end
                
                
       splitresult = struct;
       tempopt = options;         
       tempopt.optsplits = 1;
       tempopt.compsplits = 1;
       splitresult.splitminsolution = fmindist(result,model,tempopt,1);
       splitresult.splitmaxsolution = fmindist(result,model,tempopt,-1);
       splitresult.expsplitvalues = computesplits(result.expvalues,0); 
       splitresult.splitminvalues = computesplits(splitresult.splitminsolution,model_id);
       splitresult.splitmaxvalues = computesplits(splitresult.splitmaxsolution,model_id);
       tempopt.compsplits = 1;
       splitresult.splitmindistance = compdist(splitresult.splitminsolution,tempopt,1);
       splitresult.splimaxdistance = compdist(splitresult.splitmaxsolution,tempopt,1);
       splitresult.splitrawmindistance = compdist(splitresult.splitminsolution,options,1);
       splitresult.splitrawmaxdistance = compdist(splitresult.splitmaxsolution,options,1);
       result.splitresult = splitresult;
       %Not currently working.
        end
    
    options.optsplits = 0;    
    result.Fmin_minsol = fmindist(result,model,options,1);
    result.Fmin_maxsol = fmindist(result,model,options,-1);
    
    tempoptions = options; %Make a copy of the options structure
    tempoptions.excludereactions =0; %In calculating the distance of the optimal solution, do not exclude any reactions
    
    result.Fmin_mindistance= compdist(result.Fmin_minsol,tempoptions,1);
    result.Fmin_maxdistance = compdist(result.Fmin_maxsol,tempoptions,1);
    result.Fmin_minsol_exp = extractflux(result.Fmin_minsol,options);
    result.Fmin_maxsol_exp = extractflux(result.Fmin_maxsol,options);
    end
    
    %Use fmincon iteratively with a specified number of random starting
    %points
    if iterations ~=0
        if verbflag == 1
            disp('Running fmincon with multiple randoms starting points...')
            disp('Number of iterations:')
            disp(num2str(iterations))
        end
        
        k = iterations;
        randresultmatrix = zeros(size(model.S,2),k);
        randminsolutions = zeros(size(model.S,2),k);
        randmaxsolutions = zeros(size(model.S,2),k);
        randmindistance = zeros(1,k);
        randmaxdistance = zeros(1,k);
        %Make a copy of the model, so the objective can be changed without
        %affecting other operations:
        randmodel = model;
        
        for i = 1:k
            randmodel.c = rand(size(randmodel.c,1),1); %Generate a random objective
        randresult = optimizeCbModel(randmodel); %Perform FBA with the random objective
        randresultmatrix(:,i) = randresult.x; %Record the resulting flux vector
        
        randminsolutions(:,i) = fmindist(randresult,model,options,1); %Use fmincon to minimize the distance
        randmaxsolutions(:,i) = fmindist(randresult,model,options,-1); %Use fmincon to maximize the distance
       
        result.randmindistance(i) = compdist(randminsolutions(:,i),options,1); %Record the minimal distances found
        result.randmaxdistance(i) = compdist(randmaxsolutions(:,i),options,1); %Record the maximal distances found
        
        end
        %Compute further outputs:
        result.meanrandmindistance = mean(result.randmindistance);
        result.meanrandmaxdistance = mean(result.randmaxdistance);
        result.stddevrandmindistance = std(result.randmindistance);
        result.stddevrandmaxdistance = std(result.randmaxdistance);
        [result.bestrandmindistance,bestminindex] = min(result.randmindistance); %Find the best of the minimal distances found
        [result.bestrandmaxdistance,bestmaxindex] = max(result.randmaxdistance); %Find the best of the maximal distances found
        
        result.bestrandminsolution = randminsolutions(:,bestminindex);
        result.bestrandmaxsolution = randmaxsolutions(:,bestmaxindex);
        
        if makeplots == 1
            if verbflag == 1
                disp('Generating plot: Minimal and maximal distance found with fmincon using starting points generated by FBA using random objective functions.')
            end
            
            %t =1:k;
            %plot(t,randmindistance,t,randmaxdistance)
            figure('name','Minimal distances found using random starting points')
            bar(randmindistance)
            xlabel('Objective #')
            ylabel('Distance')
            figure('name','Maximal distances found using random starting points')
            bar(randmaxdistance)
            xlabel('Objective #')
            ylabel('Distance')
        end
    end   
    %END of fmincon iterative routine.
    if usegurobi == 1
        
        if verbflag == 1
         fprintf(1,'%s \n','Minimizing distance of raw fluxes using Gurobi:');
        end
        
       gurobi_output  = QPmindist(model,result,options,1);
       
       if isfield(gurobi_output,'x') == 1
       gurobi_minsol = gurobi_output.x;   
       
        tempoptions = options; %Make a copy of the options structure
        tempoptions.excludereactions =0; %In calculatig the distance of the optimal solution, do not exclude any reactions
       
       gurobi_mindist= compdist(gurobi_minsol,tempoptions,1);
      
       result.gurobi_output = gurobi_output;
       result.gurobi_minsol = gurobi_minsol;
       result.gurobi_mindist = gurobi_mindist;
       result.gurobi_minsol_exp = extractflux(result.gurobi_minsol,options);
        
       else warning('runsim:solver','Gurobi was unable to produce a solution.')
       end
       
       if options.maxwithgurobi ~=0
       gurobi_maxsolution = QPmindist(model,result,options,-1);   
       result.guobi_maxsolution =gurobi_maxsolution;
       result.gurobi_maxdistance = eucdistFCG(gurobi_maxsolution.x,options,1);
       end
        
    end
end
%end of "if FBAonly~=1..."
    
if verbflag == 1
    fprintf(1,'%s \n','Summary:');
    fprintf(1,'%s %.1f \n','Relative optimality required is:',optreq);
    
        if usefmincon == 1
        fprintf('%s %.2f \n','Minimal distance found using fmincon (starting point: FBA result):',result.Fmin_mindistance);
       % disp('Minimzal achievable distance found using fmincon (starting point: FBA result):')
        %disp(result.Fmin_mindistance)
        fprintf('%s %.2f \n','Maximal distance found using fmincon (starting point: FBA result):',result.Fmin_maxdistance);
        %disp('Maximal achievable distance found using fmincon (starting point: FBA result) is:')
        %disp(result.Fmin_maxdistance)
        end
        
        if isfield(result,'gurobi_mindistance') == 1
            disp('Minimzal achievable distance found using Gurobi is:')
            disp(result.gurobi_mindistance)
        end
        
        if isfield(result,'gurobi_maxdistance') == 1
            disp('Maximal achievable distance found using Gurobi is:')
            disp(result.gurobi_maxdistance)
        end
        
        if  iterations ~=0
            disp('Minimzal achievable distance found using fmincon (Random starting points) is:')
            disp(result.bestrandmindistance)
            disp('Maximal achievable distance found using fmincon (Random starting points) is:')
            disp(result.bestrandmaxdistance)
            disp('Number of random objectives used to generate starting points:')
            disp(num2str(iterations))
        end
end
%end of "if verbflag ==1..."
    
%Compute correlation coefficient between computed and experimental fluxes:
if computecorrcoef == 1
        if isfield(result,'Fmin_minsolution') == 1
        result.fmincorrcoef = corrcoef(result.expvalues(:,1),extractflux(result.Fmin_minsol,options));
        result.fmincorrcoef = result.fmincorrcoef(2,1);
        end
        if isfield(result,'gurobi_minsolution') == 1
        result.gurobimincorrcoef = corrcoef(result.expvalues(:,1),extractflux(result.gurobi_minsol,options));
        result.gurobimincorrcoef = result.gurobimincorrcoef(2,1);
        end
        if isfield(result,'bestrandminsolution')==1
        result.bestrandmincorrcoef = corrcoef(result.expvalues(:,1),extractflux(result.bestrandminsolution,options));
        result.bestrandmincorrcoef = result.bestrandmincorrcoef(2,1);
        end
    
        if verbflag == 1
        disp('Correlation coefficients:')
        
            if isfield(result,'fmincorrcoef') == 1
            disp('Fmincon solution (starting point FBA result):')
            disp(num2str(result.fmincorrcoef))
        
            end
            if isfield(result,'bestrandmincorrcoef') == 1
            disp('Fmincon solution (best of solutions generated from random starting points):')
            disp(num2str(result.bestrandmincorrcoef))
            end  
            if isfield(result,'gurobimincorrcoef') == 1
            disp('Gurobi solution:')
            disp(num2str(result.gurobimincorrcoef))
            end
        end
end
%End of "if computecorrcoef..."


switch exp_id
    case 1
           expfluxvalues = expdata.perrenoud.abs.batch.aerobe.fluxvalues;
           expsolution = expfluxvalues(:,1);
           experrors = expfluxvalues(:,2);

    case 2
    case 3
end

%generate plot:
%Plot the best fit solution(s) together with experimental values and
%uncertainties:
%Also include FBA solution?
    
   if makeplots == 1
    if verbflag == 1
       disp('Generating result plot...')
    end
   
           
           figure('name','Feasability analysis')
           title('Experimental and computed reaction rates')
           errorbar(expsolution,experrors,'.k')
           hold on
           legendstring = {'Experimental values'};
           if isfield(result,'Fmin_minsol') == 1
           fluxresult =extractflux(result.Fmin_minsol,options);
           scatter(1:size(expfluxvalues,1),fluxresult);
           legendstring{length(legendstring)+1} = 'Fmincon solution';
           end
           if isfield(result,'gurobi_minsol') == 1
           fluxresult =extractflux(result.gurobi_minsol,options);
           scatter(1:size(expfluxvalues,1),fluxresult);
           legendstring{length(legendstring)+1} = 'Quadratic minimization of distance';
           end
           if options.plotFBA == 1
           fluxresult = extractflux(result.x,options);   
           scatter(1:size(expfluxvalues,1),fluxresult);
           legendstring{length(legendstring)+1} = 'FBA solution';
           end
           xlabel('Experimental reaction #')
           ylabel('Flux value (mmol/g*h)')
           legend(legendstring);
   end
   %end of "if makeplots ==1..." 
   %gtext('test') Use this command to place text on the plot manually.
   %text(1,12,'Optimality requirement') Place text on graph using
   %coordinates.



%Compute objective function values of generated solutions:
if isfield(result,'Fmin_minsol') == 1
    result.Fmin_minsol_objval = model.c'*result.Fmin_minsol;
end

if isfield(result,'gurobi_minsol') == 1
    result.gurobi_minsol_objval =model.c'*result.gurobi_minsol;
end



%generate flux report:
if options.fluxreport ==1
    if isfield(result,'Fmin_minsol') == 1
        
        tempopt = options;
        tempopt.solvername = 'fimincon';
       fluxresult =extractflux(result.Fmin_minsol,tempopt);
        result.fluxresult =fluxresult;
        
        result.fluxreport.fmincon = fluxreport(fluxresult,expfluxvalues,options);
   
    end
    
    if isfield(result,'gurobi_minsol') == 1
        tempopt = options;
        tempopt.solvername = 'gurobi';
        gurobifluxresult = extractflux(result.gurobi_minsol,options);
         result.fluxreport.gurobi = fluxreport(gurobifluxresult,expfluxvalues,tempopt);
    end
end
   
%Compare Fmincon and Gurobi results:
if isfield(result,'Fmin_minsol') == 1 && isfield(result,'gurobi_minsol') == 1
   Fminexp = extractflux(result.Fmin_minsol,options);
   Gurobiexp = extractflux(result.gurobi_minsol,options);
   diff = Fminexp-Gurobiexp;
   totdiff = sum(diff);
   result.compare = [Fminexp Gurobiexp diff];
   %number = size(result.compare,1);
    if options.verbflag == 1    
  
    fprintf(1,'%s \r\n','Comparison of Fmincon and Gurobi solutions:');
    fprintf(1,'%s','Excluded reactions (reaction IDs): ');
    fprintf(1,'%d ',options.excludereactions(1,:));
    fprintf(1,'\r\n');
    
    textheader = {'Reaction #','Fmincon','Gurobi','Difference'};
    fprintf(1, '%10s \t %10s \t %10s \t %10s', textheader{1,:});
    
    for row = 1:size(result.compare,1)
    fprintf(1,'\n %10d \t %10.4f \t %10.4f \t %10.4f',row,result.compare(row,:));
    end
    fprintf(1,'\n');
    
    %fprintf(1,'%s \t %10.4f \ t %10.4f \t %10.4f \r\n','sum',sum(result.compare(:,1)),sum(result.compare(:,2)),sum(result.compare(:,3)));
    fprintf(1,'%s %.4f \r\n','Total difference:',totdiff);
    
    
    if makeplots > 0
    figure('name','Differences in reaction rates predicted by Fmincon and gurobi')    
    bar(result.compare(:,3))
    title('Differences in predicted reaction rates')
    xlabel('Experimental reaction #')
    ylabel('Difference (mmol/g*h)')
    set(gca,'XTick',1:size(result.compare,1))
    set(gca,'YLim',[-3 3])
    end
    
    
    end
    
    
    
end
   
%Select the best (most distance-optimal) solution of those calculated:

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
     
end