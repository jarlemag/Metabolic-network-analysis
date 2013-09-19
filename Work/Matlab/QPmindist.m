%QPmindistanceFCG.m
%Finds the shortest distance between model and experimental data subject to
%FBA constraints.
%
%The quadratic programming problem solved is on the form:
% c'*x + x'*Q*x + alpha
%
%Syntax: solution = QPmindist(model,options)
%model is a COBRA model strucutre
%result is a COBRA FBA result structure
%"options" is a structure, with the following relevant fields.
%
%model_id: model identifier
%exp_id: experiment identifier
%
%sense: Set to 1 for mininization, -1 for mininization
%optreq: Requirement for optimality, relative to objective function value in FBA result.  
function solution = QPmindist(model,result,options,sense)

solution = struct;
if isfield(options,'verbflag') == 0
   options.verbflag = 0;
end

verbflag = options.verbflag;

if nargin < 4
    sense = 1;
    if verbflag == 1
    disp('QPmindist: Input variable "sense" not defined. Set to 1 (minimization) by default')
    end
end

if isfield(options,'debugmode') == 0
    options.debugmode = 0;
end

if isfield(options,'model_id') == 0
    options.model_id = 1;
end

if isfield(options,'exp_id') ==0
    options.exp_id = 1;
end

if isfield(options,'optreq') == 0
    options.optreq = 0;
end

if (options.optreq == 0)&& (options.model_id == 3);
    options.optreq = 1e-3; %Numerical problem is encountered with optreq = 0 for model iJO1366
end

if isfield(options,'dropquad') == 0
    options.dropquad = 0; %option for ignoring quadratic terms. Off by default.
end

if isfield(options,'droplin') == 0 
    options.droplin = 0; %option for ignoring linear terms. Off by default.
end

if isfield(options,'usetokens') == 0
    options.usetokens = 1;
end

if isfield(options,'osensestr') == 0
    options.osesensetr = 'max';
end

model_id = options.model_id;
exp_id = options.exp_id;
optreq = options.optreq;
debugmode = options.debugmode;
dropquad = options.dropquad;
droplin = options.droplin;
usetokens = options.usetokens;

if isfield(options,'excludereactions') == 1
    excludereactions = options.excludereactions;
    if debugmode == 1
    disp('Excludereactions:')
    disp(excludereactions)
    end
else excludereactions=0;
end

%Specify if optimality requirement is a less than or larger than
%inequality, depending on the sense of the objective:
if strcmp(options.osensestr,'max') == 1
    objectivesense = 1; 
else objectivesense = -1;
end

if isfield(options,'QPsenseoverride') == 1    
    fprintf(1,'%s \n %s %d \n %s %d \n','QPmindist: Objective sense value overriden.','Old value:',objectivesense,'New value:',options.QPsenseoverride);
    objectivesense=options.QPsenseoverride;
end

if debugmode == 1
    disp('QPmindist: objectivesense:')
    disp(objectivesense)
end

%For use if adding support for non-linear objectives:
%if isfield(options,'objectivetype') == 0
%    objectivetype = 'linear';
%else objectivetype=options.objectivetype;
%end


load('expdata.mat')
if verbflag == 1
    disp('QPmindist: Loading experimental data...')
end


load('reactionmaps.mat')
if verbflag == 1
    disp('QPmindist: Loading reaction maps...')
end
%The reaction map maps experimental reactions to their corresponding model
%reactions. The format is as follows
%Column 1: Experimental reaction #
%Remaining Columns: Corresponding model reaction. (0 if
%blank)


n = size(model.S,2);
m = size(model.S,1);

%Make a vector of same length as the complete flux result vector, to
%contain the reference fluxes:
Y = zeros(1,n);
B = zeros(1,n);
%Load experimental fluxes

if verbflag == 1
    disp('QPmindist.m: Reading experimental fluxes...')
end

switch exp_id
    case 1
        expflux = expdata.perrenoud.abs.batch.aerobe.fluxvalues;
    case 2
        error('Experimental data not available')
        
    case 3
        error('Experimental data not available')
end


if verbflag == 1
    disp('QPmindist.m: Reading reaction map...')
end

if usetokens == 1 
        switch model_id
        case 1 
            reactionmap = reactionmaps.Fmap2;
        case 2
            reactionmap = reactionmaps.Cmap2;
        case 3
            reactionmap = reactionmaps.Gmap2;
        end
else
        switch model_id
        case 1
            reactionmap = reactionmaps.Fmap;
        case 2
            reactionmap = reactionmaps.Cmap;
        case 3
            reactionmap = reactionmaps.Gmap;
        end
end

modreactlist = zeros(1,size(expflux,1)); %Make a list of all reactions considered in the calculation.

if debugmode > 0
    disp('Reaction map:')
    mapheader = {'Experimental reaction #','Model reaction #'};
    fprintf(1,'%15s \t %15s \r\n',mapheader{1,:});
    for i = 1:size(reactionmap,1)
    fprintf(1,'%15d \t %15d \r\n',reactionmap(i,:));
    %disp(reactionmap)
    end
    
    %Check that no model reaction appears more than once:
    hist = histc(reactionmap(:,2),sort(reactionmap(:,2)));
    disp('# of times model reactions appearing in reaction map:')
    disp(hist)
    
end


%For all experimental reactions...
for i = 1:size(reactionmap,1)
    %...set the reference flux to the flux in the first reaction specified by
    %the reaction map, or the negative of the experimental flux, if the experimental
    %reaction is defined as the reverse of the model reaction:
    modreactlist(i) = abs(reactionmap(i,2));
    Y(modreactlist(i)) = expflux(i)*sign(reactionmap(i,2));
    B(modreactlist(i)) = 1;
   
end
    %Example: For experimental reaction 2 (i = 2), the corresponding SCHUETZR model
    %reaction is model reaction #33. Thus, reactionmap(2,2) = 33, and Y(33)
    %=expflux(2)
    %Example 2: For ECME model, experimental reaction #11 is the reverse of
    %model reaction #77. Thus, Y(77) = expflux(11)*-1


if debugmode ==1 || debugmode == 2
    disp('Y:')
    disp(Y)
end

alpha = sum(Y.^2);

Q=zeros(n,n);

%Consider fluxes which are covered by the experimental data:
if verbflag == 1
    disp('qpmindist.m: Constructing Q matrix...')
    disp('Model reactions used:')
    disp(modreactlist);
end

for i = 1:length(modreactlist)
            Q(modreactlist(i),modreactlist(i)) = 1;
           
end

%Exclude selected reactions:
if verbflag == 1
    disp('QPmindist.m: Excluding specified experimental reactions (if any)...')
end

if nnz(excludereactions)~=0 %If any reactions to be excluded are specified
    for i = 1:length(excludereactions) %for every element in the specification vector
        expreact = excludereactions(i); %Read the number of the experimental reaction to be excluded
        for j = 2:size(reactionmap,2)  %Read the corresponding reactions in the reaction map.
            if reactionmap(expreact,j)~=0 %For every model reaction associated with that experimental reaction
        modreact = abs(reactionmap(expreact,j)); % Read the number of the model reaction
        Q(modreact,modreact) = 0; % Set the quadratic objective coefficient for that reaction to zero.
         if debugmode > 0
         fprintf(1,'%s %d %s \r\n','Removing model reaction # ',modreact,' from quadratic objective matrix Q');
         end
         
        B(modreact) = 0; %Set the linear objective coefficient for that reactio to zero
          if debugmode > 0
         fprintf(1,'%s %d %s \r\n','Removing model reaction # ',modreact,' from linear objective vector B');
          end           
            end
        end
    end
else
    if verbflag == 1
      disp('No reactions excluded')
    end
end

if debugmode > 0 
   disp('Number of non-zero entries in Q:')
   nonzero = nnz(Q); %Finds number of non-zero entries in Q
   disp(nonzero);
   I = find(Q); %Finds linear indices of non-zero entries in Q
   [rows columns] = ind2sub(size(Q),I);
   indices = [rows columns];
   disp('Indices of non-zero entries in Q:')
   disp(indices)
   
   disp('Number of non-zero entries in B:')
   nonzeroB = nnz(B);
   disp(nonzeroB);
 
   disp('Indices of non-zero entries in B:')
   disp(find(B))
   
   quadrecs = find(diag(Q));
   linrecs = find(B)';
  
   fprintf(1,'%s \r\n','Model reactions considered in linear and quadratic objective terms:');
   fprintf(1,'% \r\n','(Non-zero entries in B and Q');
   header = {'Linear','Quadratic'};
   fprintf(1,'%10s \t %10s \r\n',header{1,:});
   for i = 1:length(quadrecs) 
       fprintf(1,'%10d \t %10d \r\n',linrecs(i),quadrecs(i));
   end   
end

 
if verbflag == 1
    fprintf(1,'%s','Qpmindist.m: Excluding following experimental reactions: ');
    for i = 1:length(options.excludereactions)
        if i < length(options.excludereactions)
        fprintf(1,'%d%s',options.excludereactions(i),',');
        else
        fprintf(1,'%d%s \n',options.excludereactions(i),'.');    
        end
    end

        
C = -2*B.*Y; %Defining linear objective vector      

if debugmode > 0
    num = 1:length(B);
    t = [num' B' Y' C'];
    disp('B, Y, C:')
    disp(t)
end

if dropquad == 1
    Q = zeros(size(Q,1),size(Q,2));
    disp('Dropquad flag set. Dropping quadratic term.');
end

if droplin == 1
    C = zeros(size(C,1),size(C,2));
    disp('Droplin flag set. Dropping linear term.');
end

if verbflag == 1
    disp('QPmindist.m: Creating Gurobi model structure...')
end

gurobimodel.modelname='QPmodel';
gurobimodel.A = model.S;
gurobimodel.obj = C;
gurobimodel.lb = model.lb;
gurobimodel.ub = model.ub;
gurobimodel.objcon = alpha; %Constant offset = alpha
gurobimodel.Q=sparse(Q);


%Sense (minimization/maximization)
if sense == -1
%gurobimodel.modelsense = 'max';
%Maximization does not work currently. Throws an error about Q matrix not
%positive semi-definite. Ignore for now.
end


%FBA steady state constraints:
gurobimodel.sense = '=';   %Sense of linear constraints. Set to =, as S*v = 0.
gurobimodel.rhs = zeros(1,m); %Righ-hand side of S* v =0   

%Objective optimality constraint (quadratic constraints):
gurobimodel.quadcon(1).Qc=sparse(zeros(n,n)); %No quadratic constraints

if objectivesense == 1
gurobimodel.quadcon(1).q = -model.c; %objective function of the model
gurobimodel.quadcon(1).rhs = -result.f*optreq; %Requirement for optimality
elseif objectivesense == -1
gurobimodel.quadcon(1).q = model.c; %objective function of the model
gurobimodel.quadcon(1).rhs = result.f*(1/optreq);    
else error('QPmindist: Objective sense not specified')
end
disp('QPmindist: Objetive sense:')
disp(objectivesense)

if debugmode == 1
    disp('QPmindist: Target objective value:')
    disp(gurobimodel.quadcon(1).rhs)
end

if verbflag == 1
    disp(gurobimodel)
end

if verbflag == 1
    disp('Performing optimization with Gurobi...')
    fprintf(1,'%s %.1f \n','Objective-optimality requirement:',optreq);
    fprintf(1,'%s %.1f \n','Reference maximal FBA objective value',result.f*optreq);
    
end

params = struct;

if verbflag == 0
    params.outputflag = 0;
end

params.BarHomogeneous =1;
params.OptimalityTol = 1e-3; %default is 1e-6. Range 1e-9 - 1e-2.

%Perform optimization:
solution = gurobi(gurobimodel,params);
%Pass additional technical information to the output:
solution.Y = Y;
solution.Q = Q;
solution.alpha = alpha;
solution.modreactlist = modreactlist;
solution.expflux = expflux;
solution.B=B;
solution.C=C;
solution.reactionmap = reactionmap;

if debugmode == 1 || debugmode == 2

    disp('Alpha:')
    disp(alpha)
    if isfield(solution,'objval') == 1
    objectivevalue = solution.objval;
    disp('Gurobi objectivevalue:')
    disp(num2str(objectivevalue))
    disp('Distance check:')
    CHK = compdist(solution.x,[],1);
    disp(num2str(CHK))
    disp('FBA objective value of gurobi solution')
    disp(num2str(model.c'*solution.x))

    end
end

end



