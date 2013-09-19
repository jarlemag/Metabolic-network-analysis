%fmindistanceFCG.m
%
%uses the fmincon function to  seek a solution of a given FBA problem which
%minimizes (or maximizes) the distance between a set of fluxes within the solution and a
%set of experimental fluxes
%
%syntax: distance = fmindist(result,model,options,sense)
%
%Input:
%result:
%model: the COBRA model structure
%
%options: a structure with the following fields 
%model_id: a identifier for the model
%exp_id: Identifier for the experimental data
%
%sense: Set to 1 for  maximization, -1 for minimization
%Output:
%
%Solution: Allowed flux distribution (of all fluxes) giving minimal (or maximal) distance between the computed and
%experimental flux vectors.


function solution = fmindist(result,model,options,sense)

%solution = struct;

if isfield(options,'verbflag') == 0
    options.verbflag = 0;
end

verbflag = options.verbflag;

if nargin < 4
    sense = 1;
    if verbflag == 1
    disp('fmindist: Input variable "sense" not defined. Set to 1 (minimization) by default')
    end
end


if isfield(options,'debugmode') == 0
    options.debugmode = 0;
end

if isfield(options,'model_id') == 0
    options.model_id = 1;
end

if isfield(options,'exp_id') == 0
    options.exp_id = 1;
end

if isfield(options,'optreq') == 0
    options.optreq = 0;
end

if isfield(options,'optsplits') == 0
    options.optsplits = 0;
    %Option to optimize split ratios instead of raw fluxes
end
options.compsplits = options.optsplits;

if isfield(options,'osensestr') == 0
    options.osesensetr = 'max';
end


optreq = options.optreq;
start = result.x; %Use the FBA result as starting point for the fmincon algorithm
nullvector = zeros(size(model.S,1),1);
optsplits = options.optsplits;

%Minimize/mazimize distance to experimental values

fopt = optimset; %initialize fmincon options
if options.debugmode ~=1
    fopt.Display='off'; %Hide fmincon process output if debug mode is off
    
end

fopt.Algorithm ='active-set'; %Use active-set algorithm in fmincon

switch sense 
    case 1
        sensename ='Minimizing';
    case -1
        sensename = 'Maximizing';
end

if verbflag == 1
        switch optsplits
            case 0
            fprintf(1,'%s %s %s \n','Using fmincon with active-set algorithm.',sensename,'distance between raw flux vectors') ;
            case 1
            fprintf(1,'%s %s %s \n','Using fmincon with active-set algorithm.',sensename,'distance between split ratio vectors'); 
        end
            %disp('Using fmincon with active-set algorithm (minimization of distance)...')
end

%Perform optimization with fmincon:
%fmincon syntax: X = fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB)
solution = fmincon(@(x) compdist(x,options,sense),start,-model.c',-optreq*result.f,full(model.S),nullvector,model.lb,model.ub,[],fopt);
    

end