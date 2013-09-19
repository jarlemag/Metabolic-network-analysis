%singleparametertest.m
%
%Use a single-parameter test to check if the model/objective/constraints
%are consistent with the experimental data. Reaction rates are constrained
%to experimentally determined fluxes. 
%
%syntax:
%output = singleparamtest(options,m,stepsize)
%options: a structure containing the following fields:
%model_id 
%exp_id
%stepsize: How much the tolerance should be varied in each step. Smaller
%stepsize gives more iterations.
%m
%To test while constraining only some reactions, include a field
%fluxconstrainset in the options structure, containing a vector of the
%experimental reaction IDs of the reactions to be constrained.

function  output = singleparamtest(options)


if isfield(options,'debugmode') ==0
    options.debugmode =0;
end

if isfield(options,'maxtol') == 0
    options.maxtol = 10;
end

if isfield(options,'stepsize') == 0
    options.stepsize = 1;
end

if isfield(options,'model_id') == 0
    options.model_id = 1;
end

if isfield(options,'exp_id') == 0
    options.exp_id = 1;
end


debugmode = options.debugmode;
model_id = options.model_id;
maxtol = options.maxtol;
stepsize = options.stepsize;

load('modeldata.mat')

switch model_id
    case 1
        model = modeldata.SCHUETZR;
    case 2
        model = modeldata.ECME;
    case 3
        model = modeldata.iJO1366b;
end

modelname = model.description;
%For varying flux deviation tolerance, optimize the model with respect to
%the FBA objective:

k = 0;
%results = zeros(1,maxtol);
for i = 1:stepsize:maxtol
    k = k+1;
    tolerance(k) =i; 
    options.tolerance = tolerance(k);
disp('Constraining model...')
model = constrainfluxes(model,options);
disp('iteration:')
disp(num2str(k))
result = optimizeCbModel(model);
results(k) = result.f;
end

if debugmode > 0
    disp('Tolerance:')
    disp(tolerance)
    disp('Objective values:')
    disp(results)
end

figure('Name','Single parameter analysis')
scatter(tolerance,results)
title('Single parameter analysis')
ylabel('Objective value')
xlabel('Maximal flux deviation')
output = [tolerance; results];

str(1) ={'Model:'};
str(2) = {modelname};
ylimits = ylim;
text(1.2,ylimits(2)*0.95,str)

performFVA = 0;
if max(results) > 1e-2
    performFVA = 1;
end

if performFVA == 1
%Perform FVA to find which reactions are limiting:
optPercentage = 10;
osenseStr = 'max';
[minFlux,maxFlux] = fluxVariability(model,optPercentage,osenseStr);


range = maxFlux - minFlux;
values_diff = [minFlux range];
figure('name','Flux Variability Analysis (all reactions)')
bh = bar(values_diff,'stacked');
set(bh(1),'FaceColor','none','EdgeColor','none')
title('Flux Variability Analysis (all reactions)')
ylabel('Flux (mmol/g*h)')
xlabel('Model reaction #')
set(gca,'YLim',[-100 100])

minFluxexp = extractflux(minFlux,options);
maxFluxexp = extractflux(maxFlux,options);

exprange = maxFluxexp -minFluxexp;
expvalues_diff = [minFluxexp exprange];
figure('name','Flux Variability Analysis (experimental reactions)')
bh = bar(expvalues_diff,'stacked');
set(bh(1),'FaceColor','none','EdgeColor','none')
title('Flux Variability Analysis (experimental reactions)')
ylabel('Flux (mmol/g*h)')
xlabel('Experimental reaction #')

str(1) ={'Optimality requirement:'};
str(2) = {optPercentage/100};
ylimits = ylim;
text(1.2,ylimits(2)*0.90,str)

end

end
%Attempt to fit the model growth rate to the experimental growth rate:
%solution = fmincon(@(x) paramfit(x,options,sense),start,-model.c',-optreq*result.f,full(model.S),nullvector,model.lb,model.ub,[],fopt);