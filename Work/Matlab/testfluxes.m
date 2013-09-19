
function output = testfluxes(options)


load('reactionmaps.mat')


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


if isfield(options,'model_id') == 0
options.model_id = 1;
end

if isfield(options,'tolerance') == 0
    options.tolerance = 1;
end

if isfield(options,'consecutive') == 0
    options.consecutive = 1;
end

model_id = options.model_id;

switch model_id 
    case 1
        model = readCbModel('SCHUETZR');
        reactionmap = reactionmaps.Fmap2;
    case 2
        model = readCbModel('ECME');
        reactionmap = reactiomaps.Cmap2;
    case 3
        model = readCbModel('iJO1366b');
        reactionmap = reactiomaps.Gmap2;
 end

rawmodel = model;
k = size(reactionmap,1);
results = zeros(k,1);
for i = 1:k
    options.fluxconstrainset = i;
   switch options.consecutive
       case 1
    model = constrainfluxes(rawmodel,options);
       case 0
    model = constrainfluxes(model,options);
   end
    result = optimizeCbModel(model);
    results(i) = result.f;

    
end
figure('name','fluxtest')
bar(results)

output = 1;
end

