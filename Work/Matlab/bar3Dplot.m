%bar3D

function output = bar3Dplot(options)

if nargin<1
    options = struct;
end

tic()
if isfield(options,'objectives') == 0
    options.objectives = [1 2 4 6];   
end

if isfield(options,'constraints') == 0
    options.constraints = [0 3 4 5 6 7];
end

if isfield(options,'compareconstraints') == 0
    options.compareconstraints = 1;
end

if isfield(options,'comparemodels') == 0
    options.comparemodels = 0;
end

if isfield(options,'measure') == 0
    options.measure = 'distance';
end

objectives = options.objectives;
constraints = options.constraints;

if isfield(options,'models') == 0
    options.models = 1;
end

if isfield(options,'scatterplot') == 0
    options.scatterplot = 1;
end

models = options.models;
scatterplot = options.scatterplot;


a = length(objectives);
b = length(constraints);
c = length(models);



objectivelabels =  cell({});
constraintslabels =  cell({});
modellabels = cell({});

if isfield(options,'model_id') == 0
options.model_id = 1;
end

if isfield(options,'optreq') ==0
options.optreq =1;
end

if isfield(options,'verbflag') == 0
options.verbflag=1;
end

gurobidists = zeros(a,b);
fmindists = zeros(a,b);
diffs = zeros(a,b);
objvals = zeros(a,b);
fminobjvals = zeros(a,b);
FBAobjvals = zeros(a,b);
optfraq = zeros(a,b);
fminoptfraqs = zeros(a,b);

if options.compareconstraints == 1
    output = cell(a,b);
for i = 1:a
    for j =1:b
    options.objective = objectives(i);
    options.constraints = constraints(j);
    result = runsim(options);
    gurobidists(i,j) = result.gurobi_mindist;
    
    objectivelabels{i} = result.objectivename;
    constraintslabels{j} = result.constraintsdescription;
    objvals(i,j) = result.gurobi_minsol_objval;
    FBAobjvals(i,j) = result.f;
    optfraq(i,j) = objvals(i,j)/result.f;
    output{i,j} = result;
    
    if isfield(result,'Fmin_mindistance') == 1
        fmindists(i,j) = result.Fmin_mindistance;
        diffs(i,j) = gurobidists(i,j)-result.Fmin_mindistance;
        fminobjvals(i,j) = result.Fmin_minsol_objval;
        fminoptfraqs(i,j) = fminobjvals(i,j)/result.f;
    end
    
    end
  
end
figure('name','Minimal distance for objective/constraints combinations')
if strcmp(options.measure,'distance') == 1
bar3(gurobidists);
zlabel('(mmol/g*h)');
elseif strcmp(options.measure,'inversedistance') == 1
bar3(1./gurobidists)
zlabel('(mmol/g*h)^-1');
end
title('Minimal distance for different objective/constraints combinations')
xlabel('Constraints');
ylabel('Objective');
set(gca, 'XTickLabel',constraintslabels);
set(gca, 'YTickLabel',objectivelabels);


if scatterplot == 1
    figure('name','Minimal distance for objective/constraints combinations')
    hold all
    for f = 1:a %For every objective
    values = gurobidists(f,:); %Save the distance values for each constraint in a vector
    plot(1:b,values,':o')
    %plot(1:length(b),values,':o')  %Plot the values.
    xlabel('Constraints');
    ylabel('Distance (mmol/g*h)');
    set(gca, 'XTickLabel',constraintslabels);
    set(gca,'XTick',1:b)
    labels = objectivelabels;
    legend(labels)
    end
    str(1) ={'Model:'};
    modelname = result.model.description;
    str(2) = {modelname};
    ylimits = ylim;
    text(1.2,ylimits(2)*0.95,str)
    
    
end



end

if options.comparemodels == 1
output = cell(a,c);
for i = 1:a
    for j =1:c
    options.objective = objectives(i);
    options.model_id = models(j);
    result = runsim(options);
    gurobidists(i,j) = result.gurobi_mindist;
    
    objectivelabels{i} = result.objectivename;
    modellabels{j} = result.model.description;
    objvals(i,j) = result.gurobi_minsol_objval;
    FBAobjvals(i,j) = result.f;
    optfraq(i,j) = objvals(i,j)/result.f;
    output{i,j} = result;
    
    if isfield(result,'Fmin_mindistance') == 1
        fmindists(i,j) = result.Fmin_mindistance;
        diffs(i,j) = gurobidists(i,j)-result.Fmin_mindistance;
        fminobjvals(i,j) = result.Fmin_minsol_objval;
        fminoptfraqs(i,j) = fminobjvals(i,j)/result.f;
    end
    
    end
            
end
figure('name','Minimal distance for objective/model combinations')
if strcmp(options.measure,'distance') == 1
bar3(gurobidists);
zlabel('(mmol/g*h)');
elseif strcmp(options.measure,'inversedistance') == 1
bar3(1./gurobidists)
zlabel('(mmol/g*h)^-1');
end
title('Minimal distance for different objective/constraints combinations')
xlabel('Model');
ylabel('Objective');
set(gca, 'XTickLabel',constraintslabels);
set(gca, 'YTickLabel',objectivelabels);

if scatterplot == 1
    figure('name','Minimal distance for objective/model combinations')
    hold all
    for f = 1:a %For every objective
    values = gurobidists(f,:); %Save the distance values for each constraint in a vector
    plot(1:b,values,':o')
    %plot(1:length(b),values,':o')  %Plot the values.
    xlabel('Model');
    ylabel('Distance (mmol/g*h)');
    set(gca, 'XTickLabel',modellabels);
    set(gca,'XTick',1:b)
    labels = objectivelabels;
    legend(labels)
    end
    
    
    
end








end
toc()


disp('diffs:')
disp(diffs)

disp('gurobidists:')
disp(gurobidists)

disp('fmindists:')
disp(fmindists)

disp('objvals:')
disp(objvals)

disp('fmin objvals:')
disp(fminobjvals)

disp('FBA objvals:')
disp(FBAobjvals)

disp('Optfraq:')
disp(optfraq)

disp('fminoptfraqs:')
disp(fminoptfraqs)


end