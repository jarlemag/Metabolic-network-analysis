

models = [1 3];
objectives = [1 2 4 6];

a = length(models);
b = length(objectives);


objectivelabels =  cell({});
modellabels = cell({});


gurobidists = zeros(a,b);
objvals = zeros(a,b);
opt = struct;
opt.makeplots = 0;
FBAobjvals = zeros(a,b);
optfraq = zeros(a,b);
for i = 1:a
    
    for j = 1:b
       opt.model_id = models(i);
        opt.objective = objectives(j);
        result = runsim(opt);
        gurobidists(i,j) = result.gurobi_mindist;
        objvals(i,j) = result.gurobi_minsol_objval;
        FBAobjvals(i,j) = result.f;
        optfraq(i,j) = objvals(i,j)/result.f;
        modellabels{i} = result.model.description;
        objectivelabels{j} = result.objectivename;
    end
end

 figure('name','Minimal distance for different objectives')
    hold all
    for f = 1:a %For every model
    values = gurobidists(f,:); %Save the distance values for each constraint in a vector
    plot(1:b,values,':o')
    %plot(1:length(b),values,':o')  %Plot the values.
    xlabel('Objective');
    ylabel('Distance (mmol/g*h)');
    set(gca, 'XTickLabel',objectivelabels);
    set(gca,'XTick',1:b)
    labels = modellabels;
    legend(labels)
    end

%bar3(gurobidists);
