%Testobjectives.m:

objectives = [1 2 4 6];

%Max BM, %Max ATP, %


options.model_id = 3;
options.optreq =1;
options.verbflag=1;
options.constrainglucose = 0;

results ={};
distances = zeros(1,length(objectives));

for i = 1:length(objectives)
    options.objective = objectives(i);
    result = runsim(options);
    distances(i) = result.gurobi_mindist;
end

