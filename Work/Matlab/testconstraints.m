%Testconstraints.m:

constraints = [1 3 4 5 6 7];

options.model_id = 1;
options.optreq =1;
options.verbflag=1;

results ={};
distances = zeros(1,length(constraints));

for i = 1:length(constraints)
    options.constraints = constraints(i);
    result = runsim(options);
    distances(i) = result.gurobi_mindist;
end
