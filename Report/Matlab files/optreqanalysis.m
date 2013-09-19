%optreq.m:
%
%Analyze how the minimal achievable distance as determined through runsim varies with requirement of
%optimality.
%
%output = optreq(options)
%options is a structure with the following fields:
%models: a vector with one or more model IDs
%objectives: a vector with one or more objective IDs
%experiments: a vector with one or more experiment IDs
%stepsize: a number between zero and one


function output = optreqanalysis(options)


 output = figure('name','Distance-optimality vs. Objective-optimality');
 xlabel('min(Z/Zmax)')
 ylabel('Distance (mmol/g*h)')
 %title('Distance-optimality vs. Objective-optimality')
 hold all
 if nargin <1
    options = struct();
    options.verbflag = 1;
end
 
 if isfield(options,'models') == 0
     options.models = [1 2 3];
 end
 
 if isfield(options,'objectives') == 0
     options.objectives = [1 2];
 end
 
 if isfield(options,'experiments') == 0
     options.experiments = 1;
 end
 
 if isfield(options,'steps') == 0
     options.steps = 10;
 end
 
 if isfield(options,'usegurobi') == 0
     options.usegurobi = 1;
 end
 
 if isfield(options,'usefmincon') == 0
     options.usefmincon = 0;
 end
 
 if isfield(options,'verbflag') == 0
     options.verbflag = 1;
 end
 
 if isfield(options,'makeplots') ==0
     options.makeplots =0;
 end
 
 models = options.models;
 objectives = options.objectives;
 experiments = options.experiments;
 steps = options.steps;
 

 %Determine the number of graphs to be generated
 graphnumber = length(models)*length(objectives)*length(experiments);
 
 labels = cell(graphnumber,1);


 DS = zeros(1,steps);
 req = zeros(1,steps);

 currentcombination=0;
for n = 1:length(models) 
    
    for l = 1:length(objectives)
        
        for o = 1:length(experiments)
        
        options.model_id = models(n);
        options.objective = objectives(l);
        options.exp_id = experiments(o);
        
            for i = 1:steps+1
            fprintf(1,'%s %d \r\n','Iteration:',i);
            %disp('Iteration:')
            %disp(num2str(i))
            req(i) = (i-1)/steps;
            options.optreq = req(i);
            
            if (options.optreq == 0) && (options.model_id == 3); 
                %Gurobi does not complete if optimality requirement is exactly zero when using iJO1366.
                options.optreq = 0.01;
            end
            
            if options.optreq == 1
                options.optreq = 0.99;
            end
            
            res = runsim(options);
             
                if options.usefmincon == 1
            
                DS(i) = res.Fmin_mindistance;
                elseif options.usegurobi == 1
                DS(i) = res.gurobi_mindist;   
                end
            end
        currentcombination=currentcombination+1;
            
            labels{currentcombination} = [res.modelname,' ',res.objectivename];
        %scatter(req,DS)
        plot(req,DS,':o')
         end
       
     end
    
end

legend(labels) %Apply labels to the plot

end     