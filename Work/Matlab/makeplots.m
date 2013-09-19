%Makeplots.m
%
%Generates a selection of relevant plots for evaluation of objective
%functions
%
%Plot best fit distances as function of optimality requirement: 
% 
 options = struct;
 
 options.objective = 1;
 options.model_id = 1;
 options.exp_id = 1;
 options.usegurobi = 0;
 options.iterations=0;
 options.computecorrcoef=0;
 options.verbflag = 0;
 options.makeplots=0;
 
 m = 20;
 
 k = 1/m;
 
 DS = zeros(1,m);
 req = zeros(1,m);
 
 for i = 1:m+1
     disp('Iteration:')
     disp(num2str(i))
     req(i) = (i-1)/m;
     options.optreq = req(i);
     res = runsim(options);
     DS(i) = res.Fmin_mindistance;
 end
 
 figure('name','Distance-optimality vs. Objective-optimality')
 scatter(req,DS)
 xlabel('min(Z/Zmax)')
 ylabel('Distance (mmol/g*h)')
 
  
 