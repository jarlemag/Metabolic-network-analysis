%fluxreport.m
%Generates a report comparing computed and experimental fluxes
%syntax
%report = fluxreport(fluxresult,expflux)
%
%input: 
%fluxresult: vector of model result fluxes corresponding to experimental fluxes
%expflux: matrix of experimental fluxes and experimental errors in two
%columns.
%column 1: experimental fluxes
%column 2: experimental errors
%
%output:
%Returns a matrix with the following content:
%Column 1: Experimental reaction # (not reaction # in the model!)
%Column 2: Solution flux
%Column 3: Experimental flux (either absolute, or scaled to the computed flux by glucose flux)
%Column 4: Difference between solution and experimental flux
%Column 5: Difference betweens solution and experimental flux divided by
%experimental uncertainty



function report = fluxreport(fluxresult,expvalues,options)

if isfield(options,'usescaledfluxes') == 0
    options.usescaledfluxes =0;
end

%Implement the option to use scaled fluxes later.

if isfield(options,'verbflag') == 0
    options.verbflag =0;
end

if isfield(options,'printfile') == 0
    options.printfile =0;
end

if isfield(options,'modelname') == 0
    modelname = 'unknown';
else modelname=options.modelname;
end

if isfield(options,'objectivename') == 0
    objectivename = 'unknown';
else objectivename=options.objectivename;
end

if isfield(options,'constraintsname') == 0
    constraintsname = 'unknown';
else constraintsname = options.constraintsname;
end

if isfield(options,'filename') == 0
    filename = 'fluxreport.txt';
else filename = options.filename;
end

if isfield(options,'solvername') == 0
    options.solvername ='?';
end
    

expfluxes = expvalues(:,1);
experrors = expvalues(:,2);

difference = expfluxes(:,1) - fluxresult;
uncertfraction = difference./experrors;
number =(1:size(fluxresult,1))';

%if usescaledfluxes ==1
%scalingfactor = expflux(1)/fluxresult(1);
%difference =  (expflux(:,1)/scalingfactor)-fluxresult;
%uncertfraction = difference./(expflux(:,2)/scalingfactor);
%number =(1:size(fluxresult,1))';
%end

%report = [number fluxresult expflux/scalingfactor difference  uncertfraction];
rawreport = [number fluxresult expfluxes difference experrors uncertfraction];

textheader = {'Ex rxn #','Comp. flux','Exp. flux', 'Diff.', 'Exp. uncert.', 'Diff./uncert.'}; 

report = {textheader; rawreport};

timedate = fix(clock);

outsidefluxes =0;
for i = 1:size(uncertfraction)
    if (abs(uncertfraction(i)) > 1) &&  (abs(uncertfraction(i))~= inf)
        outsidefluxes = outsidefluxes +1;
    end
end

if options.verbflag == 1
    
fprintf(1,'%s %s \r\n','Flux report generated at:',datestr(timedate));
fprintf(1,'%s %s \r\n','Solver:',options.solvername);
fprintf(1,'%s %s \r','Model:',modelname);
fprintf(1,'%s %s \r','Objective:',objectivename);
fprintf(1,'%s %s \r\n','Constraints:',constraintsname);
fprintf(1,'%s %s \r\n','Number of computed fluxes outside experimental uncertainty bounds:',num2str(outsidefluxes));
fprintf(1, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s', textheader{1,:});
for row = 1:size(rawreport,1)
fprintf(1,'\n %10d \t %10f \t %10f \t %10f \t %10f \t %10f',rawreport(row,:));
end
fprintf(1,'\n');
end

if options.printfile == 1
    fileID = fopen(filename,'a');
    
    fprintf(fileID,'%s %s \r\n','Flux report generated at:',datestr(timedate));
    fprintf(fileID,'%s %s \r','Model:',modelname);
    fprintf(fileID,'%s %s \r','Objective:',objectivename);
    fprintf(fileID,'%s %s \r\n','Constraints:',constraintsname);
    fprintf(fileID,'%s %s \r\n','Number of computed fluxes outside experimental uncertainty bounds:',num2str(outsidefluxes));
    fprintf(fileID, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s', textheader{1,:});
    for row = 1:size(rawreport,1)
    fprintf(fileID,'\r\n %10d \t %10f \t %10f \t %10f \t %10f \t %10f',rawreport(row,:));
    end
    fprintf(fileID,'\r\n');
    %fclosestatus = fclose(fileID);
end


end


