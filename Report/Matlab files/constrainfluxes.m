%constrainreactions.m
%Applies experimentallt determined fluxes as constraints on a model:
%Syntax:
%model = constrainfluxes(model,options)
%model is a COBRA model structure
%options is a structure with the following fields
%model_id Model identifier
%exp_id Experimental data identifier
%tolerance: The distance between the reported experimental flux value and
%the upper/lower bounds, relative to the experimental uncertainty.
%if tolerance = 1, the upper and lower bounds are equal to the reported
%flux plus/minus the reported uncertainty, respectively.
function model = constrainfluxes(model,options)

if nargin < 2
    options = struct;
end

if isfield(options,'model_id') == 0
    options.model_id = 1;
end

if isfield(options,'exp_id') == 0;
    options.exp_id = 1;
end

if isfield(options,'tolerance') == 0
    options.tolerance = 1;
end

if isfield(options,'fluxconstrainset') == 0 
    %Option to constrain only a subset of fluxes. By default, all fluxes
    %are constrained.
    options.fluxconstrainset = 0;
    disp('constrainfluxes.m: No flux constrain set given. Constraining all fluxes.')
end

model_id = options.model_id;
exp_id = options.exp_id;
tolerance = options.tolerance;
fluxconstrainset = options.fluxconstrainset;

load('reactionmaps.mat')
load('expdata')

if isfield(model,'description') == 1
    modelname = model.description;
else modelname = '?';
end

fprintf(1,'%s %d%s \n','contrainfluxes.m:Using reaction map for model with model ID',model_id,'.');
switch model_id
    case 1
         reactionmap = reactionmaps.Fmap2;
    case 2
        reactionmap = reactionmaps.Cmap2;
    case 3
        reactionmap = reactionmaps.Gmap2;
end
        
fprintf(1,'%s %d%s \n','constrainfluxes.m:Using experimental fluxes from experiment',exp_id,'.');
switch exp_id
    case 1
        expflux = expdata.perrenoud.abs.batch.aerobe.fluxvalues;
    case 2
        error('Data not available')
    case 3
        error('Data not available')
end

if fluxconstrainset == 0
    fprintf(1,'%s %d%s \n','Constraining fluxes. Tolerance factor',tolerance,'.'); 
    for i = 1:size(reactionmap,1)
        if nnz(reactionmap(i,:)) < 3
        
        fprintf(1,'%s%d %s%d%s','Constraining experimental reaction #',i,'(Model reaction #',abs(reactionmap(i,2)),').');       
        if sign(reactionmap(i,2)) == 1
        
        model.lb(abs(reactionmap(i,2))) = (expflux(i,1)-expflux(i,2)*tolerance); 
        model.ub(abs(reactionmap(i,2))) = expflux(i,1)+expflux(i,2)*tolerance;
        elseif sign(reactionmap(i,2)) == -1 
                
        model.lb(abs(reactionmap(i,2))) = -(expflux(i,1)+expflux(i,2)*tolerance); 
        model.ub(abs(reactionmap(i,2))) = -(expflux(i,1)-expflux(i,2)*tolerance);
        end
            
        fprintf(1,' %s %.1f %s %.1f \n','LB:',model.lb(abs(reactionmap(i,2))),'UB:', model.ub(abs(reactionmap(i,2))));
        end
    end
else
    %Constrain only a subset of fluxes:
    for i = 1:length(fluxconstrainset)
        j = fluxconstrainset(i);
        fprintf(1,'%s%d %s%d%s','Constraining experimental reaction #',j,'(Model reaction #',abs(reactionmap(j,2)),').'); 
        if sign(reactionmap(j,2)) == 1
        
        model.lb(abs(reactionmap(j,2))) = expflux(j,1)-expflux(j,2)*tolerance; 
        model.ub(abs(reactionmap(j,2))) = expflux(j,1)+expflux(j,2)*tolerance;
        elseif sign(reactionmap(j,2)) == -1
            model.lb(abs(reactionmap(i,2))) = -(expflux(j,1)+expflux(j,2)*tolerance); 
        model.ub(abs(reactionmap(i,2))) = -(expflux(j,1)-expflux(j,2)*tolerance);
        end
        fprintf(1,' %s %.1f %s %.1f \n','LB:',model.lb(abs(reactionmap(j,2))),'UB:', model.ub(abs(reactionmap(j,2))));
    end
end