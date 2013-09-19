           figure('name','Feasability analysis')
           title('Experimental and computed reaction rates')
           errorbar(expsolution,experrors,'.k')
           hold on
           legendstring = {'Experimental values'};
           if isfield(result,'Fmin_minsol') == 1
           fluxresult =extractflux(result.Fmin_minsol,options);
           scatter(1:size(expfluxvalues,1),fluxresult);
           legendstring{length(legendstring)+1} = 'Fmincon solution';
           end
           if isfield(result,'gurobi_minsol') == 1
           fluxresult =extractflux(result.gurobi_minsol,options);
           scatter(1:size(expfluxvalues,1),fluxresult);
           legendstring{length(legendstring)+1} = 'Quadratic minimization result';
           end
           if options.plotFBA == 1
           fluxresult = extractflux(result.x,options);   
           scatter(1:size(expfluxvalues,1),fluxresult);
           legendstring{length(legendstring)+1} = 'FBA solution';
           end
           xlabel('Experimental reaction #')
           ylabel('Flux value (mmol/g*h)')
           legend(legendstring);