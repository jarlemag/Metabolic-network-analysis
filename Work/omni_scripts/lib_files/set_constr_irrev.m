function [sel_ir,b_ir,csense_ir] = set_constr_irrev(sel_rxn,rxns,rxns_ir,b,rev2irrev,csense)
%SET_CONSTR_IRREV - Set constraints for a subset of rxns (converts reversible to irreversible)
%
% [sel_ir,b_ir,csense_ir] = set_constr_irrev(sel_rxn,rxns,rxns_ir,b,rev2irrev,csense)
% 
% sel_rxn   Reaction selection cell array (for rev representation)
% rxns      Rxn names in the order of columns of the S matrix
% rxns_ir   Irreversible reaction names
% b         RHS vector ordered as sel_rxn
% rev2irrev Mapping from reversible to irreversible
% csense    Constraint senses ordered as sel_rxn (opt)
%
% sel_ir    Selection index (e.g. [2 4 5 9 10])
% b_ir      Correctly ordered rhs
% csense_ir Correctly ordered 
%
% 6/9/05 Changed this so that it allows multiple occurences of the same
% rxn
%
% Markus Herrgard 10/14/03

sel_ir = [];
b_ir = [];
csense_ir = '';
% Bit vector describing the processing status of this rxn
proc_rxns = zeros(length(sel_rxn),1);
for i = 1:length(sel_rxn)
    
    % Has this rxn already been processed?
    if (~proc_rxns(i))
        
        % Get reaction name
        sr = sel_rxn{i};
        % Find out if there are more than one constraint for this rxn
        match_ids = find(strcmp(sel_rxn,sr));
        % How many matching rxns
        n_rxn = length(match_ids);
        % Mark as processed
        proc_rxns(match_ids) = 1;
        
        % Get constraint values
        bvals = b(match_ids);
        % Get constraint directions
        senses = csense(match_ids);
        
        % Find reaction index in reversible reaction set
        rxn_ind = find(strcmp(rxns,sr));
        
        if (~isempty(rxn_ind))
            % Find reaction index in irreversible reaction set
            map = rev2irrev{rxn_ind}.conv;
            
            % Add reaction indices and constraint values into the list
            if (length(map) == 1) % Irreversible
                % Value of map is directly the index of the reaction
                for j = 1:n_rxn
                    sel_ir = [sel_ir; map];
                    b_ir = [b_ir;bvals(j)];
                    csense_ir(end+1) = senses(j);
                end
            else % Reversible
                % Only one constraint
                if (n_rxn == 1 | (bvals(1) == bvals(2)))
                    if (n_rxn == 2)
                        senses = 'E';
                    end
                    % Map contains both forward and reverse reaction indices
                    sel_ir = [sel_ir; map(1); map(2)];
                    if (bvals(1) >= 0)
                        b_ir = [b_ir;bvals(1);0];
                        % Deal with different directions of constraints
                        if (strcmp(senses,'E'))
                            csense_ir(end+1) = 'E';    
                            csense_ir(end+1) = 'E';
                        elseif (strcmp(senses,'G'))
                            csense_ir(end+1) = 'G';
                            csense_ir(end+1) = 'E';
                        elseif (strcmp(senses,'L'))
                            csense_ir(end+1) = 'L';
                            csense_ir(end+1) = 'G';
                        end 
                    else
                        b_ir = [b_ir;0;-bvals(1)];
                        % Deal with different directions of constraints
                        if (strcmp(senses,'E'))
                            csense_ir(end+1) = 'E';    
                            csense_ir(end+1) = 'E';
                        elseif (strcmp(senses,'G'))
                            csense_ir(end+1) = 'G';
                            csense_ir(end+1) = 'L';
                        elseif (strcmp(senses,'L'))
                            csense_ir(end+1) = 'E';
                            csense_ir(end+1) = 'G';
                        end 
                    end
                else % More than one constraint (the only case that makes sense is a <= v <= b)
                    low_c = min(bvals);        
                    hi_c = max(bvals);
                    
                    if ((low_c > 0) & (hi_c > 0))  % Both positive
                        sel_ir = [sel_ir; map(1); map(1); map(2)];
                        b_ir = [b_ir;low_c;hi_c;0];
                        csense_ir(end+1:end+3) = 'GLE';
                    elseif ((low_c < 0) & (hi_c < 0)) % Both negative
                        sel_ir = [sel_ir; map(1); map(2); map(2)];
                        b_ir = [b_ir;0;-hi_c;-low_c];
                        csense_ir(end+1:end+3) = 'EGL';
                    else % Low positive, hi negative
                        sel_ir = [sel_ir; map(1); map(2)];
                        b_ir = [b_ir;hi_c;-low_c];
                        csense_ir(end+1:end+2) = 'LL';
                    end
                end
            end
        end
    end
end
