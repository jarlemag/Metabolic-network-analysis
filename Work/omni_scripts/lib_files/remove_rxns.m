function model_out = remove_rxns(model,rxn_remove_list,nrev_fl,met_fl)
%REMOVE_RXNS Remove reactions from a model
%
% model_out = remove_rxns(model,rxn_remove_list,nrev_fl,met_fl)
%
% model             Input model structure
% rxn_remove_list   Cell array of reaction names to be removed
% nrev_fl           Irreverseble (true) or reversible (false) reaction
%                   format
% met_fl            Remove unused metabolites too (true)
%
% Markus Herrgard 7/22/05

if (nargin < 3)
  nrev_fl = false;
end
if (nargin < 4)
  met_fl = true;
end

[n_mets,n_rxns] = size(model.S);

% Find indices to rxns in the model
[is_valid_rxn,remove_index] = ismember(rxn_remove_list,model.rxns);
remove_index = remove_index(is_valid_rxn);

% Remove reversible tag from the reverse reaction if the reaction to be
% deleted is reversible
if (nrev_fl)
  for i = 1:length(remove_index)
    rem_id = remove_index(i);
    if (model.match(rem_id) > 0)
      rev_rxn_id = model.match(rem_id);
      model.rev(rev_rxn_id) = 0;
      model.rxns{rev_rxn_id} = model.rxns{rev_rxn_id}(1:end-2);
    end
  end
end

% Construct vector to select rxns to be included in the model rapidly
select_rxns = (ones(n_rxns,1) == 1);
select_rxns(remove_index) = false;

% Construct new model
model_out.S = model.S(:,select_rxns);
model_out.rxns = model.rxns(select_rxns);
model_out.lb = model.lb(select_rxns);
model_out.ub = model.ub(select_rxns);
model_out.rev = model.rev(select_rxns);
if (isfield(model,'c'))
    model_out.c = model.c(select_rxns);
end
if (isfield(model,'int'))
    model_out.int = model.int(select_rxns);
end
% Reconstruct the match list
if (nrev_fl)
  model_out.match = reassign_fb_match(model.match,select_rxns);
end

% Remove metabolites that are not used anymore
if (met_fl)
  sel_mets = any(model_out.S ~= 0,2);
  model_out.S = model_out.S(sel_mets,:);
  model_out.mets = model.mets(sel_mets);
else
  model_out.mets = model.mets;
end
