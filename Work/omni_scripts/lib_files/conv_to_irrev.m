function [S_ir,c_ir,rxns_ir,lb_ir,ub_ir,fb_match,sel_int_ir,rev2irrev,RM] = conv_to_irrev(S,c,rxns,lb,ub,rev,sel_int)
% CONV_TO_IRREV Convert reversible rxns to irreversible
% 
% [S_ir,c_ir,rxns_ir,lb_ir,ub_ir,fb_match,sel_int_ir,rev2irrev,RM] = conv_to_irrev(S,c,rxns,lb,ub,rev,sel_int)
%
% S         Stoichiometric matrix   
% c         Objective
% rxns      Rxn list
% lb        Lower bounds
% ub        Upper bounds
% rev       Reversible (1/0)
% sel_int   Selectable reactions
%
% S_ir      S-matrix (irrev)
% c_ir      Objective (irrev)
% rxns_ir   Rxn list (irrev)
% lb_ir     Lower bounds (irrev)
% ub_ir     Upper bound (irrev)
% fb_match  Forward/backward matches
% sel_int_ir    Selectable reaction (irrev)
% rev2irrev Mapping from reversible to irreversible rxns
% 
% Markus Herrgard 10/14/03

[m,n] = size(S);
% Find reversible reactions
rev_ids = find(rev == 1);
% Number of reversible
nr = length(rev_ids);
% Find irreversible rxns
irrev_ids = find(rev == 0);
% Number of irrevs
nir = length(irrev_ids);

% Mapping
rev2irrev = [];
RM = zeros(nir+2*nr,n);
fb_match = [];

% Irreversible rxns
S_ir = S(:,irrev_ids);
rxns_ir = rxns(irrev_ids);
lb_ir = lb(irrev_ids);
ub_ir = ub(irrev_ids);
c_ir = c(irrev_ids);
if (nargin > 6) 
    sel_int_ir = sel_int(irrev_ids);
else
    sel_int_ir = [];
end

% Loop through irreversible rxns
for i = 1:nir
    rev2irrev{irrev_ids(i)}.conv = i;
    RM(i,irrev_ids(i)) = 1;
end

% Loop through reversible rxns
for i = 1:nr
    rxn_id = rev_ids(i);
    rev2irrev{rxn_id}.conv = [nir+2*i-1 nir+2*i];
    RM(nir+2*i-1,rxn_id) = 1;
    RM(nir+2*i,rxn_id) = 1;
    % S matrix
    s_col = S(:,rxn_id);
    S_ir = [S_ir s_col -s_col];
    % Reaction list
    rxns_ir{nir+2*i-1} = [rxns{rxn_id} '_f'];
    rxns_ir{nir+2*i} = [rxns{rxn_id} '_b'];
    % Lower/upper bounds
    lb_val_b = 0;
    lb_val_f = 0;
    ub_val_b = 0;
    ub_val_f = 0;
    if (lb(rxn_id) <= 0)
        % Negative lower bound -> Upper bound for reverse direction
        ub_val_b = -lb(rxn_id);
    else
        % Positive lower bound -> Lower bound for forward direction
        lb_val_f = lb(rxn_id);
    end
    if (ub(rxn_id) >= 0)
        % Positive ub -> ub for forward dirn
        ub_val_f = ub(rxn_id);
    else
        % Negative ub -> lb for rev dirn
        lb_val_b = -ub(rxn_id);
    end
    % Assemble bounds
    lb_ir = [lb_ir;lb_val_f;lb_val_b];
    ub_ir = [ub_ir;ub_val_f;ub_val_b];
    % Forward/backward match
    if (nargin > 6)
        sel_int_ir = [sel_int_ir; sel_int(rxn_id); sel_int(rxn_id)];
    end
    fb_match(i,1) = nir+2*i-1;
    fb_match(i,2) = nir+2*i;
    % Objective function
    if (c(rxn_id) >= 0)
        c_ir = [c_ir; c(rxn_id); 0];
    else
        c_ir = [c_ir; 0; c(rxn_id)];
    end
end