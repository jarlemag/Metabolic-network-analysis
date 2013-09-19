function omni_sol = run_omni(model,opt,verb_fl)
%RUN_OMNI Run Optimal Metabolic Network Identification
%
% omni_sol = run_omni(model,opt)
%
% model     Structure containing all necessary variables to described a
%           stoichiometric model and allowed integer variables
%   rxns    Rxns in the model
%   mets    Metabolites in the model
%   S       Stoichiometric matrix (sparse)
%   b       RHS of Sv = b (usually zeros)
%   c       Objective coefficients
%   lb      Lower bounds for fluxes
%   ub      Upper bounds for fluxes
%   rev     Reversibility of fluxes
%   int     Rxns that are allowed to be deleted from the model
%
% opt       OMNI options
%
%   - General options -
%   K            Number of deletions allowed
%   ksense       Direction of # of deletions constraint (G/E/L)
%   H            Max flux through any rxn
%   solve_omni   Solve OMNI problem in GAMS
%   create_gams  Create GAMS input file
%   gams_file    GAMS input file name (optional)
%   sel_int_prev Previously selected integer variables (optional)
%
%   - Constraints for the OMNI problem (set E in paper) - 
%   sel_icc_rxn  Names for rxns with equality/minimum constraints
%   b_icc_rxn    Values for exns with eq/min constraints
%   csense_icc   Constraint senses for eq/min constraints (G/E/L)
%
%   - Measured target fluxes to be included in the objective (set M in
%   paper) - 
%   sel_meas_rxn Names for reactions to be included in the objective 
%   b_meas_rxn   Target values for fluxes in the objective (i.e. measured
%   fluxes)
%   wt_meas_rxn  Weights for fluxes in the objective
%
% verb_fl   Verbose mode
%
% Markus Herrgard 2/7/06

warning off MATLAB:divideByZero;

if (nargin < 3)
    verb_fl = 0;
end

% Unpack model
[rxns,mets,S,b,c,lb,ub,rev,sel_int] = unpack_fba_model(model);

% Set options
K = opt.K;
ksense = opt.ksense;
H = opt.H;
sel_icc_rxn = opt.sel_icc_rxn;
b_icc_rxn = opt.b_icc_rxn;
csense_icc = opt.csense_icc;
solve_omni = opt.solve_omni;
create_gams = opt.create_gams;
if (isfield(opt,'gams_file'))
    gams_file = opt.gams_file; 
end
if (isfield(opt,'sel_int_prev'))
    sel_int_prev = opt.sel_int_prev; 
else
    sel_int_prev = [];
end

%******************************
% Convert to irreversible form
%******************************

% Convert to irreversible reactions
[S_ir,c_ir,rxns_ir,lb_ir,ub_ir,fb_match,sel_int_ir,rev2irrev] = conv_to_irrev(S,c,rxns,lb,ub,rev,sel_int);
if (~isempty(sel_int_prev))
    [tmp,n_exko] = size(sel_int_prev);
    for i = 1:n_exko
        [S_ir,c_ir,rxns_ir,lb_ir,ub_ir,fb_match,sel_int_prev_tmp,rev2irrev] = conv_to_irrev(S,c,rxns,lb,ub,rev,sel_int_prev(:,i));
        sel_int_prev_ir(:,i) = sel_int_prev_tmp;
    end
else
    sel_int_prev_ir = [];
end
[m,n] = size(S_ir);

% Create matchings for reversible reactions in the set selected for KOs
sel_int_f = find(sel_int_ir == 1);
cnt = 0;
for i = 1:length(sel_int_f)
    fwr = sel_int_f(i);
    if (sum(fb_match(:,1) == fwr) > 0)
        cnt = cnt + 1;
        rev_int(cnt,1) = i;
        rev_int(cnt,2) = i+1;
    end
end

%******************************************
% Set constraints and measured flux values
%******************************************

% Set constraints
[sel_ic,b_ic,csense_ic] = set_constr_irrev(sel_icc_rxn,rxns,rxns_ir,b_icc_rxn,rev2irrev,csense_icc);

% Set measured reaction in objective
sel_meas_rxn = opt.sel_meas_rxn';
b_meas_rxn = opt.b_meas_rxn';
wt_meas_rxn = opt.wt_meas_rxn';
n_m = length(sel_meas_rxn);

% Create selection vector in the decoupled representation
% This is to ensure that the objective function for measured reversible
% reactions is constructed correctly
sel_m = zeros(n,1);
ord_ir = [];
b_meas_tmp = [];
wt_meas_tmp = [];
for i = 1:n_m
    rxn_name = sel_meas_rxn{i};
    rxn_id = find(strcmp(rxns,rxn_name));
    if (~isempty(rxn_id)) % Protect against measured fluxes that are not part of the model
        b_meas_tmp = [b_meas_tmp;b_meas_rxn(i)];
        wt_meas_tmp = [wt_meas_tmp;wt_meas_rxn(i)];
        % Reversible rxns
        if (rev(rxn_id))
            rxn_id_ir = rev2irrev{rxn_id}.conv(1);
            sel_m(rxn_id_ir) = 1;
            sel_m(rxn_id_ir+1) = -1;
        else
            % Irrev rxns
            rxn_id_ir = rev2irrev{rxn_id}.conv;
            sel_m(rxn_id_ir) = 1;
        end
        % Figure out ordering in decoupled representation
        ord_ir = [ord_ir rxn_id_ir];
    end
end
% Get ordering indices
[tmp,ord_ind] = sort(ord_ir);
% Reorder or create weights
if (sum(wt_meas_rxn) == 0)
    wt_m = ones(n_m,1);
else
    wt_m = wt_meas_tmp(ord_ind);
end
% Reorder measured flux values
b_m = b_meas_tmp(ord_ind);

%****************************
% Create OMNI model matrices
%****************************

% Set objectives for continuous and integer variables (currently all zeros)
% This can be used to set an OptKnock style objective
% Note that these objectives will be added to the actual OMNI objective
% (minimization of the distance between measured and predicted fluxes)
clp = zeros(n,1);
cip = zeros(sum(sel_int_ir),1)/sum(sel_int_ir);

% Create the matrices for the OMNI bilevel MILP
[A_bl,b_bl,c_bl,csense_bl,lb_bl,ub_bl,vartype_bl,sel_cont_sol,sel_int_sol] = ...
    omni(clp,cip,S_ir,c_ir,ub_ir,sel_int_ir,rev_int,sel_ic,b_ic,csense_ic,...
    sel_m,b_m,wt_m,H,K,ksense,sel_int_prev_ir);

% The variables created by the omni function can be used to solve
% the OMNI optimization problem using any solver capable of solving
% MILP problems. The constraint senses etc are written in a form
% compatible with the LINDO solvers.
% Here we use the CPLEX solver called through GAMS.

% Create GAMS input file
if (create_gams)
    
    % Clean up names so that GAMS accepts them
    intvarnames = strrep(rxns_ir(sel_int_ir == 1),'/','');
    cnt = 0;
    tmp = [rxns_ir(sel_m==1);rxns_ir(sel_m==1)];
    for i = 1:(length(c_bl)-length(intvarnames))
        if (i <= length(rxns_ir))
            contvarnames{i} = rxns_ir{i}; 
        elseif (i > (length(c_bl)-length(intvarnames))-2*length(b_m))
            cnt = cnt + 1;
            if (cnt <= length(b_m)) 
                contvarnames{i} = ['eu_' tmp{cnt}];
            else
                contvarnames{i} = ['el_' tmp{cnt}];
            end
        else
            contvarnames{i} = ['x' num2str(i)];
        end
    end
    contvarnames = strrep(contvarnames,'/','');
    contvarnames = strrep(contvarnames,'''','a');
    contvarnames = strrep(contvarnames,'(','_');
    contvarnames = strrep(contvarnames,')','');
    contvarnames = strrep(contvarnames,'.','');
    % Print out GAMS input file
    create_gams_input(A_bl,b_bl,c_bl,lb_bl,ub_bl,vartype_bl,csense_bl,gams_file,intvarnames,contvarnames);    
end

% Solve problem using CPLEX through GAMS
if (solve_omni)
  if (copyfile(gams_file,'omni_model.gms'))
    % Run GAMS
        system('gams omni_driver.gms');  
        if (exist('omni.stat','file'))
            % Parse solution files
            gams_omni_sol = parse_gams_omni_output('omni');
            [rxn_found,rxn_ind] = ismember(gams_omni_sol.flux_rxns,rxns_ir);
            omni_sol.x = zeros(size(rxns_ir));
            omni_sol.x(rxn_ind(rxn_found)) = gams_omni_sol.flux_vector(rxn_found);
            kos = gams_omni_sol.kos;
            kos = strrep(kos,'_f','');
            kos = strrep(kos,'_b','');
            kos = unique(kos);
            omni_sol.kos = kos; 
            omni_sol.rxns = rxns_ir;
            omni_sol.stat = gams_omni_sol.status;
            omni_sol.f = gams_omni_sol.obj_val;
            % Remove temporary files
            delete('omni_model.gms');
            delete('omni_kos.out');
            delete('omni_fluxes.out');
            delete('omni.stat');
            delete('omni_driver.lst');
            delete('omni_example.gms');
        else
            fprintf('Error: GAMS failed to produce an output file\n');
        end
    else
        fprintf(['Error: Copy file ' gams_file ' failed\n']);
    end
else 
  omni_sol = [];
end

