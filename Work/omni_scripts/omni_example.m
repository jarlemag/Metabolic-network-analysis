% OMNI example for E. coli evolved tpiE1 strain
%
% Markus Herrgard 2/8/06

% Load metabolic model
load iJR904_aer_gmm_omni.mat

% Experimental flux data (in this case only growth rate and exchange
% fluxes)
exp_flux_rxns = {'BiomassEcoli','EX_glc(e)','EX_ac(e)','EX_pyr(e)'};
exp_fluxes = [0.51 -7.8 1 0]';

% Number of iterations to find alternate optimal solutions
max_omni_iter = 1;

% Add the function library to the path
path(path,[pwd '/lib_files']);

% Remove the TPI reaction (corresponding to a tpi knockout)
model = remove_rxns(model,'TPI');

% Apply the exchange flux constraints
for i = 2:length(exp_fluxes)
    rxn_id = find(strcmp(model.rxns,exp_flux_rxns{i}));
    model.lb(rxn_id) = exp_fluxes(i);
    model.ub(rxn_id) = exp_fluxes(i);
end

% Store options for OSI
opt.K = 1; % Number of bottlenecks sought
opt.ksense = 'E'; % Type of K constraint (E/G/L)
opt.H = 100; % Maximum flux
opt.create_gams = true; % Create GAMS output
opt.gams_file = 'omni_example.gms'; % Name of GAMS input file
opt.solve_omni = true; % Solve problem using GAMS

%  Set constraining fluxes for the overall problem 
opt.sel_icc_rxn = exp_flux_rxns; % Names of constraining fluxes
opt.b_icc_rxn = exp_fluxes; % Flux constraints
opt.b_icc_rxn(1) = 0.3; % Minimum biomass constraint
opt.csense_icc = 'GEEE'; % Type of constraint (G/L/E)

% Set measured reaction fluxes for the OMNI objective function 
% (in this case just the growth rate)
opt.sel_meas_rxn = {exp_flux_rxns{1}}; % Names of measured fluxes
opt.b_meas_rxn = exp_fluxes(1); % Target value for measured fluxes
opt.wt_meas_rxn = 1; % Weights for measured fluxes

% Run OMNI iteratively to search for alternative solutions
previous_solutions = [];
for omni_iter = 1:max_omni_iter

  % Set up and solve OMNI problem
  omni_sol = run_omni(model,opt,false);

  % Print out results
  fprintf('Bottlenecks identified:\n');
  for i = 1:length(omni_sol)
    fprintf('%s\n',omni_sol.kos{i});
  end
  fprintf('OMNI objective value: %6.3f\n',omni_sol.f);
  
  % Store results
  omni_status(omni_iter) = omni_sol.stat(1);
  obj_values(omni_iter) = omni_sol.f;
  all_solutions{omni_iter} = omni_sol.kos;
  
  % Add to previously found results
  [tmp,ko_rxn_ind] = ismember(omni_sol.kos,model.rxns);
  previous_solutions(:,omni_iter) = zeros(length(model.rxns),1);
  previous_solutions(ko_rxn_ind,omni_iter) = 1;

  opt.sel_int_prev = previous_solutions;
  
end
