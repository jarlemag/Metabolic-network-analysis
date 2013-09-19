function omni_sol = parse_gams_omni_output(file_base)
%PARSE_GAMS_OMNI_OUTPUT Parse GAMS output files for the OMNI method
%
% file_base     Base name for the output files
%
% Markus Herrgard 2/7/06

% Read solution status and objective values
fid=fopen([file_base '.stat']);
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    fields = split(tline,' ');
    if (~isempty(regexp(tline,'model status')))
        model_status = str2num(fields{3});
    elseif (~isempty(regexp(tline,'solver status')))
        solver_status = str2num(fields{3});
    elseif (~isempty(regexp(tline,'objective')))
        omni_sol.obj_val = str2num(fields{2});
    end
end
fclose(fid);

omni_sol.status = [model_status solver_status];

% Read KOs identified
fid=fopen([file_base '_kos.out']);
cnt = 0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    fields = split(tline,' ');
    if (~isempty(fields))
        cnt = cnt + 1;
        rxn_name = fields{1};
        ko_lst{cnt} = rxn_name;
    end
end
fclose(fid);

% Read fluxes and errors between measured and predicted fluxes
fid=fopen([file_base '_fluxes.out']);
cnt = 0;
cnt_e = 0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    fields = split(tline,' ');
    if (~isempty(fields))
        rxn_name = fields{1};
        if (regexp(rxn_name,'^EX_.*_e_[fb]'))
          rnx_name = strrep(rxn_name,'_e','(e)');
        end
        flux_val = str2num(fields{2});
        if (isempty(regexp(rxn_name,'^x\d')))
            if (~isempty(regexp(rxn_name,'^eu_')))
                cnt_e = cnt_e + 1;
                rxn_name = strrep(rxn_name,'eu_','');
                flux_err_rxns{cnt_e} = rxn_name;
                flux_errors(cnt_e) = -flux_val;
            elseif (~isempty(regexp(rxn_name,'^el_')))
                cnt_e = cnt_e + 1;
                rxn_name = strrep(rxn_name,'el_','');
                flux_err_rxns{cnt_e} = rxn_name;
                flux_errors(cnt_e) = flux_val;
            else
                cnt = cnt + 1;
                flux_rxns{cnt} = rxn_name;
                flux_vector(cnt) = flux_val;
            end
        end
    end
end
fclose(fid);

omni_sol.kos = ko_lst;
omni_sol.flux_rxns = flux_rxns';
omni_sol.flux_vector = flux_vector';
omni_sol.flux_err_rxns = flux_err_rxns';
omni_sol.flux_errors = flux_errors';

