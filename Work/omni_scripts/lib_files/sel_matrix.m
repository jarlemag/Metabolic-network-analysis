function sel_mat = sel_matrix(sel_vec)
% SEL_MATRIX Create selection matrix from a selection vector
%
% sel_mat = sel_matrix(sel_vec)
%
% If sel_vec = [1 0 0 1 0 0]
%
% sel_mat = [1 0 0 0 0 0
%            0 0 0 1 0 0]
%
% For reversible selections
%
% If sel_vec = [1 0 0 1 -1 0]
%
% sel_mat = [1 0 0 0  0 0
%            0 0 0 1 -1 0]
%
% Markus Herrgard 3/28/03

nvar = length(sel_vec);
if (sum(sel_vec == -1) == 0)
    
    nsel = sum(sel_vec);
    isel = [1:nsel];
    jsel = find(sel_vec);
    sel_mat = sparse(isel,jsel,ones(nsel,1),nsel,nvar);
    
else
    
    sel_fw_ind = find(sel_vec == 1);
    sel_mat = sparse(length(sel_fw_ind),nvar);
    for i = 1:length(sel_fw_ind)
        sel_fw_id = sel_fw_ind(i);
        if (sel_vec(sel_fw_id+1) == -1)
            sel_mat(i,sel_fw_id) = 1;
            sel_mat(i,sel_fw_id+1) = -1;
        else
            sel_mat(i,sel_fw_id) = 1;    
        end
    end
    sel_mat = sparse(sel_mat);
    
end