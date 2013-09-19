function [rxns,mets,S,b,c,lb,ub,rev,sel_int] = unpack_fba_model(model)
%UNPACK_FBA_MODEL Unpack a FBA model structure
%
% [rxns,mets,S,b,c,lb,ub,rev,sel_int] = unpack_fba_model(model)
%
% Markus Herrgard 5/24/05

% Unpack the model
rxns = model.rxns;
mets = model.mets;
S = model.S;
b = zeros(size(model.mets));
c = model.c;
lb = model.lb;
ub = model.ub;
rev = model.rev;

if (isfield(model,'int'))
    sel_int = model.int;
else
    sel_int = [];
end
    