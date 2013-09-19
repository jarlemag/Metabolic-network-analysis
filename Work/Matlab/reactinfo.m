function result = reactinfo(model)
% reactinfo - extracts the reaction name, lower bound, and upper bound 
%             for each reaction in the FBA model.
%
% Usage:
%  result = reactinfo(model)
% 
%  model       a COBRA formatted metabolic model
%  result      a matrix with
%                       1st column = short reaction name
%                       2nd column = corresponding lower bound
%                       3rd column = corresponding upper bound
%
%  Eivind Almaas 8/4/2010

[nrows,ncols]=size(model.rxns);
result = cell(nrows,3);
for i=1:nrows
    result(i,:)=[model.rxns(i) model.lb(i) model.ub(i)];
end
