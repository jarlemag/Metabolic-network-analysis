function result = EXreactinfo(model)
% reactinfo - extracts the name, lower bound, and upper bound 
%             for each exchange reaction in the FBA model.
%
% Usage:
%  result = EX_reactinfo(model)
% 
%  model       a COBRA formatted metabolic model
%  result      a matrix with
%                       1st column = short reaction name
%                       2nd column = corresponding lower bound
%                       3rd column = corresponding upper bound
%
%  Eivind Almaas 8/4/2010

[nrows,ncols] = size(model.rxns);
test1 = cell(nrows,3);
for i=1:nrows
    test1(i,:) = [model.rxns(i) model.lb(i) model.ub(i)];
end

test2 = findExcRxns(model,0,0);
EX    = find(test2 > 0.5);
EXno  = size(EX,1);

result = cell(EXno,3);
for i=1:EXno
    result(i,:)=test1(EX(i),:);
end
