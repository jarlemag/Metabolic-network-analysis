function result = printShadowPrice(model,shadow)
% printShadowPrice - prints the name of the metabolite for which 
%                    the shadow price (aka. dual variable) is associated
%
% USAGE:
%  result = printShadowPrice(model,shadow)
%
% INPUTS:
%  model        a COBRA formatted metabolic model
%  shadow       the output variable with flux results from optimization
%
% OUTPUTS:
%  result       a cell matrix with metabolite names and their shadow prices
%
%
%  Eivind Almaas 24. March 2011
%

[nrows,ncols] = size(model.mets);
result = cell(nrows,2);
for i=1:nrows
    result(i,:) = [model.mets(i) shadow.y(i)];
end
