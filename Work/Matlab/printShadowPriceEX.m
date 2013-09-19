function result = printShadowPriceEX(model,shadow)
% printShadowPrice - prints the name of the exchange metabolite for which 
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

nrows = size(model.mets,1);
test = cell(nrows,2);
j=0;
for i=1:nrows
    eMets = regexp(model.mets(i),'\[e\]');
    hasEs = ~isempty([eMets{:}]);
     if (hasEs == 1)
         j=j+1;
         test(j,:) = [model.mets(i) shadow.y(i)];
     end
end

result = test(1:j,:);
