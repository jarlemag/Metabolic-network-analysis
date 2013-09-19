%Bounds
%Set max values for all fluxes in fixed model to 200% of maximum glucose uptake rate
%(Slightly different than described by Schuetz et al. 


limit = -2*model.lb(2);
for i = 12:85 %Exclude exchange reactions
model.ub(i) = limit; 
if model.lb(i) ~=0
    model.lb(i) = -limit;
end


%To constrain to a maximum of the *actual* glucose uptake rate, it might be possible to run
%an interative script: Run FBA -> Check if any fluxes are above 200 % of
%glucose uptake rate, if yes set max flux to 200 % of current glucose
%uptake -> Repeat until constraint is obeyed.

limit = -2*model.lb(2); %Limit all internal fluxes to max 200 % of GLC uptake. 
for i = 12:85 %Exclude exchange reactions

model.ub(i) = limit; 
if model.lb(i) ~=0
    model.lb(i) = -limit;
end
constraintment = false;
while constraintmet == false
    result = optimizeCbModel(model)
    m = 1;
    for j = 12:size(model.rxns)
    if abs(result.x(j)) > limit
        m = 0
    else
    end
    if m ==1
        constraintmet = true;
    else
    end
    end
    finalresult = result;
    
end