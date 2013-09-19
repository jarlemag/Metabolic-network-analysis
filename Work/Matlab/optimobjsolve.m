%optimobjsolve.m

function distance = optimobjsolve(c,model,model_id,exp_id)

testmodel=model;
testmodel.c = c;
result = optimizeCbModel(testmodel);
distance = eucdistFCG(result.x,model_id,exp_id);
end