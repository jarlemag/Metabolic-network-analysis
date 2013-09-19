%optimobj.m

function cmin = optimobj(model,model_id,exp_id,startflag)


switch startflag
    case 0
c0 = zeros(size(model.c,1),1);
    case 1        
c0 = model.c;
end
        cmin = fminunc(@(c) optimobjsolve(c,model,model_id,exp_id),c0);

end