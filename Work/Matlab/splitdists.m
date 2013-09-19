%comparesplits.m:
%Compare results produced when optimizing with respect to split ratios and
%raw fluxes.

options = struct;
options.objective = 2;
options.model_id =1;
options.usesplits = 1;

F= runsim(options);

opt = struct;
opt.compsplits = 1;
opt.model_id = 1;

splitsF = F.splitresult;
Fsplitdists = [splitsF.splitmindistance compdist(F.Fmin_minsol,opt)  compdist(F.gurobi_minsol,opt)];

Fmindists = [splitsF.splitrawmindistance F.Fmin_mindistance F.gurobi_mindist];

 
 options.model_id = 2;
 
 B = runsim(options);
 
 opt = struct;
 opt.compsplits = 1;
 opt.model_id = 2;
 
splitsB = B.splitresult;
Csplitdists = [splitsB.splitmindistance compdist(B.Fmin_minsol,opt)  compdist(B.gurobi_minsol,opt)]

 Cmindists = [splitsB.splitrawmindistance B.Fmin_mindistance B.gurobi_mindist]

 
 options.model_id = 3;
options.usesplits = 0;

options.usefmincon = 0;
C = runsim(options)
 
opt.model_id =3;
Gdist = C.gurobi_mindist
Gsplitdist = compdist(C.gurobi_minsol,opt)
Gsplitdists = [0 0 Gsplitdist];
Gmindists = [0 0 C.gurobi_mindist];

mat=vertcat(Fsplitdists,Csplitdists,Gsplitdists)
figure('name','Split ratio distance')
bar(mat)
xlabel('Model')
ylabel('Distance (dimensionless)')
legend('Minimization of split ratios','Minimization of raw fluxes (fmincon)','Minimization of raw fluxes (Gurobi)')
set(gca,'XTickLabel',{'SHUZETR','ECME','iJO1366b'})

mat2 = vertcat(Fmindists,Cmindists,Gmindists)
figure('name','Raw flux distance')
bar(mat2)
xlabel('Model')
ylabel('Distance (mmol/g*h)')
legend('Minimization of split ratios','Minimization of raw fluxes (fmincon)','Minimization of raw fluxes (Gurobi)')
set(gca,'XTickLabel',{'SHUZETR','ECME','iJO1366b'})