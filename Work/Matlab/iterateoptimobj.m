%iterateoptimobj.m

%Make a bar plot of euclidean distance to experimental fluxes (aerobe batch
%culture) for all three models, with following objectives:
%100 random objectives (average of distance), BM. ATP max.


%PROBLEM: ALL objectives are identified as local minima. Possibly because
%of the simplex-algorithm used in solving FBA?

%Possible to use trust-region algorithm instead?

%Generate locally optimized objective functions from random objective
%functions

Fmodel = readCbModel('fixed_schuetz_v4');
Cmodel = readCbModel('combinedcore');
Gmodel = readCbModel('iJO1366');

%n = 86;

m = 100; %number of trials
objvectors = zeros(n,1);


%Expected run time: 2+ minutes

%Test random objectives, find the euclidean distance to experimental data,
%then perform optimization of the random objective with respect to
%euclidean distance:

Fresult = optimizeCbModel(Fmodel);
startdistanceF_BM = eucdistFCG(Fresult.x,1,1,1);

Cresult = optimizeCbModel(Cmodel);
startdistanceC_BM = eucdistFCG(Cresult.x,2,1,1);

Gresult = optimizeCbModel(Gmodel);
startdistanceG_BM = eucdistFCG(Gresult.x,3,1,1);


%Set objective to max ATP in fixed model:
atpOBJrxns ={'pgk','pykAF','sucCD','atp','ackAB_tdcD_purT'};
  Fmodel = changeObjective(Fmodel,atpOBJrxns)
  Fmodel = changeRxnBounds(Fmodel,'biomass',0.6,'b')
Fresult = optimizeCbModel(Fmodel);
startdistanceF_ATP = eucdistFCG(Fresult.x,1,1,1);



%Set objective to maximization of ATP producing fluxes in ECM model:
  atpOBJrxns ={'ATPS4r','pyk'};
  Cmodel = changeObjective(Cmodel,atpOBJrxns)
  Cmodel = changeRxnBounds(Cmodel,'Biomass_Ecoli_core_w_GAM',0.6,'b')  
  Cresult = optimizeCbModel(Cmodel);
startdistanceC_ATP = eucdistFCG(Cresult.x,2,1,1);


%Set objective to max ATP in iJO1366:
atpOBJrxns ={'AP5AH','ATPS4rpp','PYK'};
Gmodel = changeObjective(Gmodel,atpOBJrxns);
startdistanceG_ATP = eucdistFCG(Gresult.x,3,1,1)


startdistancesFrand=0;
for i = 1:m
    c0 = rand(size(Fmodel.c,1),1);
    Fmodel.c = c0;
    result = optimizeCbModel(Fmodel);
    startdistancesFrand = horzcat(startdistancesFrand,eucdistFCG(result.x,1,1,1)); %Calculate the euclidean distance between experimental data and the FBA result for the starting objective vectors
    %minobj = optimobj(model,1,1,1);
    %objvectors = horzcat(objvectors,minobj);
end
    
startdistancesCrand=0;
for i = 1:m
    c0 = rand(size(Cmodel.c,1),1);
    Cmodel.c = c0;
    result = optimizeCbModel(Cmodel);
    startdistancesCrand = horzcat(startdistancesFrand,eucdistFCG(result.x,1,1,1)); %Calculate the euclidean distance between experimental data and the FBA result for the starting objective vectors
    %minobj = optimobj(model,1,1,1);
    %objvectors = horzcat(objvectors,minobj);
end

 startdistancesGrand=0;
for i = 1:m
    c0 = rand(size(Gmodel.c,1),1);
    Gmodel.c = c0;
    result = optimizeCbModel(Cmodel);
    startdistancesGrand = horzcat(startdistancesFrand,eucdistFCG(result.x,1,1,1)); %Calculate the euclidean distance between experimental data and the FBA result for the starting objective vectors
    %minobj = optimobj(model,1,1,1);
    %objvectors = horzcat(objvectors,minobj);
end

mindistances = 0;

%Use the locally optmized objectives to perform FBA and calculate the
%euclidean distance to the experimental data:

%{
for i = 1:m
    model.c = objvectors(:,i+1)
    result = optimizeCbModel(model)
    mindistances = horzcat(mindistances,eucdistFCG(result.x,1,1,1))
end
%}
     


%plot(mindistances)
hold on



bar([1:9],[mean(startdistancesFrand) mean(startdistancesCrand) mean(startdistancesGrand) startdistanceF_BM startdistanceC_BM startdistanceG_BM startdistanceF_ATP startdistanceC_ATP startdistanceG_ATP]) 
set(gca,'XTickLabel',{'SHUETZR Random' 'ECM Random' 'iJO1366 Random' 'SCHUETZR BM' 'ECM BM' 'iJO1366 BM' 'SCHUETZR ATP' 'ECM ATP' 'iJO1366 ATP'})


%bar[1 2 3
% bar([1 2 3],[average BMdis ATPdis])
%set(gca,'XTickLabel',{'Average, random objectives' 'biomass' 'atp'})