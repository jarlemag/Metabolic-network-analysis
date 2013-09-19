%AnalyseHFB.m

isHFB = zeros(28:1);


fixedmodel = readCbModel('fixed_schuetz_v4');
combinedmodel = readCbModel('combinedcore');
genomemodel = readCbModel('iJO1366');


fixedresult = optimizeCbModel(fixedmodel);
combinedresult = optimizeCbModel(combinedmodel);
genomeresult = optimizeCbModel(combinedmodel);

[numb res names] = HighFluxBackbone(fixedmodel,fixedresult.x)
isHFB = horzcat(isHFB,abs(extractfluxFCG(res,1)))

[numb res names] = HighFluxBackbone(combinedmodel,combinedresult.x)
isHFB = horzcat(isHFB,abs(extractfluxFCG(res,2)))

%[numb res names] = HighFluxBackbone(genomemodel,genomeresult.x)
%isHFB = horzcat(isHFB,extractfluxFCG(res,3))


%Change objective of fixed model to max ATP:
 atpOBJrxns ={'pgk','pykAF','sucCD','atp','ackAB_tdcD_purT'};
 fixedmodel = changeObjective(fixedmodel,atpOBJrxns);
 fixedmodel = changeRxnBounds(fixedmodel,'biomass',0.6,'b');
  
 
%Change objective of ECM model to max ATP:

 atpOBJrxns ={'ATPS4r','pyk'};
 combinedmodel = changeObjective(combinedmodel,atpOBJrxns);
combinedmodel = changeRxnBounds(combinedmodel,'Biomass_Ecoli_core_w_GAM',0.6,'b')  ;


 %Change objective of iJO1366 model to max ATP:
 
 atpOBJrxns={'AP5AH','ATPS4rpp','pyk'};
genomemodel = changeObjective(genomemodel,atpOBJrxns);
genomemodel = changeRxnBounds(genomemodel,'Ec_biomass_iJO1366_core_53p95M',0.6,'b');


%Repeat analysis with ATP as objective


ATPfixedresult = optimizeCbModel(fixedmodel);
ATPcombinedresult = optimizeCbModel(combinedmodel);
ATPgenomeresult = optimizeCbModel(combinedmodel);



[numb res names] = HighFluxBackbone(fixedmodel,ATPfixedresult.x)
isHFB = horzcat(isHFB,abs(extractfluxFCG(res,1)))

[numb res names] = HighFluxBackbone(combinedmodel,ATPcombinedresult.x)
isHFB = horzcat(isHFB,abs(extractfluxFCG(res,2)))

%[numb res names] = HighFluxBackbone(genomemodel,genomeresult.x)
%isHFB = horzcat(isHFB,extractfluxFCG(res,3))

 
