%extractflux.m. Extracts/calculates chosen experimental fluxes from a flux vector. 
%Supecercedes, contains and may call on extractfluxFCG.m
%fluxes = extractflux(fluxvector,options)


function fluxes = extractflux(fluxvector,options)


if exist('options','var') == 0
    options = struct; 
    %Make an options structure if it was not given in the function call
end

load('reactionmaps.mat')
%Load reaction maps

if isfield(options,'extractmode') == 0
    options.extractmode = 'default';
    %Set the default operation mode.
end

if isfield(options,'usealtmap') == 0
    options.usealtmap =0;
end

if strcmp(options.extractmode,'default') == 1
    options.extractmode = 3;
elseif strcmp(options.extractmode,'hardcoded') == 1
    %operation mode is 'hardcoded', passing the function call on to
    %the original extractfluxFCG function located below this functio
    options.extractmode = 1;
elseif strcmp(options.extractmode,'simplemap') ==1
     options.extractmode = 2;
elseif strcmp(options.extractmode,'tokens') ==1
    options.extractmode =3;
end

if isfield(options,'model_id') == 0
    options.model_id = 1;
    %Default model
    disp('WARNING: Extractflux.m: Model not specified! Model ID set to 1 by default.')
end

if isfield(options,'verbflag') == 0
    options.verbflag =0;
end

if isfield(options,'debugmode') ==0
    options.debugmode =0;
end

debugmode = options.debugmode;
verbflag = options.verbflag;
model_id = options.model_id;

if debugmode >0
    disp('Extractflux.m: Debugmode:')
    disp(debugmode)
    disp('Extractflux.m: extractmode:')
    disp(options.extractmode)
end



if  options.extractmode == 1
    
    if debugmode > 0
        disp('extractflux.m: Calling extractfluxFCG')
    end
        
    fluxes = extractfluxFCG(fluxvector,options.model_id);
    if verbflag == 2
        disp('extractflux.m: Using hardcoded reaction relations.')
    end
    return
elseif options.extractmode == 2
    
     if verbflag > 1 || debugmode >0
        disp('extractflux.m: Using simple reaction map.')
    end
    
    switch model_id 
        case 1
            reactionmap = reactionmaps.Fmap;
        case 2
            reactionmap = reactionmaps.Cmap;
        case 3
            reactionmap = reactionmaps.Gmap;
    end
    
elseif options.extractmode == 3
    
    if verbflag >1 || debugmode >0
        disp('extractflux.m: Using token reactions.')
    end
    
     switch model_id 
        case 1
            reactionmap = reactionmaps.Fmap2;
        case 2
            reactionmap = reactionmaps.Cmap2;
        case 3
            reactionmap = reactionmaps.Gmap2;
    end
    
end

if options.usealtmap == 1
    if isfield(options,'altmap') == 1
        disp('extractflux: Using provided alternate reaction map')
        reactionmap = options.altmap;
   else
       warning('extractflux:map','Alternate reaction map not found. Default reaction map in use.')
   end
end



n = size(reactionmap,1); %Number of experimental fluxes 
fluxes = zeros(n,1); %Make vector to store output values in
j = size(reactionmap,2); %number of columns in reaction map matrix
modreactlist = zeros(n,j-1);
indfluxes = zeros(n,j-1);
for i = 1:n
    modreactlist(i,:) = reactionmap(i,2:j);
    for k = 2:j
        if reactionmap(i,k) ~= 0
            indfluxes(i,k-1) = fluxvector(abs(reactionmap(i,k)))*sign(reactionmap(i,k));
        end
    end
fluxes(i) = sum(indfluxes(i,:));
end
    
end
   
%extractfluxFCG.m:
%Extract a subset of fluxes from FBA result for comparison with
%experimentally determined values.

%APPLICABLE TO:
%Schuetz "fixed" model (model_id = 1)
%Ecoli combined core model (model_id = 2);
%iJO1366 (model_id = 3);

%Input:

%fluxvector:FBA flux result vector from COBRA toolbox:
%model_id: model identifier 

%Output: A vector of 28 fluxes corresponding to the 28 experimentally
%determined fluxes used in the article by Schuetz et al. (2007).

%Column 1: Flux value
%Column 2: Reaction # in model

%syntax: fluxFCG = extractfluxFCG(fluxvector,model_id)



function fluxFCG = extractfluxFCG(fluxvector,model_id)

fluxFCG = zeros(28,1);

switch model_id
    
    case 1
        

fluxFCG(1,1) = fluxvector(12);        
%fluxFCG(1,1) = fluxvector(12)+fluxvector(72); %glk reaction (GLC + ATP -> G6P) + ptsGHI reaction

fluxFCG(2,1) = fluxvector(28); %zwf reaction (directly coupled to pgl reaction)
                     
fluxFCG(3,1) = fluxvector(30); %gnd reaction (6PGC -> NADPH + CO2 + RL5P)

fluxFCG(4,1) = fluxvector(13); %pgi reaction

fluxFCG(5,1) = fluxvector(69); %edd reaction

fluxFCG(6,1) = fluxvector(15); %fbaAB


%Experimental reaction: F6P + ATP -> 2T3P

%Model reactions
%pfkA/pfkB: F6P + ATP -> FDP

%fbp & glpX: FDP (FBP) -> F6P

%fbaA/fbB: FDP <-> DHAP + GA3P

%The experimental flux described is the net conversion of
%Fructose 6-phosphate (F6P)to Glyceraldehyde 3-phosphate (GAP
%or GA3P). The conversion of F6P to Fructose-1,6-bisphosphate (FBP or FDP)
%by phosphofructokinase (pfkA/pfkB)is an intermediate step in this process. 
%fbp and glpX reactions convert FDP back to F6P, so these fluxes must be substracted.   
%FBP is split by fructose bisphosphate aldolase (fbaA/fbaB) to Dihydroxyacetone (DHAP or DHP) and
%GAP. pfkA + fpkB - ( fbp + glpX) should equal (fbaA + fbaB), as there are no other reactions involving FBP. DHP is then converted to GAP by isomerization, catalyzed by TpiA.
%However, DHP can also be consumed by the MgsA reaction. If this reaction
%is active, the comparison between the model fluxes and experimental fluxes becomes more problematic.
%However, the mgsA reaction is inactive in the experiments considered, and also appears to be inactive when optimizing for biomass production. It is therefore ignored here.
%If mgsA reaction is zero, all DHAP should be converted to GA3P. (or vice
%versa)


fluxFCG(7,1) = fluxvector(33); %tktAB reaction

fluxFCG(8,1) = fluxvector(34); %tktAB_R2

fluxFCG(9,1) = fluxvector(35); %talAB

fluxFCG(10,1) = fluxvector(17); %reaction gapA

fluxFCG(11,1) = fluxvector(19);   %reactions gpmAB


%Model reaction (gpmA/gpmB) refers to conversion of 3PG to 2PG. Experimental reaction
%refers to conversion of 3PG ("PGA") to PEP. (assumed to be directly
%linked).

fluxFCG(12,1) = fluxvector(21)-fluxvector(27); %reactions pykF, pykA (pyruvate kinase) and pps (PEP synthase)

fluxFCG(13,1) = fluxvector(22); %reaction aceEF (pyruvate + coA -> acetyl-CoA + NADH + CO2)

fluxFCG(14,1) = fluxvector(36); %gltA/prpC (Acetyl-CoA + OAA -> CoA + Citrate)

fluxFCG(15,1) = fluxvector(39); %reaction icd (isocitrate dehydrogenase)


fluxFCG(16,1) = fluxvector(40); %reaction sucAB

%Have a closer look at this one. Not a direct match between model and
%experimental reactions.
%experimental reaction: OGA -> FUM + CO2 + 1.5ATP + 2NADH

%In Shuetz model, Alpha-ketoglutarate (OGA, or AKG) participates only in
%three reactions. It is produced from isocitrate in the (reversible) reaction icd and reacts
%with coenzyme A (CoA) to form CO2, NADH and succinyl-CoA in reaction
%sucAB. It is also drained in the biomass reaction.

%Succinyl-CoA can react in reaction sucCD to produce ATP, CoA and
%Succinate. Succinyl-CoA only appears in these two reactions, so the
%reactions sucAB and sucCD are directly coupled.

%Succinate can react in reaction frdABCD or sdhAB (identical reactions) to produce fumarate and FADH.
%The reaction rates of frdABCD/sdhAB are not coupled to sucAB/sucCD, as
%succinate can also be produced several other ways, including directly from isocitrate by the glyxocylate
%shunt. The reactions sdhAB and sdhABCD also produce/consume succinate.

%The coupling of reactions sucAB, sucCD and frdABCD thus gives the reaction
%equation: AKG -> Fumarate + CO2 + ATP + FADH

%Conclusion: Comparing the experimental flux to that of sucAB/sucCD seems
%most reasonable.

fluxFCG(17,1) = fluxvector(43); % fumABC

fluxFCG(18,1) = fluxvector(44)+fluxvector(45); % reactions mdh & mqo

fluxFCG(19,1) = fluxvector(23)+fluxvector(24); %reactions maeB & maeA (sfcA)

fluxFCG(20,1) = fluxvector(25); %reaction pck (PEP carboxykinase)

fluxFCG(21,1) = fluxvector(48); %reaction ppc

fluxFCG(22,1) = fluxvector(64)-fluxvector(66); %reactions pta and acs

%Experimental reaction: AcCoA -> Acetate + ATP
%pta catalyzes conversion of acetyl-CoA to acetyl-phosphate (ACTP), while reactions ackA, ackB, tdcD and purT are coupled to pta, consuming the ACP produced)
%Acetyl-CoA synthase (acs) catalyzes the synthesis of Acetyl-CoA from
%acetate and ATP (anc CoA), so that reaction should be substracted. Or
%should it?

%The reaction ACKr in ECME is analogous to the reverse of
%ackA/ackB/tdcD/puRT in SCHUETZR.

fluxFCG(23,1) = -fluxvector(49)+fluxvector(50); %reactions pntAB and udhA. Reaction pntAB is the reverse of the experimental reaction.

fluxFCG(24,1) = fluxvector(74); %reaction O2 (transport). Experimental reaction is respiration.

fluxFCG(25,1) = fluxvector(85); %biomass reaction

fluxFCG(26,1) = fluxvector(46); %linked reactions aceA/aceB 

%Glyxocylate shunt. Experimental reaction is ICT + AcCoA -> MAL + FUM +
%NADH. Model reactions are ICIT <-> GLX + SUCC (aceA) and GLX + ACCOA ->
%MAL + COA (aceB). Equivalence between model and experimental flux is only
%valid if all succinate produced is converted to fumarate. If succinate for example is
%instead exported, the equivalence fails. The equivalence most likely holds
%in most cases.


%fluxF(27,1) = fluxvector(67)+fluxvector(68) %reactions dld & ldhA, together they equal mgsA.

%Experimental reaction is described as "DHP -> Pyr", but does not specify
%intermediates.

%Decided to use mgsA flux directly:

fluxFCG(27,1) = fluxvector(71); %mgsA

fluxFCG(28,1) = fluxvector(77); %ethanol secretion

%END OF CASE 1
    case 2
        
        fluxFCG(1,1) = fluxvector(96)+fluxvector(50); 
%reactions glk and pts

fluxFCG(2,1) = fluxvector(48); %G6PDH2r reaction

fluxFCG(3,1) = fluxvector(57); %gnd reaction

fluxFCG(4,1) = fluxvector(74); %pgi reaction

fluxFCG(5,1) = fluxvector(98); %edd reaction

fluxFCG(6,1) = fluxvector(40); %fba reaction

fluxFCG(7,1) = fluxvector(93); %tkt1

fluxFCG(8,1) = fluxvector(94); %tkt2

fluxFCG(9,1) = fluxvector(91); %talA

fluxFCG(10,1) = fluxvector(49); %gapD

fluxFCG(11,1) = -fluxvector(77); %pgm

fluxFCG(12,1) = fluxvector(83)-fluxvector(81); %pyk - pps

fluxFCG(13,1) = fluxvector(71); %PDH

fluxFCG(14,1) = fluxvector(15); %Citrate Synthase

fluxFCG(15,1) = fluxvector(59); %ICDhyr

fluxFCG(16,1) = fluxvector(8); %AKgdh

fluxFCG(17,1) = fluxvector(46); %FUM

fluxFCG(18,1) = fluxvector(64); %MDH

fluxFCG(19,1) = fluxvector(65)+fluxvector(66); %ME1 + ME2

fluxFCG(20,1) = fluxvector(80); %ppck

fluxFCG(21,1) = fluxvector(79); %ppc

fluxFCG(22,1) = fluxvector(82) - fluxvector(101); %PTAr, %acs

fluxFCG(23,1) = fluxvector(68)-fluxvector(92); %NADTRHD, THD2

fluxFCG(24,1) = fluxvector(70); %O2t

fluxFCG(25,1) = fluxvector(13); %biomass

fluxFCG(26,1) = fluxvector(60); %icl

fluxFCG(27,1) = fluxvector(100); %mgsA

fluxFCG(28,1) = fluxvector(24); %ethanol exchange

    case 3

fluxFCG(1,1) = fluxvector(1500)+fluxvector(1353);
%Reaction "HEX1" (#1500) of iJO1366 ~ reaction "glk" of Schuetz model
%+ reaction glcptspp (#1353)

fluxFCG(2,1) = fluxvector(1283);

%reaction "R_G6PDH2r" (#1283) (coupled to "R_PGL" (#2079) of iJO1366 ~ R_zwf in the Schuetz model.


fluxFCG(3,1) = fluxvector(1423); %gnd reaction

fluxFCG(4,1) = fluxvector(2077); %pgi reaction

fluxFCG(5,1) = fluxvector(1095); %edd reaction

fluxFCG(6,1) = fluxvector(1151); %fba reaction

%Note that in iJO366, Fructose 1,6-bisphosphate (FBP) is also produced from fructose 1-phosphate by
%fructose-1-phosphate kinase. The experimental reaction described by
%Schuetz et al is "F6P + ATP -> 2T3P", where the intermediate steps
%involve FBP. Thus any reactions involving FBP may affect the results. I
%have not investigated exactly how the experimentally found flux is
%defined/calculated. This should be done for a rigorous analysis.


fluxFCG(7,1) = fluxvector(2442); %tkt1 reaction (schuetz: tktA/tktB)

fluxFCG(8,1) = fluxvector(2443); %tkt2 reaction (schuetz tktA_R2/tktB_R2)

fluxFCG(9,1) = fluxvector(2398); %tala (transaldolase) reaction 

fluxFCG(10,1) = fluxvector(1315); %GAPD reaction

fluxFCG(11,1) = -fluxvector(2081); %pgm reaction (shuetz: gpmA/gpmB)


%Model reaction (pgm) refers to interconversion of 3PG/2PG. Experimental reaction
%refers to conversion of 3PG ("PGA") to PEP. (assumed to be directly
%linked).

fluxFCG(12,1) = fluxvector(2266)-fluxvector(2198); %reactions "pyk" (schuetz: pykA/pykF) and PPS

fluxFCG(13,1) = fluxvector(2047); %reaction PDH (Schuetz: aceEF)

fluxFCG(14,1) = fluxvector(878); %reaction CS (Schuetz: gltA/prpC)

fluxFCG(15,1) = fluxvector(1548); %reaction ICDHyr

fluxFCG(16,1) = fluxvector(608); %reaction akGDH (Schuetz: sucAB)
%Possible alternative: R_SUCOAS (Schuetz: sucCD)
%See notes in extractfluxF.m

fluxFCG(17,1) = fluxvector(1245); %reaction FUM

fluxFCG(18,1) = fluxvector(1758)+fluxvector(1759)+fluxvector(1760); %reactions MDH, MDH_2 and MDH_3

fluxFCG(19,1) = fluxvector(1761)+fluxvector(1762); %reactions ME1 and ME2 (Schuetz: maeA and maeB)

fluxFCG(20,1) = fluxvector(2183); %reaction PPCK (Schuetz: pck)

fluxFCG(21,1) = fluxvector(2181); %reaction ppc

fluxFCG(22,1) = fluxvector(2240)-fluxvector(547);     %reaction PTAr and ACS

fluxFCG(23,1) = fluxvector(1871)-fluxvector(2418);  %reactions NADTRHD (Schuetz: udhA) and THD2PP (Schuetz: pntAB)

%Experimental reaction: NADPH->NADH

%R_THD2pp: NAD(P) transhydrogenase (periplasm)
%NAD + NADPH -> NADH + NADP+

%R_NADTRHD: NAD transhydrogenase
%NADH + NADP+ + H+(p) -> NADPH + NAD+ + H+(c) 

%For a discussion of udhA and pntAB, see http://www.ncbi.nlm.nih.gov/pubmed/14660605

%fluxG(24,1) = -fluxvector(252) ;%Oxygen uptake.

%Experimental reaction: O2 + 2NADH -> 2P/OxATP
%Simplification used here: All oxygen is used for the above reaction.
%This approximation should be evaluated, or the reaction above replaced
%with a more appropriate reaction after investigating how respiration is
%modelled in JO1366.

%Done below:

fluxFCG(24,1) = fluxvector(914)+fluxvector(915)+fluxvector(916); %CYTBD2pp,CYTBDpp,CYTBO3_4pp



%fluxF(25,1) = fluxvector(7); %biomass reaction (Wildtype)
fluxFCG(25,1) =fluxvector(8); %biomass reaction (Core)

fluxFCG(26,1) = fluxvector(1552); %reaction ICL (Schuetz: AceA)

%Glyxocylate shunt

fluxFCG(27,1) = fluxvector(1790); %reaction MGSA

%Experimental reaction: DHAP -> Pyruvate

%MgsA catalyzes conversion of DHAP to methylglyoxal. Here it is assumed
%that all the methylglyoxal produced is further converted to lactate.


fluxFCG(28,1) = fluxvector(124); %Ethanol secretion

end


end

%Abbreviations used in supplementary material by Schuetz et al. 2007:

%FDP = FBP = fructose 1,6-diphosphate (fructose-1,6-bisphosphate) 
%2PG = 2-phosphoglycerate
%3PG (=PGA) = 3-phosphoglycerate 
%
%ICT = Isocitrate
% P5P = ribose 5-phosphate
% T3P = GA3P = Glyceraldehyde Triphosphate
% S7P =  sedoheptulose 7-phosphate
% X5P = xylulose 5-phosphate
% E4P = 
% F6P =
% OGA = oxoglutarate (2-oxoglutarate = alpha-ketoglutarate)
% AKG = alpha-ketoglutarate
%DHP/DHAP = Dihydroxyacetone phosphate
%6GP = 6-phosphogluconate
%KDG = 2KD6PG = 2-keto-3-deozy-6-phosphogluconate




