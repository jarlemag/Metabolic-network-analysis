%computesplits.m:
%Combines and supercedes computesplitratiosFCG.m and
%computesplitratiosEXP.m
%
%syntax: splitratios = computesplits(fluxvector,model_id)
%model IDs follows scheme in runsim.m
%use model_id = 0 to compute split ratios from experimental data.
%
%Function for computing split ratios from raw flux
%ONLY APPLICABLE TO "SHUETZR", "ECME" & iJO1366 models, and experimental
%data from Schuetz et al. 2007
%Syntax:
% splitratios = computesplits(fluxvector,model_id)
%Default name using COBRA toolbox is result.x
%model_id: model identifier
%options:
%0: Experimental data
%1: SCHUETZR
%2: ECM
%3: iJO1366
%
%Output
% splitratios: Vector with computed split ratios
%
%First column: Split ratios calculated without taking biomass drain into
%account.
%
%Second column: Split ratios calculated after substracting the amount of
%metabolite being drained in the biomass reaction, from the sum of
%producing fluxes.


function splitratios = computesplits(fluxvector,model_id)

R = zeros(10,2);

switch model_id
    
    case 0
        
%See Schuetz et al. 2007, Supplementary files 2 (experimental flux values)
%and 3 (calculated split ratios).

%NB: Each split ratio is the ratio of one or more reactions consuming a
%metabolitee, to the total *producing* fluxes for that metabolite. In that
%regard, table 1 in Schuetz et al. 2007 is misleading and the equations
%there cannot be applied directly!

%For example, the split ratio R1 is defined there as:

% R1 = pgi/(Glk+Pts+Zwf+Pgi)

%Here, pgi is a reaction consuming the metabolite in question, Glucose
%6-phosphate (G6P), while Pts and Glk are producing reactions, and zwf is
%a reversible reaction either consuming or producing G6P. However, when the
%zwf flux is positive (consuming G6P), it should not be included in the
%denominator. Neither should pgi, as it is also a consuming reaction, and
%the calculation simplifies to R1 = pgi/(glk+pts)


%In short, a reaction should only be present in the denominator if it is a
%producing reaction for that metabolite.

% ---> Arrows denote the association between in silico and experimental
% fluxes in the supplementary material of Schuetz et al. 2007

R = zeros(10,1);

R(1) = fluxvector(4)/(fluxvector(1));

%From Table 1, Schuetz et al. 2007:
%R1 = pgi/(glk + pts +zwf + pgi) 


%Pts is not included (explicitly) in the experimental data. Considered
%together with glk. With zwf and pgi non-producing fluxes, we then get 
%R1 = pgi/(Glk+Pts)



R(2) = fluxvector(5)/fluxvector(2);

%R2 = edd/pgl

%edd -> experimental flux 5: 6PG -> T3P + PYR
%pgl -> experimental flux 2: G6P -> 6PG + NADPH

%In the schuetz model, the reactions edd and eda are directly coupled and
%together give the above reaction equation.

R(3) = fluxvector(27)/fluxvector(6);

%From Table 1, Schuetz et al. 2007: R3 = mgsA/(fbaA+fbaB+tpiA)

%mgsA -> experimental flux 27: DHAP -> PYR
%fbaA+fbaB -> experimental flux 6: F6P + ATP -> 2T3P
%TpiA: No experimental flux. See discussion of flux #6 in extractfluxF.m

%The calculation of this split ratio from the experimental values is
%problematic, but in the paper by Schuetz it is always zero, so they 'get
%away' with it!

R(4) = fluxvector(12)/(fluxvector(11)+fluxvector(20));

%Equation in Table 1, Schuetz et al. 2007:
%R4 = (pykA+pykF+pts)/(eno+pps+pckA)

%pykA/pykF (PEP -> PYR + ATP) -> experimental flux 12: PEP -> PYR + ATP
%pts (PEP + GLC(ex) -> PYR + G6P) -> experimental flux 1: GLC + ATP -> G6P
%This association seems somewhat dubious.

%eno (2PG <-> PEP) -> experimental flux 12: PEP -> PYR + ATP
%pps (PYR + ATP -> PEP) -> experimental flux 12: PEP -> PYR + ATP
%pckA ->experimental flux 20: OAA + ATP -> PEP + CO2

%R4 is the ratio of (PEP converted to pyruvate) to the amount of PEP
%producing fluxes.

%In the supplementary material, both eno - a producing flux - and pykA/F, a
%consuming flux, is referred to the same experimental flux! Is this an
%error? Otherwise it seems like circular logic.

%Better to take 
%eno -> experimental reaction 11: 3PG <-> 2PG
%In the model by Schuetz, the reactions 3PG <-> 2PG (gpmA/gpmB) and 2PG <->
%PEP (eno) are directly coupled.

%The reactions pykA/pykF and pps are exact opposites. It is therefore not
%possible to treat them separately. Only net fluxes can be considered.

%Thus we get the equation

%R(4) = fluxvector(12)/(fluxvector(11)+fluxvector(20))
% (2PG <-> PEP)/((PGA -> PEP)+(OAA + ATP -> PEP + CO2))

%This gives the same value as reported by Schuetz et al. for experimental fluxes from batch aerobe
%culture.

R(5) = fluxvector(13)/(fluxvector(12)+fluxvector(19)+fluxvector(27)+fluxvector(5));

%R5 = (aceE+pflB+tdcE)/(pykA+pykF+MaeA+MaeB+Dld+ldhA+eda+pts+pflb+tdcE)

%aceE-> experimental reaction 13: PYR -> AcCoA + CO2 + NADH
%pflB -> no experimentally determined flux
%tdcE -> no experimentally determined flux
%pykA/pykF -> experimental flux 12
%maeA/maeB -> experimental flux 19: MAL -> PYR + CO2 + NADH
%dld/ldhA -> experimental flux 27: DHAP -> PYR (Questionable association? Certainly not valid if lactate is taken up as a growth substrate.)
%eda -> experimental flux 5: 6PG -> T3P + PYR
%pts -> experimental flux 1 (or: no experimentally determined flux)

%Simplified equation:

%R5 = Ex13/(Ex12 +Ex 19 + Ex 27 + Ex5) 
%This gives the same value as reported by Schuetz et al.

%The lack of experimental reactions for pflB/tdcE raises a question about
%how the split ratio should be calculated from the computed fluxes. Note
%that pflB/tdcE is a reversible reaction.

R(6)= fluxvector(14)/fluxvector(13);

%From Table 1, Schuetz et al. 2007: 
%R6 = (gltA+prpC)/(aceE+pflB+TdcE+Acs+AdhE+MhpF) 

%gltA/prpC -> experimental flux 14: OAA + AcCoA -> ICT
%aceE-> experimental reaction 13: PYR -> AcCoA + CO2 + NADH
%pflB -> no experimentally determined flux
%tdcE -> no experimentally determined flux
%acs -> experimental flux 22: AcCoA -> AC + ATP (not a producing flux)
%adhE -> no experimentally determined flux
%MhpF -> no experimentally determined flux

%Simplified equation:
% R6 = ex14/(ex13)

%This gives the same value as reported by Schuetz et al. for aerobic batch
%culture.

R(7) = fluxvector(26)/fluxvector(14);
%Glyxocylate shunt


%AceA -> experimental flux 26: ICT + AcCoA -> MAL + FUM + NADH
%acnA+ancB -> experimental flux 14: OAA + AcCoA -> ICT

R(8) = fluxvector(20)/(fluxvector(18)+fluxvector(21));

%R8 = pckA/(mdh+mqo+ppc)

%pckA -> experimental flux 20: OAA + ATP -> PEP + CO2
%mdh+mqo -> experimental flux 18: MAL -> OAA + NADH
%ppc -> experimental flux 21: PEP + CO2 -> OAA

R(9) = fluxvector(22)/fluxvector(13);

%pta -> experimental flux 22: AcCoA -> AC + ATP
%aceE-> experimental reaction 13: PYR -> AcCoA + CO2 + NADH
%pflB -> no experimentally determined flux
%tdcE -> no experimentally determined flux
%acs -> experimental flux 22: AcCoA -> AC + ATP ( not a producing flux)
%adhE -> no experimentally determined flux
%MhpF -> no experimentally determined flux



R(10) = fluxvector(28)/fluxvector(13);

%R10 = (adhE+MhpF)/(aceE+pflB+TdcE+Acs+AdhE+MhpF)

%adhE -> no experimentally determined flux
%MhpF -> no experimentally determined flux

% No experimentally available fluxes for adhE/mhpF. Instead, use
% experimental flux for ethanol secretion.
%aceE-> experimental reaction 13: PYR -> AcCoA + CO2 + NADH
        
         case 1

%R1 (Flux into glycolysis)
%Focus metabolite: G6P
R(1,1) = fluxvector(13)/(subplus(fluxvector(12))+subplus(fluxvector(72))+subplus(-fluxvector(28)));
R(1,2) = fluxvector(13)/(subplus(fluxvector(12))+subplus(fluxvector(72))+subplus(-fluxvector(28))-0.33*fluxvector(85));

%Definition in Table 1, Schuetz et al. 2007:

%R1 = pgi/(glk+pts+zwf+pgi)

%the term '-0.33*fluxvector(85)' substracts the amount of G6P drained in
%the biomass reaction from the sum of fluxes producing G6P.

%Relevant reactions: pgi (#13), glk (#12), pts (#72), zwf (#28) 

%zwf catalyzes formation of 6GPL from G6P and is reversible, while pgl
%catalyzes formation of 6PGC from 6PGL, and is irreversible. The two
%reactions are directly coupled in the model, thus in practice zwf is
%irreversible as well.

%Subplus(-fluxvector(28)) returns the absolute value of the flux of zwf, if
%the flux of zwf is negative. As the positive direction of the zwf reaction
%is defined as consuming G6P, this ensures that the flux value of zwf is
%only counted if it is a producing reaction for G6P. Using subplus on the
%fluxes which are defined as producing reactions when their flux is
%positive also ensures that these are only counted when they are actually
%producing reactions

%If the net flux of the reaction(s) in the numerator is negative, the split
%ratio gives the proportion of metabolite production

% A possibly "absurd" result of the modelling process: Glucose 6-phosphate "disappears" to biomass
% before it can be used for glycolysis. 

%Issue: should the drain to biomass be taken into account (substracted) when
%calculating split ratios?


%R2 (Flux into Entner-Doudoroff pathway):

%Focus Metabolite: 6PG
R(2,1) = fluxvector(69)/fluxvector(29);
R(2,2) = R(2,2); %6PG is not drained in the biomass reaction.

%R2 = edd/pgl
%Both edd and pgl are irreversible, so no need to use subplus.

%R3 (Flux into methylgoyoxal pathway):
%Focus metabolite: DHP (DHAP)
R(3,1) = fluxvector(71)/(subplus(fluxvector(15))+subplus(-fluxvector(16)));
R(3,2) = R(3,1); %DHP is not drained in the biomass reaction.

%R3 = mgsA/(fbaAB+TpiA)
%Relevant reactions: mgsA (71), fbaAB (15), TpiA (16)


%R4 (PEP to pyruvate flux):
%Focus metabolite: PEP 
R(4,1) = (fluxvector(21)-fluxvector(27)+fluxvector(72))/(subplus(fluxvector(20))+fluxvector(25));
R(4,2) = (fluxvector(21)-fluxvector(27)+fluxvector(72))/(subplus(fluxvector(20))+fluxvector(25)-0.77*fluxvector(85)); %PEP is drained with stoichiometric coefficient 0.77 in the biomass reaction

%Definition in Table 1, Schuetz. et al. 2007:
%R4 = (pykAF (#21)+pts)/(eno+pps+pckA)

%Equation used here.

%R4 = (pykAF-pps+pts)/(eno+pck)
%See notes about R4 in computesplitratiosEXP.m 

%R5 (Pyruvate to acetyl-CoA):

R(5,1) = (fluxvector(22)+fluxvector(60))/(fluxvector(21)+fluxvector(23)+fluxvector(24)+fluxvector(67)+subplus(fluxvector(68))+fluxvector(70)+fluxvector(72)+subplus(-fluxvector(60)));
R(5,2) = (fluxvector(22)+fluxvector(60))/(fluxvector(21)+fluxvector(23)+fluxvector(24)+fluxvector(67)+subplus(fluxvector(68))+fluxvector(70)+fluxvector(72)+subplus(-fluxvector(60))-2.94*fluxvector(85));
%Pyruvate is drained in the biomass reaction with stoichiometric
%coefficient 2.94

%Definition in Table 1, Schuetz et al. 2007:

%R5 = (AceE + pflB + tdcE)/(PykA+pykF+MaeA+MaeB+Dld+ldhA+Eda+pts+pflB+tdcE)

%Consider if the use of Dld+ldA should be somehow changed to give a closer
%mapping to their corresponding experimental reaction.

%Equation used here:

%R5 = (Ace + pflB/tdcE)(pykAF+MaeA+MaeB+Dld+subplus(ldhA)+Eda+pts+pflB/tdcE)

%R6 (Flux into TCA cycle):
%Metabolite: Acetyl-CoA
R(6,1) = fluxvector(36)/(fluxvector(22)+subplus(fluxvector(60))+fluxvector(66)+subplus(-fluxvector(62)));
R(6,2) = fluxvector(36)/(fluxvector(22)+subplus(fluxvector(60))+fluxvector(66)+subplus(-fluxvector(62))-2.41*fluxvector(85));
%Acetyl-CoA is drained in the biomass reaction with stoichiometric
%coefficient 2.41


%Definition in Table 1, Schuetz et al. 2007:

%R6 = (gltA+prpC)/(AceE+pflB+TdcE+Acs+AdhE+MhpF)

%Equation used here: 

%R6 = [gltA/prpC]/(acEF+[pflB/tdcE]+Acs+subplus(-adhE/mhpF))

%Relevant reactions:

%gltA/prpC, Citrate synthase (#36)
%AceEF, pyruvate dehydrogenase (#22)
%pflB/tdcE, pyruvate formate lyase (#60)
%Acs, Acetyl-CoA Synthase (#66)
%adhE/mhpF, Acetaldehyde dehydrogenase (#62)

%R7 (Flux into glyoxylate shunt)
%Focus metabolite: Isocitrate
%Definition in Table 1, Schuetz et al. 2007:

%R7 = AceA/(acnA+acnB)

%Equation used here:

%R7 = aceA/acnAB

R(7,1) = fluxvector(46)/fluxvector(37);
R(7,2) = R(7,1);


%R8 (Oxaloacetate to PEP flux)
%Focus metabolite: Oxaloacetate
%Definition in Table 1, Schuetz et al. 2007:

%R8 = pckA/(Mdh+Mqo+Ppc)

%Equation used here:

%R8 = pckA/(Mdh+Mqo+Ppc)

R(8,1) = fluxvector(25)/(subplus(fluxvector(44))+subplus(fluxvector(45))+fluxvector(48));
R(8,2) = fluxvector(25)/(subplus(fluxvector(44))+subplus(fluxvector(45))+fluxvector(48)-1.65*fluxvector(85));
%Oxalocatete is drained in biomass reaction with stoichiometric coefficient
%1.65

%Relevant reactions

%pck (#25):OAA + ATP -> PEP + CO2
%mdh (#44): MAL <-> NADH + OA
%mqo (#45): MAL + Q -> QH2 + OA
%Ppc (#48), phosphoenolpyruvate carboxylase: PEP + CO2 -> OA 

%Note that pck and Ppc are reverse reactions, except that pck requires ATP.
%Thus running both reactions at once represent a futile cycle. When
%optimizing cell behavior in the model, we can expect pck to be zero.



%R9 (Acetate secretion):
%Focus metabolite: Acetyl-CoA
%Definition in Table 1, Schuetz et al. 2007:

%R9 = pta/(aceE+pflB+tdcE+Acs+AdhE+MhpF)

R(9,1) = fluxvector(64)/(fluxvector(22)+subplus(fluxvector(60))+fluxvector(66)+subplus(-fluxvector(62)));
R(9,1) = fluxvector(64)/(fluxvector(22)+subplus(fluxvector(60))+fluxvector(66)+subplus(-fluxvector(62))-2.41*fluxvector(85));
%Acetyl-CoA is drained in the biomass reaction with stoichiometric
%coefficient 2.41


%R10 (Ethanol secretion)

%Definition in Table 1, Schuetz et al. 2007:

%R10 = (adhE+Mhpf)/(aceE+pflB+tdcE+Acs+AdhE+MhpF)

R(10,1) = fluxvector(62)/(fluxvector(22)+subplus(fluxvector(60))+fluxvector(66)+subplus(-fluxvector(62)));
R(10,2) =fluxvector(62)/(fluxvector(22)+subplus(fluxvector(60))+fluxvector(66)+subplus(-fluxvector(62))-2.41*fluxvector(85));

%If Split ratio is not defined due to mathematical singularity, set the
%split ratio to zero:
for i = 1:size(R,1)
    for j = 1:size(R,2)
    if isnan(R(i,j)) == 1
        R(i,j) = 0;
    else
    end
    end
end

    case 2
        R(1,1) = fluxvector(74)/(fluxvector(50)+fluxvector(48)+fluxvector(74));
R(1,2) = fluxvector(74)/(fluxvector(50)+fluxvector(48)+fluxvector(74)-0.205*fluxvector(13));

%Glucose uptake is only by pts pathway in the E. coli core model. An
%equivalent of the glk (glucokinase) reaction in Schuetz. model does not
%exist. There is thus one less flux in the denominator. 


%The reaction "PGI" (#74) in ECM is the equivalent of "Pgi" in the Schuetz
%model.

%The reaction "GLCpts" (#50) in ECM is the equivalent of "Pts" in the
%Schuetz model.

%reaction "R_G6PDH2r" (#48) in the ECM is the equivalent of R_zwf in the Schuetz model.


R(2,1) = fluxvector(98)/fluxvector(76);
R(2,2) = R(2,1);

%Reaction "R_PGL" (#76) in ECM is the equivalent of "R_PGL" (#34) in the
%Schuetz model. The species D6PGL in Schuetz model is called 6PGL in ECM.


R(3,1) = fluxvector(100)/(fluxvector(40)+subplus(-fluxvector(95)));
R(3,2) = R(3,1);
%mgsA (#100)
%FBA (#40)



%R4 (PEP to pyruvate flux):
%Focus metabolite: PEP 
R(4,1) = (fluxvector(83)+fluxvector(74))/(fluxvector(18)+fluxvector(81)+fluxvector(80));
R(4,2) = (fluxvector(83)+fluxvector(74))/(fluxvector(18)+fluxvector(81)+fluxvector(80)-0.5191*fluxvector(13));

%The ECM reaction "PYK" (#83) is the equivalent of the Schuetz reactions
%"PykA" and "PykF", which are identical.

%The ECM reaction "ENO" (#18) is the equivalent of the Schuetz reaction
%"eno" (#23). 

%The ECM reaction "PPS" (#81) is the equivalent of the Schuetz reaction
%"pps" (#32).

%The ECM reaction "PPCK" (#80) is the equivalent of reaction "pck" (#29) in the Schuetz model. 


%R5 (Pyruvate to acetyl-CoA):
R(5,1) = (fluxvector(71)+fluxvector(73))/(fluxvector(83)+fluxvector(65)+fluxvector(66)+fluxvector(61)+fluxvector(74)+fluxvector(73));
R(5,2) = (fluxvector(71)+fluxvector(73))/(fluxvector(83)+fluxvector(65)+fluxvector(66)+fluxvector(61)+fluxvector(74)+fluxvector(73)-2.8328*fluxvector(13));

%The ECM reaction "R_PDH" (pyruvate dehydrogenase) (#71) is the equivalent
%of the reaction "R_AceEF" in the Schuetz model.

%The ECM reaction "R_PFL"  (#73) is the equivalent of reaction "pflB", which is identical to reaction "tdcE", in
%the Schuetz model.

%The ECM reactions "ME1" (#65) and "ME2" (#66) are equivalent to the
%reactions "maeA" and "maeB" in the schuetz model.

%The ECM reaction "LDH_D" (#61) is the equivalent of reaction "dld" in the
%Schuetz model (different cofactors). ECM does not have a reaction for
%l-lactate dehydrogenase.

%The reaction "eda" in the Schuetz model does not appear to be modelled in
%the ECM.

%R6 (Flux into TCA cycle):
%Metabolite: Acetyl-CoA
R(6,1) = fluxvector(15)/(fluxvector(71)+fluxvector(73)+fluxvector(1));
R(6,2) = fluxvector(15)/(fluxvector(71)+fluxvector(73)+fluxvector(1)-3.7478*fluxvector(13));

%The ECM reaction "CS" (Citrate synthase) (#15)is equivalent to the reaction "GltA", which is identical to reaction "prpC", in the
%Schuetz model. 

%The reaction "R_acs" in the Schuetz model (acetyl-CoA synthase) does not
%appear to be present in the ECM model

%The ECM reaction "R_ACALD" (#1) is equivalent to the reaction "R_adhE", which is identical to "R_mhpF", in the
%Schuetz. model


%R7 (Flux into glyoxylate shunt)
%Focus metabolite: Isocitrate
R(7,1) = fluxvector(60)/fluxvector(4);
R(7,2) = R(7,1);

%The ECM reaction "ICL" (#60) (isocitrate lyase) is the equivalent of reaction "aceA" (#61) in
%the Schuetz model.


%The ECM reaction "R_ACONTa" (#4) is the equivalent of reaction "acnA". which is identical to reaction "acnB", in
%the Schuetz model.



%R8 (Oxaloacetate to PEP flux)
%Focus metabolite: Oxaloacetate
R(8,1) = fluxvector(80)/(fluxvector(64)+fluxvector(79));
R(8,2) = fluxvector(80)/(fluxvector(64)+fluxvector(79)-1.7867*fluxvector(13));

%The ECM reaction "MDH" (#64) is equivalent to the reaction "mdh" in the
%schuetz model.

%The reaction "R_mqo" in the Schuetz model does not appear to exist in the
%ECM.

%The ECM reaction "R_PPC" (PEP carboxylase) (#79) is equivalent to the reaction "R_ppc" in
%the schuetz model


%R9 (Acetate secretion):
%Focus metabolite: Acetyl-CoA
R(9,1) = fluxvector(82)/(fluxvector(71)+fluxvector(73)+fluxvector(1));
R(9,2) = fluxvector(82)/(fluxvector(71)+fluxvector(73)+fluxvector(1)-3.7478*fluxvector(13));
%The ECM reaction "R_PTAr" (phosphotransacetylase) (#82) is the equivalent of the reaction "pta" (phosphate acetyltransferase) in the Schuetz model. 

%The reaction "R_acs" in the Schuetz model (acetyl-CoA synthase) does not
%appear to be present in the ECM model


R(10,1) = fluxvector(1)/(fluxvector(71)+fluxvector(73)+fluxvector(1));

    case 3


%R1 (Flux into glycolysis)
%Focus metabolite: G6P

R(1,1) = fluxvector(2077)/(fluxvector(1500)+fluxvector(1353)+subplus(fluxvector(1283))+subplus(-fluxvector(2077)));
R(1,2) = R(1,1); %G6P is not drained in biomass reaction

%Definition in Table 1, Schuetz et al. 2007:

%R1 = pgi/(glk+pts+zwf+pgi)


%Reaction "pgi" (#2077) of iJO1366 ~ reaction "pgi" of Schuetz model

%Reaction "HEX1" (#1500) of iJO1366 ~ reaction "glk" of Schuetz model

%Reaction "R_GLCptspp" (#1353) of iJO1366 ~ reaction "ptsGHI" of Schuetz
%model

%reaction "R_G6PDH2r" (#1283) of iJO1366 ~ R_zwf in the Schuetz model.
%(reversible)

%reaction "PGI" (#2077) of iJO1366 ~ "pgi" in the Schuetz model


%R2 (Flux into Entner-Doudoroff pathway):

%Focus Metabolite: 6PG

R(2,1) = fluxvector(1095)/fluxvector(2079);
R(2,2) = R(2,1); %6PG is not drained in biomass reaction

%Reaction "R_EDD" (#1095) of iJO1366 ~ reaction "edd" in the Schuetz model

%Reaction "R_PGL" (#2079) of iJO1366 ~ reaction "pgl" in the Schuetz model


%R3 (Flux into methylgoyoxal pathway):
%Focus metabolite: DHP (DHAP)

R(3,1) = fluxvector(1790)/(subplus(fluxvector(1151))+subplus(-fluxvector(2456)));
R(3,2) = R(3,1); %DHAP is not drained in biomass reaction

%Reaction "MGSA" ´(#1790) of iJO1366 ~ reaction "mgsA" in the Schuetz model

%Reaction "R_FBA" (#1151) of iJO1366 ~ reactions "fbaA"/"fbaB" in the Schuetz
%model.

%Reaction "R_TPI" (#2456) of iJO1366 = reaction "R_tpiA" in the Schuetz
%model.


%R4 (PEP to pyruvate flux):
%Focus metabolite: PEP 

R(4,1) = (fluxvector(2266)+fluxvector(1353)-fluxvector(2198))/(fluxvector(1102)+fluxvector(2183));
R(4,2) = R(4,1); %PEP is not drained in biomass reaction 

%Definition in Table 1, Schuetz. et al. 2007:
%R4 = (pykAF (#21)+pts)/(eno+pps+pckA)

%Equation used here:

%R4 = (pykAF-pps+pts)/(eno+pck)


%Reaction "R_PYK" (2266) of iJO1366 ~ reactions "pykA"/"pykF" in the
%Schuetz model.

%Reaction "R_GLCptspp" (#1353) of iJO1366 ~ reaction "ptsGHI" of Schuetz
%model

%Reaction "R_ENO" (1102) of iJO1366 ~ reaction "Eno" of Schuetz model

%Reaction "R_PPS" (2198) of iJO1366 ~ reaction "pps" of Schuetz model

%Reaction "R_PPCK" (2183) of iJO1355 ~ reaction "pck" of Schuetz model

%See notes on R4 in computesplitratiosEXP.m 

%R5 (Pyruvate to acetyl-CoA):
%Focus metabolite: Pyruvate

R(5,1) = (fluxvector(2047)+fluxvector(2067))/(fluxvector(2266)+fluxvector(1761)+fluxvector(1762)+subplus(fluxvector(1622))+fluxvector(1623)+fluxvector(1601)+fluxvector(1602)+fluxvector(1094)+fluxvector(1353)+fluxvector(2067));
R(5,2) = R(5,1); %Pyruvate is not drained in biomass reaction

%Definition in Table 1, Schuetz et al. 2007:

%R5 = (AceE + pflB + tdcE)/(PykA+pykF+MaeA+MaeB+Dld+ldhA+Eda+pts+pflB+tdcE)

%Reaction "R_PDH" (2047) (irreversible) of  iJO1366 ~ reaction "AceEF" in the Schuetz model

%Reaction "R_PFL" (2067) (irreversible) of iJO1366 ~ reactions "pflB"/"tdcE" in the
%Schuetz model

%Reaction "R_PYK" (2266) of iJO1366 ~ reactions "pykA"/"pykF" in the
%Schuetz model.

%Reactions "R_ME1" (1761) and "R_ME2" (1762) (both irreversible) of iJO1366 ~ reactions "maeA"
%and "maeB" of Schuetz model, respectively. 

%CONVERSION OF LACTATE TO PYRUVATE:

%Reactions in iJO1366:

%"R_LDH_D" (#1622): D-lactate+NAD+ -> NADH + pyruvate (E.C 1.1.1.28)
%(reversible)
%"R_LDH_D2" (#1623): D-lactate + Q8 -> pyruvate + Q8H2 (E.C 1.1.2.4)
%(irreversible)
%"R_L_DASH_LACD2" (#1601) : L-lactate + Q8 -> pyruvate + Q8H2 (E.C 1.1.2.3)
%(irreversible)
%"R_L_DASH_LACD3" (#1602): L-lactate + mqn8 -> pyruvate + mql8 (E.C
%1.1.2.3) (irreversible)

%Reactions in Schuetz model:

%dld: lactate -> FADH + pyruvate
%ldhA: lactate -> pyruvate + NADH

%Schuetz model does not discern between D- and L-lactate.

%The closest equivalents are reaction "R_LDH_D" in iJO1366 and "ldhA" in
%Schuetz modell

%CONCLUSION: "R_LDH_D", "R_LDH_D2", "R_L_DASH_LACD" and "R_L_DASH_LACD3"
%are used in place of "dld" and "ldhA" in calculating R5.


%reaction "R_EDA" (#1094) (irreversible) of iJO1366 = "R_eda" of Schuetz model
%Reaction "R_GLCptspp" (#1353) (irreversible) of iJO1366 ~ reaction "ptsGHI" of Schuetz
%model
%Reaction "R_PFL" (2067) of iJO1366 ~ reactions "pflB"/"tdcE" in the
%Schuetz model. However, in iJO1366 it is irreversible. As such, its inclusion as a
%producing flux of pyruvate is superflous here.


%R6 (Flux into TCA cycle):
%Metabolite: Acetyl-CoA
R(6,1) = fluxvector(878)/(fluxvector(2047)+fluxvector(2067)+fluxvector(547)+subplus(fluxvector(493)));
R(6,2) = fluxvector(878)/(fluxvector(2047)+fluxvector(2067)+fluxvector(547)+subplus(fluxvector(493))-0.000279*fluxvector(7));
%Acetyl-CoA is drained with stoichiometric coefficient 0.000279 in WT
%biomass reaction

%Definition in Table 1, Schuetz et al. 2007:

%R6 = (gltA+prpC)/(AceE+pflB+TdcE+Acs+AdhE+MhpF)

%Equation used here: 

%R6 = CS/PDH+pflB+Acs+subplus(-acald))

%The iJO1366 reaction "CS" (Citrate synthase) (#878) (irreversible) is equivalent to the reaction "GltA", which is identical to reaction "prpC", in the
%Schuetz model. 

%Reaction "R_PDH" (2047) (irreversible) of  iJO1366 ~ reaction "AceEF" in the Schuetz model

%Reaction "R_PFL" (2067) (irreversible) of iJO1366 ~ reactions "pflB"/"tdcE" in the
%Schuetz model

%Reaction "R_ACS" (#547) (irreversible) of iJO1366 ~ reaction "R_acs" of the Schuetz model.

%Reaction "R_ACALD" (#493) (reversible) of iJO1366 = reaction "adhE"/"mhpF" in the Schuetz
%model

%R7 (Flux into glyoxylate shunt):
%Focus metabolite: Isocitrate
%Definition in Table 1, Schuetz et al. 2007:

%R7 = AceA/(acnA+acnB)

R(7,1) = fluxvector(1552)/(fluxvector(536));
R(7,2) = R(7,1); %Isocitrate is not drained in biomass reaction

%Reaction "R_ICL" (#1552) (irreversible) of iJO1366 ~ reaction "aceA" of Schuetz model.

%Reaction "R_ACONTa" (#536) of iJO1366 ~ reaction "acnAB" of Schuetz model.


%R8 (Oxaloacetate to PEP flux):
%Focus metabolite: Oxcaloacetate
%Definition in Table 1, Schuetz et al. 2007:

%R8 = pckA/(Mdh+Mqo+Ppc)

%Equation used here:

%R8 = pckA/(Mdh+Mqo+Ppc)

R(8,1) = fluxvector(2183)/(subplus(fluxvector(1758))+fluxvector(1759)+fluxvector(1760)+fluxvector(2181));
R(8,2) = R(8,1);
%Reaction "R_PPCK" (#2183) (irreversible) of iJO1366 ~ reaction "pck" of Schuetz model

%Reaction "R_MDH (#1758) (reversible) of iJO1366 ~ reaction "mdh" of Schuetz model

%Reactions "R_MDH2" (#1759) and "R_MDH3" (#1760) (both irreversible) ~ reaction "R_mqo" of
%Schuetz model

%Reaction "R_PPC" (#2181) (irreversible) of iJO1366 ~ reaction "ppc" of Schuetz model.



%R9 (Acetate secretion):
%Focus metabolite: Acetyl-CoA

%Definition in Table 1, Schuetz et al. 2007:

%R9 = pta/(aceE+pflB+tdcE+Acs+AdhE+MhpF)


R(9,1) = fluxvector(2240)/(fluxvector(2047)+fluxvector(2067)+fluxvector(547)+subplus(fluxvector(493)));
R(9,2) = fluxvector(2240)/(fluxvector(2047)+fluxvector(2067)+fluxvector(547)+subplus(fluxvector(493))-0.000279*fluxvector(7));

%Reaction "R_PTAr" (#2240) (reversible) of iJO1366 ~ reaction "pta" of Schuetz model. 
%Reaction "R_PDH" (2047) (irreversible) of  iJO1366 ~ reaction "AceEF" in the Schuetz model
%Reaction "R_PFL" (2067) (irreversible) of iJO1366 ~ reactions "pflB"/"tdcE" in the
%Schuetz model
%Reaction "R_ACS" (#547) (irreversible) of iJO1366 ~ reaction "R_acs" of the Schuetz model.
%Reaction "R_ACALD" (#493) (reversible) of iJO1366 = reaction "adhE"/"mhpF" in the Schuetz
%model



%R10 (Ethanol secretion):
%Focus metabolite: Acetyl-CoA


%Definition in Table 1, Schuetz et al. 2007:

%R10 = (adhE+Mhpf)/(aceE+pflB+tdcE+Acs+AdhE+MhpF)


R(10,1) = fluxvector(493)/(fluxvector(2047)+fluxvector(2067)+fluxvector(547)+subplus(fluxvector(493)));
R(10,2) = fluxvector(493)/(fluxvector(2047)+fluxvector(2067)+fluxvector(547)+subplus(fluxvector(493))-0.000279*fluxvector(7));
%Reaction "R_ACALD" (#493) (reversible) of iJO1366 = reaction "adhE"/"mhpF" in the Schuetz
%model
end

%If Split ratio is not defined due to mathematical singularity, set the
%split ratio to zero:
for i = 1:size(R,1)
    for j = 1:size(R,2)
    if isnan(R(i,j)) == 1
        R(i,j) = 0;
    else
    end
    end
end

splitratios = R;

end