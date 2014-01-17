function [numb res names] = HighFluxBackbone(model,flux_vector)

% result = HFB(model,flux_vector)
%
% Function calculates the High-Flux Backbone.  Note that if two or more 
%              reactions have the same (max) rate for a metabolite, all are
%              potentially included in the HFB. This function uses the
%              results of the COBRA Toolbox.
%
% Please cite the following paper:
%      E. Almaas, B. Kovacs, T. Vicsek, Z.N. Oltvai and A.-L. Barabási. 
%      “Global organization of metabolic fluxes in the bacterium 
%       Escherichia coli.” Nature 427, 839 (2004).
%
%VARIABLES:
% model        contains the COBRA Toolbox formatted constrain-based model
% flux_vector  contains the COBRA Toolbox formatted vector of fluxes from 
%              which the HFB will be extracted
%
% numb         returns the number of reactions participating in the HFB
% res          is a vector of all reactions, with 0 (no) or 1 (yes)
%              denoting if that reaction was part of the HFB
% names        is a vector of reaction names of only the HFB participants
%
%
%EXAMPLE USE:
%  data = optimizeCbModel(model,'max');
%
%  [numb res names] = HFB(model,data.x)
%

mets = size(model.mets,1);
rxns = size(model.rxns,1);

prod_cand = zeros(rxns,1);
pp        = 0; #Number of producing candidates found
cons_cand = zeros(rxns,1);
pc        = 0; #Number of consuming candidates found

%First, find all largest consuming and producing reactions
for i=1:mets
    j    = find(model.S(i,:)); %# Find the indices of all non-zero elements in the stoichoimetric matrix row corresponding to the current metabolite.
    disp('j') %#Corresponds to the reaction #s of all reactions that the metabolite is part of.
    disp(j)
    n    = size(j,2); %# The number of reactions that the current metabolite partakes in.
    disp(n)
    bigp = +1e-7;   %Numerical cut-off for fluxes to be considered nonzero/largest flux found so far for this metabolite
    bign = -1e-7; %Numerical cut-off/smallest flux found so far for this metabolite
    indp = zeros(1,200);
    indn = zeros(1,200);
    fp   = 0;
    fn   = 0;
    if n > 0
        %Loop over all the reactions metabolite i participates in.
        for k=1:n %#For all reactions the metabolite participates in
            value = model.S(i,j(k))*flux_vector(j(k)); #The flux of the metabolite in the current reaction
            if (value >= bigp) %If the flux is the largest found yet for this metabolites
                %Metabolite is a product
                if (value == bigp)
                    fp       = fp + 1; #Increase the number of largest (tied) producing reactions by one
                    indp(fp) = j(k); #Record the number of (tied for) largest producing reactions found in the vector indp. The reaction # j(k) is assigned to
                    		     #position fp. Thus, the reaction number of the first reaction to be found is saved in indp(1). If additional reactions with the
                    		     #same flux is found, their reaction #s are stored in position 2, 3, etc. up to a maximum of 200 tied reactions.
                else
                    fp       = 1; #Set the number of largest producing reactions to one (no tie-breaks)
                    bigp     = value; #Set the new largest value 
                    indp(fp) = j(k); #Save the reaction id of the largest value in indp(1).
                end
            elseif (value < bign) %If the flux is the most negative found yet for this metabolite
                %Metabolite is a substrate
                if (value == bign)
                    fn       = fn + 1 #Increase the number of largest consuming reactions by one;
                    indn(fn) = j(k);
                else
                    fn       = 1; #Set the number of largest producing reactions to one
                    bign     = value;
                    indn(fn) = j(k);
                end
            end
        end
        
        %Take care of any producing-reaction candidates identified
        if (fp > 0) #If any largest-producing reactions were found (true if at least one flux is producing)
            for l=1:fp #For every largest-producing reaction (several reactions can be tied, giving more than one largest-producing reaction)
                if (ismember(indp(l),prod_cand)) #If the candidate reaction # is recorded in the "producing candidates" list, do nothing.
                else #If it is not already a member
                    pp = pp + 1; #Increment the number of candidate reactions identified
                    prod_cand(pp) = indp(l); #Add the reaction # of the candidate reaction to the list of candidate reactions.
                end
            end
        end
        
        %Take care of any consuming-reaction candidates identified
        if (fn > 0)
            for l=1:fn
                if (ismember(indn(l),cons_cand))
                else
                    pc = pc + 1;
                    cons_cand(pc) = indn(l);
                end
            end
        end
    end
end

res = zeros(rxns,1);

%Consolidate the two candidate lists: for a reaction to be member of the
%   HFB, it must be a maximal production and consumption reaction, not
%   necessarily for same metabolite.

for i=1:pp
    if (ismember(prod_cand(i),cons_cand))
        res(prod_cand(i)) = 1;
    end
end

list  = find(res);
numb  = size(list,1);
names = model.rxns(list);
