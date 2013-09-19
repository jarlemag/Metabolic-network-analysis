%compdist.m
%Computes a distance measure (euclidean or taxicab) (weighed or unweighed) between experimental and computed fluxes/split ratios. 
%Syntax:
%ds = compdist(fluxvector,options,sense)
%Output: Euclidean or taxicab norm.
%Distance can be computed on basis of raw fluxes or split ratios. If using
%split ratios, set options.compsplits = 1;
%options:
%model_id
%exp_id
%normchoice ('euclidean' or 'taxicab')

function ds = compdist(fluxvector,options,sense)

load('expdata')

if isfield(options,'verbflag') ==0
    options.verbflag =0;
end

verbflag = options.verbflag;

if nargin <3
    sense = 1;
    if verbflag == 1
    disp('compdist: Input variable "sense" not defined. Set to 1 (returns positive value of distance) by default')
    end
end

if isfield(options,'model_id') == 0
    options.model_id = 1;
end

if isfield(options,'exp_id') == 0
    options.exp_id = 1;
end

if isfield(options,'usescaledfluxes') == 0
    options.usescaledfluxes = 0;
end

if isfield(options,'useweights') == 0
    options.useweights = 0;
end

if isfield(options,'normchoice') == 0
    options.normchoice = 'euclidean';
end

if isfield(options,'excludereactions') == 1
    excludereactions=options.excludereactions;
else excludereactions = 0; %By default, do not exclude any reactions from being used  calculations.
end

if isfield(options,'debugmode') == 0
    options.debugmode = 0;
end

if isfield(options,'compsplits') == 0
    options.compsplits = 0;
end

if isfield(options,'substractBM') == 0
    options.substractBM = 0;
end

exp_id = options.exp_id;
usescaledfluxes = options.usescaledfluxes;
useweights = options.useweights;
debugmode = options.debugmode;
substractBM = options.substractBM;
compsplits = options.compsplits;

 if verbflag == 2
    disp('compdist.m: Excluding following experimental reactions:')
    disp(excludereactions)
 end

%Extract/calculate vector with which the experimental will be compared:
switch compsplits
    case 0
    Fcomp = extractflux(fluxvector,options); % Use raw fluxes
    if debugmode > 1
        disp('compdist.m: Using raw fluxes.')
    end
    
    case 1
         
         Fcomp = computesplits(fluxvector,options.model_id); %Use split ratios
        switch substractBM
            case 0
                if debugmode > 1
                disp('compdist.m: Using split ratios fluxes.')
                end
            Fcomp = Fcomp(:,1);
            case 1
                 if debugmode > 1
                disp('compdist.m: Using split ratios fluxes (BM substracted).')
                end
            Fcomp = Fcomp(:,2);
        end
end

switch exp_id
    %Select which experimental data to use:
    case 1
        %Fexp = Fexp_batch_aerobe;
        Fluxvalues = expdata.perrenoud.abs.batch.aerobe.fluxvalues;
        Fexp = Fluxvalues(:,1);
        %Ferrors = Fluxvalues(:,2); 
    case 2
        error('Data not available')
        
        
    case 3
        error('Data not available')
      
end

%If split ratios are used, replace experimental raw fluxes with the
%corresponding split ratios.
if compsplits == 1
    Fexp = computesplits(Fexp,0);
end
    

%Update this file to accept both absolute (Perrenoud 2005) and
%glucose-scaled (Schuetz 2007) data as options.


%Optional: Exclude some reactions:
if excludereactions ~= 0 %if any reactions to be excluded are specified
    for i = 1: length(excludereactions)
        Fexp(excludereactions(i)) = 0;
        Fcomp(excludereactions(i))= 0;
    end 
end

if options.debugmode == 1
    disp('Fexp:')
    disp(Fexp)
    disp('Fcomp:')
    disp(Fcomp)
end


if usescaledfluxes == 1;
scalingfactor = Fexp(1)/Fcomp(1);  %Scale to glucose flux
e = (Fexp/scalingfactor) - Fcomp;
else
    
    e = Fexp - Fcomp;
end

if length(Fexp) ~= length(Fcomp)
    error('compdist.m: Vectors not same length.')
else
    n = length(Fexp);
end

W = eye(n);

if useweights == 1
   sigma = Fluxvalues(:,2);  
   inverse_sigma =(1./sigma);
   sum_inverse_sigma = sum(inverse_sigma);
   for i=1:n
   W(i,i)=(1/sigma(i))*(1/sum_inverse_sigma);
   diagonal_sum =sum(diag(W));
   %Computes the sum of the diagonal of W. Should equal 1.
        if diagonal_sum ~=1
        disp(diagonal_sum)
        disp('Compdist.m: Warning: Sanity check failed. Sum of diagonal elements in W should equal 1.')
        else
        end
   end

end

if debugmode >0
    disp('W:');
    disp(W);
end

if debugmode > 0
    disp('compdist.m. normchoice:')
    disp(options.normchoice)
end


if   strcmp(options.normchoice,'euclidean')  == 1

    switch sense
        case 1
            ds = sqrt(e'*W*e);
        case -1
            ds = -sqrt(e'*W*e); %Report the negative value of the distance. For use with fmincon.
    end
elseif  strcmp(options.normchoice,'taxi') == 1
    
       switch sense 
           case 1
               ds = sum(abs(e.*diag(W)));
           case -1
               ds =-sum(abs(e.*diag(W)));
       end
end

if debugmode >1 
    disp('compdist.m: Returned value')
    disp(ds)
end


end