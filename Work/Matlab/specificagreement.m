
%Compute specific agreement between a set of computed and experimental flux split ratios
%See Schuetz et al. 2007: "Systematic evaluation of objective functions..."

%Syntax: rho = specificagreement(R_comp_max,R_comp_min,R_exp,sigma)

function rho = specificagreement(R_comp_max,R_comp_min,R_exp,sigma)


if length(R_comp_max)~=length(R_comp_min)
    error('Error: All flux vectors must be same length')
elseif length(R_comp_min)~=length(R_exp)
    error('Error: All flux vectors must be same length')
else n = length(R_comp_max);
end


delta_R_comp = R_comp_max - R_comp_min;
%Range of computed values
%disp('Range of computed values:')
%disp(delta_R_comp)


R_comp_mean = (R_comp_max+R_comp_min)/2;
%disp('Mean computed split ratios')
%disp(R_comp_mean)
%See supplementary figure 4


rho = zeros(1,n);
for i = 1:n
    if delta_R_comp(i) == 0
    rho(i) = abs(R_exp(i)-R_comp_mean(i))*(0.05/sigma(i));
    else
    rho(i) = abs(R_exp(i)-R_comp_mean(i))*(delta_R_comp(i)/sigma(i));
    end
end
end