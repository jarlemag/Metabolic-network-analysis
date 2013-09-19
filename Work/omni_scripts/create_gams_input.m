function create_gams_input(A,b,c,lb,ub,vartype,csense,file,intvarnames,contvarnames)
%CREATE_GAMS_INPUT Create GAMS input from LINDO input matrices
%
% create_gams_input(A,b,c,lb,ub,vartype,csense,file,intvarnames,contvarnames)
%
% A         LHS matrix
% b         RHS vector
% lb        Lower bounds
% ub        Upper bounds
% vartype   Variable types
% csense    Constraint sense
% file      File name
% intvarnames   Variable names for int variables
% contvarnames  Variable names for cont variables
%
% Markus Herrgard 5/1/03

[Aeb,Aec,Agb,Agc,Alb,Alc,be,bg,bl,cb,cc,lbb,lbc,ubb,ubc] = lindo2gams(A,b,c,lb,ub,vartype,csense);
print_gams_input(file,Aeb,Aec,Agb,Agc,Alb,Alc,be,bg,bl,lbb,lbc,ubb,ubc,cb,cc,intvarnames,contvarnames);

nb = length(lbb);
nc = length(lbc);
ne = length(be);
ng = length(bg);
nl = length(bl);
lb_r = [lbb;lbc];
ub_r = [ubb;ubc];
b_r = [be;bg;bl];
A_r = [[Aeb Aec];[Agb Agc];[Alb Alc]];
c_r = [cb;cc];
%save gams_osi_tmp.mat A_r Aeb Aec Agb Agc Alb Alc be bg bl lb_r ub_r b_r c_r nb nc ne ng nl intvarnames contvarnames
