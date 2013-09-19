function [Aeb,Aec,Agb,Agc,Alb,Alc,be,bg,bl,cb,cc,lbb,lbc,ubb,ubc] = lindo2gams(A,b,c,lb,ub,vartype,csense)
%LINDO2GAMS - Convert Lindo inputs to GAMS compatible
%matrices. Primarily splits the different types of constraints into
%separate matrices
%
% [Aeb,Aec,Agb,Agc,Alb,Alc,be,bg,bl,lbb,lbc,ubb,ubc] = lindo2gams(A,b,lb,ub,vartype,csense)
%
% Markus Herrgard

sel_bin = (vartype == 'B');
sel_cont = (vartype == 'C');

sel_eq = (csense == 'E');
sel_ge = (csense == 'G');
sel_le = (csense == 'L');

Aeb = A(sel_eq,sel_bin);
Aec = A(sel_eq,sel_cont);
Agb = A(sel_ge,sel_bin);
Agc = A(sel_ge,sel_cont);
Alb = A(sel_le,sel_bin);
Alc = A(sel_le,sel_cont);

be = b(sel_eq);
bg = b(sel_ge);
bl = b(sel_le);

cb = c(sel_bin);
cc = c(sel_cont);

lbb = lb(sel_bin);
lbc = lb(sel_cont);
ubb = ub(sel_bin);
ubc = ub(sel_cont);
