* This is a driver for OMNI optimization problems *
$offdigit
$offlisting
$onempty

$Include 'omni_model.gms';

Variables

         x(jc)   cont solution
         y(jb)   bin solution
         z       Objective value;

Binary variable y(jb);

* Lower and upper bounds *

x.lo(jc) = lbc(jc);
x.up(jc) = ubc(jc);
y.lo(jb) = lbb(jb);
y.up(jb) = ubb(jb);

Equations
         Obj
         eqs(ie)
         ges(ig)
         les(il);

         Obj..                 z =e= sum(jb,cb(jb)*y(jb)) + sum(jc,cc(jc)*x(jc));
         eqs(ie)..      sum(jb,Aeb(ie,jb)*y(jb)) + sum(jc,Aec(ie,jc)*x(jc)) =e= be(ie);
         ges(ig)..      sum(jb,Agb(ig,jb)*y(jb)) + sum(jc,Agc(ig,jc)*x(jc)) =g= bg(ig);
         les(il)..      sum(jb,Alb(il,jb)*y(jb)) + sum(jc,Alc(il,jc)*x(jc)) =l= bl(il);

Model OMNI /all/;

option mip = CPLEX;

option reslim = 30000;
option optcr = 0.0;

solve OMNI using mip minimizing z;

display x.l, y.l;

file stat /omni.stat/
put stat;
put 'model status ';
put OMNI.modelstat /;
put 'solver status ';
put OMNI.solvestat /;
put 'objective ';
put z.l /;

file sol_kos /omni_kos.out/;
put sol_kos;
loop(jb, 
	if (y.l(jb) < 1,
		put @1, jb.tl:20, @20, y.l(jb):16:12/;
	);
);
put /;

file sol_flux /omni_fluxes.out/;
put sol_flux
loop(jc, 
	if (x.l(jc) > 0,
		put @1, jc.tl:20, @20, x.l(jc):16:12/;
	);
);
put /;
