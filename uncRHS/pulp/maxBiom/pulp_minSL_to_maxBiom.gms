************************* Projection using Linear Programming (PULP) ****************
*       Author: Hoang Dinh
* Minimizng the non-zero RHS terms subject to the improved growth yield 
*************************************************************************************

$INLINECOM /*  */
$set path ../GAMS/
$set nscale 1E3

options
        limrow = 10000
        limcol = 10000
        optCR = 1E-9
        optCA = 0.0
        iterlim = 100000
        decimals = 8
        reslim = 100000
        work = 50000000
        sysout = off
        solprint = on;
        
Sets

i
$include "%path%metabolites.txt"
metint(i)
$include "%path%metabolites_int.txt"
metext(i)
$include "%path%metabolites_ext.txt"
metzero(i)
$include "%path%metabolites_zero.txt"
j
$include "%path%reactions.txt"
e
$include "%path%elements.txt"
;

Parameters
S(i,j)
$include "%path%sij.txt"
Ne(i,e)
$include "%path%elements_comp.txt"
LB(j)
$include "%path%lb.txt"
UB(j)
$include "%path%ub.txt"
MW(i)
$include "%path%MW.txt"
maxSL(i)
$include "maxSL.txt"
;

Variables
z, v(j), slL(i), slU(i)
;

****************** SETTING LOWER AND UPPER BOUNDS FOR VARIABLES ******************
v.lo(j) = LB(j)*%nscale%;
v.up(j) = UB(j)*%nscale%;

v.fx('R__BIOMASS_Ec_iML1515_WT_75p37M') = 0;
v.fx('R__BIOMASS_Ec_iML1515_core_75p37M') = 924.01715914 - 1e-4;
v.lo('R__EX_glc__D_e') = -10*%nscale%; v.up('R__EX_glc__D_e') = 1000*%nscale%;
v.lo('R__EX_o2_e') = -1000*%nscale%; v.up('R__EX_o2_e') = 1000*%nscale%;
v.lo('R__EX_nh4_e') = -1000*%nscale%; v.up('R__EX_nh4_e') = 1000*%nscale%;
v.fx('R__ATPM') = 6.86*%nscale%;

slL.lo(i) = 0; slL.up(i) = maxSL(i)*%nscale%;
slU.lo(i) = 0; slU.up(i) = maxSL(i)*%nscale%;

* Disable imbalance of extracellular metab
slL.fx(i)$(metext(i)) = 0;
slU.fx(i)$(metext(i)) = 0;

***************************** EQUATION DEFINITIONS **************************************
Equations
Obj, Stoic, Elem
;

Obj..			z =e= sum(i, MW(i)*(slL(i) + slU(i)));
Stoic(i)..		sum(j, S(i,j)*v(j)) =e= slU(i) - slL(i);
Elem(e)..		sum(i, Ne(i,e)*(slU(i) - slL(i))) =e= 0;


************************ DECLARING THE MODEL with constraints****************************
Model minSL
/Obj, Stoic, Elem/;
minSL.optfile = 1;

******************* SOLVE *******************
Solve minSL using lp minimizing z;

******************* WRITE TO TXT FILE *******************
file ff /pulp_minSL_to_maxBiom.sl_vals.txt/;
put ff;
loop(i$metint(i),
	if ( (slL.l(i) gt 1E-9),
		put i.tl:0, system.tab, 'slL', system.tab, slL.l(i):0:8/;
);
);

loop(i$metint(i),
	if ( (slU.l(i) gt 1E-9),
		put i.tl:0, system.tab, 'slU', system.tab, slU.l(i):0:8/;
);
);
putclose ff;

file ff2 /pulp_minSL_to_maxBiom.modelStat.txt/;
put ff2;
put minSL.modelStat/;
putclose ff2;

file ff3 /pulp_minSL_to_maxBiom.flux.txt/;
put ff3;
loop(j,
	put j.tl:0, system.tab, v.l(j):0:8/;
);
putclose ff3;
