***************************** CMP-check ********************************
*       Author: Hoang Dinh
* Run linear programming optimization with slack variables to find coupled metabolites to the imbalance
************************************************************************

$INLINECOM /*  */
$set path ../GAMS/

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
j
$include "%path%reactions.txt"
;

Parameters
S(i,j)
$include "%path%sij.txt"
LB(j)
$include "%path%lb.txt"
UB(j)
$include "%path%ub.txt"
;

Variables
z, v(j), slL(i), slU(i)
;

****************** SETTING LOWER AND UPPER BOUNDS FOR VARIABLES ******************
v.lo(j) = LB(j);
v.up(j) = UB(j);

v.fx('R__BIOMASS_Ec_iML1515_WT_75p37M') = 0;
v.lo('R__BIOMASS_Ec_iML1515_core_75p37M') = 0;
v.up('R__BIOMASS_Ec_iML1515_core_75p37M') = 0;
v.lo('R__EX_glc__D_e') = -10; v.up('R__EX_glc__D_e') = 1000;
v.lo('R__EX_o2_e') = -1000; v.up('R__EX_o2_e') = 1000;
v.lo('R__EX_nh4_e') = -1000; v.up('R__EX_nh4_e') = 1000;
v.fx('R__ATPM') = 6.86;

slL.lo(i) = 0; slL.up(i) = 10;
slU.lo(i) = 0; slU.up(i) = 10;

* Disable imbalance of extracellular metab
slL.fx(i)$(metext(i)) = 0;
slU.fx(i)$(metext(i)) = 0;

* Set imbalance
slL.fx('M__zn2_p') = 0.01;
slU.fx('M__zn2_p') = 0;

***************************** EQUATION DEFINITIONS **************************************
Equations
Obj, Stoic
;

Obj..			z =e= sum(i, slL(i) + slU(i));
Stoic(i)..		sum(j, S(i,j)*v(j)) =e= slU(i) - slL(i);


************************ DECLARING THE MODEL with constraints****************************
Model adjustrhs
/Obj, Stoic/;
adjustrhs.optfile = 1;

******************* SOLVE *******************
Solve adjustrhs using lp minimizing z;

******************* WRITE TO TXT FILE *******************
file ff /coupling_imbal.txt/;
put ff;
loop(i,
	if ( (slL.l(i) gt 1E-8),
		put i.tl:0, system.tab, 'slL', system.tab, slL.l(i):0:8/;
);
	if ( (slU.l(i) gt 1E-8),
		put i.tl:0, system.tab, 'slU', system.tab, slU.l(i):0:8/;
);
);

file ff2 /check_infes_imbal.modelStat.txt/;
put ff2;
put adjustrhs.modelStat/;
putclose ff2;
