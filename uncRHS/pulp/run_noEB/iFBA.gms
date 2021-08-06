************************* Adjusting RHS using Linear Programming ********************
*       Author: Hoang Dinh
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
j
$include "%path%reactions.txt"
;

Parameters
S(i,j)
$include "%path%sij.txt"
b(i)
$include "imbal_adj.txt"
LB(j)
$include "%path%lb.txt"
UB(j)
$include "%path%ub.txt"
;

Variables
z, v(j)
;

****************** SETTING LOWER AND UPPER BOUNDS FOR VARIABLES ******************
v.lo(j) = LB(j)*%nscale%;
v.up(j) = UB(j)*%nscale%;

v.fx('R__BIOMASS_Ec_iML1515_WT_75p37M') = 0;
v.lo('R__BIOMASS_Ec_iML1515_core_75p37M') = 0;
v.up('R__BIOMASS_Ec_iML1515_core_75p37M') = 1000*%nscale%;
v.lo('R__EX_glc__D_e') = -10*%nscale%; v.up('R__EX_glc__D_e') = 1000*%nscale%;
v.lo('R__EX_o2_e') = -1000*%nscale%; v.up('R__EX_o2_e') = 1000*%nscale%;
v.lo('R__EX_nh4_e') = -1000*%nscale%; v.up('R__EX_nh4_e') = 1000*%nscale%;
v.fx('R__ATPM') = 6.86*%nscale%;

***************************** EQUATION DEFINITIONS **************************************
Equations
Obj, Stoic
;

Obj..			z =e= v('R__BIOMASS_Ec_iML1515_core_75p37M');
Stoic(i)..		sum(j, S(i,j)*v(j)) =e= b(i)*%nscale%;


************************ DECLARING THE MODEL with constraints****************************
Model ifba
/Obj, Stoic/;
ifba.optfile = 1;

******************* SOLVE *******************
Solve ifba using lp maximizing z;

******************* WRITE TO TXT FILE *******************
file ff /ifba.txt/;
put ff;

loop(j,
	put j.tl:0, system.tab, v.l(j):0:8/;
);
putclose ff;


file ff2 /iFBA.modelStat.txt/;
put ff2;
put ifba.modelStat/;
putclose ff2;

