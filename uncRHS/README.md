# Instructions to run uncFBA/uncRHS
The files in here are example scripts to generate RHS vectors from normal distributions accounting for metabolite pools and elemental balances (via projection using linear programming (PULP)), and find flux distributions and biomass yield correspond to the generated RHS via instationary FBA (iFBA). The output files are provided to just show the format (with only 50 samples). For proper analysis, users should download the scripts and run for 10,000 samples or more.

There are the following cases:
1) "run": run PULP and iFBA, enforcing elemental balance constraints
2) "run_noEB": run PULP and iFBA, ignoring elemental balance constraints
3) "maxBiom": run a modified formulation of PULP that set the optimization objective to be maximization of biomass yield, enforcing elemental balance constraints
4) "maxBiom_noEB": run a modified formulation of PULP that set the optimization objective to be maximization of biomass yield, ignoring elemental balance constraints

To run "run" and "run_noEB", here are the instructions:
1) Adjust the settings in "path_and_run_id.py:
2) "python A01_initiate_rhs.py": get initial randomly sampled RHS vectors from normal distribution (generate "A01_rhs_init.csv")
3) "python A02_LP_adjust.py": adjust the randomly sampled RHS vectors with PULP (generate "A02_RHS_adjust.csv")
4) "python A03_iFBA.py": run iFBA with input of adjusted RHS vectors (generate "A03_iFBA_fluxes.csv")

To run "maxBiom" and "maxBiom_noEB", here are the instructions:
1) Adjust the upper bounds of the slack variables so as to constrain the absolute values of RHS terms to be within the bounds, in the "maxSL.txt" file. Provided in the example scripts are upper bounds in mmol/gDW/h that are derived from the mass-basis amount of 0.1% of glucose uptake rate by mass.
2) "gams pulp_maxBiom.gms": get the max biomass yield under metaboltie unsteady-state
3) Copy the max biomass yield from the outfile "pulp_maxBiom.maxGR.txt" into the file "pulp_minSL_to_maxBiom.gms" at line
```
v.fx('R__BIOMASS_Ec_iML1515_core_75p37M') = 924.01715914 - 1e-4;
```
(-1e-4 is done to avoid number rounding issue of optimal solution that could affect optimization feasibility)<br>
4) "gams pulp_minSL_to_maxBiom.gms": find minimal adjustments from the steady-state subject to the max biomass yield found previously
