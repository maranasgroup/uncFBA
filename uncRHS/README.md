# Instructions to run uncFBA/uncRHS
The files in here are example scripts to generate RHS vectors from normal distributions accounting for metabolite pools and elemental balances (via projection using linear programming (PULP)), and find flux distributions and biomass yield correspond to the generated RHS via instationary FBA (iFBA). The output files are provided to just show the format (with only 50 samples). For proper analysis, users should download the scripts and run for 10,000 samples or more.

**There are the following optimization formulations**:
1) ***./cmpcheck***: CMP-check optimization to find if a metabolite is in a pool by forcing a metabolite to be under accumulation or depletion and observe if another metabolite is simultaneously in imbalance
2) ***./cmpfind***: CMP-find optimization to exhaustively find all minimal sets of metabolites whose RHS terms are coupled
3) ***./pulp***: Projection using linear programming (PULP) to adjust a randomly sampled RHS vector accounting for metabolite pools and elemental balance

**Instructions to run**:<br>
**1) Run CMP-check**<br>
Go inside ./cmpcheck directory and run "python cmpcheck.py". The output is in "cmpcheck_results.csv"<br><br>
**2) Run CMP-find**<br>
- Go inside ./cmpfind directory
- Create a subdirectory in the following format "<metabolite_id>.<*bneg* or *bpos*>" where *bneg* indicates negative RHS perturbation (depletion) and *bpos* indicates positive RHS perturbation (accumulation)
- Copy the file *cmpfind.py* from the example scripts for "2dhp_c" depletion in the subdirectory "./cmpfind/2dhp_c.bneg"
- Edit the newly copied *cmpfind.py* file at lines to your selected metabolite and imbalance (look at CMP-check results to see which metabolite and imbalance combinations you need to run for) (the below example demonstrated for tetrahydrofolate depletion)
```
# Set metabolite and imbalance
met0 = 'thf_c'
imb0 = 'bneg'
```
- Run "python cmpfind.py" and the output is in "cmpfind_solns.txt"

**3) Run PULP**<br>
Go inside the ./pulp subdirectory for detailed instruction in a separate "README.md" file<br>
