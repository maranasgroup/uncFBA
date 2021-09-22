# uncFBA example scripts repository
This GitHub repository provides example scripts for the manuscript **"Quantifying the propagation of parametric uncertainty on flux balance analysis"**.<br>
<br>
In the repository are example scripts that perform uncertainty injection and propagation to flux balance analysis with outputs for a small sample size (for demonstration purpose only). For proper analysis, user should download the scripts and run for a large sample size (e.g., 10,000 samples).<br>
<br>
If you use the scripts, please cite the following bioRxiv preprint: "Quantifying the propagation of parametric uncertainty on flux balance analysis" (https://www.biorxiv.org/content/10.1101/2021.09.16.460685v1)<br>
<br>
There are two subdirectories:
1) ***/uncFBA/uncBiom***: injection of normally distributed noise to biomass precursor coeffcients and ATP maintenance (growth-associated ATP maintenance (GAM) and non-growth associated ATP maintenance (NGAM))
2) ***/uncFBA/uncRHS***: departure from steady-state by adding noise drawn from normal distribution to the RHS terms of mass balance constraints
<br>
More detailed instructions are provided in the "README.md" inside each subdirectory and in the scripts.
