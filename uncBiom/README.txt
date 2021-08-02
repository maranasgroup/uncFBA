Description of files in uncFBA/uncBiom
The files in here are example scripts to sample, perform pFBA, and analyze uncertainty propagation from biomass composition and ATP demand uncertainty to fluxes. The output files are provided to just show the format (with only 50 samples). For proper analysis, users should download the scripts and run for 10,000 samples or more.

There are the following cases:
1) CoeffVary: Apply uncertainty to all individual biomass constituents' fractions simultaneously
2) CoeffVaryNoRescale: Apply uncertainty to all individual biomass constituents' fractions simultaneously, without rescaling the biomass molecular weight (MW) to 1 g/mmol
3) SingleCoeffVary: Apply uncertainty to only a designated single individual biomass constituent
4) MacroVary: Apply uncertainty to all macromolecular fractions for all macromolecular classes simultaneously
5) SingleMacroVary: Apply uncertainty to only a designated macromolecular fraction
6) GAMVary: Apply uncertainty to the growth-associated ATP maintenance rate (GAM)
7) NGAMVary: Apply uncertainty to the non-growth associated ATP maintenance rate (NGAM)

There are the following scripts:
- A01_sample_and_pFBA.ipynb: Sample parameters, modify model, and run pFBA
- A02_calc_fluxSum.ipynb: Calculate flux sum quantities for alternative pathways
- A03_calc_SDR.ipynb: Calculate SDR, which is used to quantify the propagated uncertainty from parameters to biomass yield and fluxes

There are the following outputs:
Sampled parameters
Columns: biomass MW and metabolites in biomass coefficient
Rows: indexes of samples. Row 0: is the original coefficients from the iML1515 model, not a part of the sampled population.
- ftype-dfCoeff_CoeffStd-0.1.csv: csv file storing all sampled parameters for cases of biomass composition types
- ftype-dfCoeff_GAMStd-0.1.csv: csv file storing all sampled biomass coefficients after GAM being sampled for GAMVary. This is necessary because GAM go into biomass reaction with biomass composition.
- For NGAM uncertainty, it is only a single value applied to the lower bound of ATPM reaction. Sampling is implemented directly in the script without outputing an intermediate file.

Calculated pFBA flux distributions:
Columns: Reaction ID
Rows: consistent with from sampled parameters. See the description there.
- ftype-dfCoeff_CoeffStd-0.1.csv or ftype-dfFlux_GAMStd-0.1.csv

Calculated flux sum quantities:
Columns / Rows: same with above for pFBA results
- ftype-dfCoeff_CoeffStd-0.1_added.csv or ftype-dfFlux_GAMStd-0.1_added.csv

Quantifying propagated uncertainty
Columns:
- Stdev_slope: slope of linear regression of relative standard deviation (SD) vs. uncertainty levels (0%, 5%, 10%, 20%, and 30%). This is the SDR mentioned in the manuscript.
- R2: R2-score of linear regression
- Stdev_0.05/0.1/0.2/0.3: absolute SD of a flux at a specific uncertainty level
- Mean_0.05/0.1/0.2/0.3: mean of a flux at a specific uncertainty level
- Stdevnorm_0.05/0.1/0.2/0.3: relative SD  = absolute SD / mean

Notes:
- A01 scripts is provided for all cases. A02 and A03 scripts are only provided for CoeffVary. Users can download and adjust the A02 and A03 scripts for a particular case that they are interested in.
- For proper results, a large enough population size from sampling is required. We recommend 10,000.
- In CoeffVary example, we provided dfFluxStat.csv analyzed from running the population size of 10,000 across uncertainty levels of 5%, 10%, 20%, and 30%.
- In all cases, we provided a sample output for A01 scripts at 10% for the population size of 50 so that the users can observe the format. Again, we recommend running for population size of 10,000.
