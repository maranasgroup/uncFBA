import pandas as pd
import numpy as np
import os
import cobra

# Set path and run_id
from path_and_runid import *

### Load items
# Model
model = cobra.io.load_json_model(path_model)
model.solver = 'cplex'
mets_model = sorted([met.id for met in model.metabolites])

biomId = 'BIOMASS_Ec_iML1515_core_75p37M'

model.reactions.EX_glc__D_e.bounds = (-10, 1000)
model.reactions.EX_o2_e.bounds = (-1000, 1000)
model.reactions.EX_nh4_e.bounds = (-1000, 1000)

model.objective = dict()
model.reactions.get_by_id(biomId).bounds = (0, 1000)
model.reactions.get_by_id(biomId).objective_coefficient = 1
fba = model.optimize()
mu_max = fba.fluxes[biomId]

# Molecular weight
df_mw = pd.read_csv(path_metMW, sep='\t')
df_mw.index = df_mw.id.to_list()

# Imbalance status
# Status of imbalance
df_stat = pd.read_csv(path_metStat, sep='\t')
df_stat.index = df_stat.id.to_list()
mets_stat = sorted(list(set([i.split('.')[0] for i in df_stat.index])))
    # Equivalent to set of metabolites that are not external

### Sampling
samp_methods = dict()
for met in mets_stat:
    bneg_infes = df_stat.modelStat[met+'.bneg'] == 4
    bpos_infes = df_stat.modelStat[met+'.bpos'] == 4
    if bneg_infes and bpos_infes:
        samp_methods[met] = 'zero'
    elif bneg_infes:
        samp_methods[met] = 'half_normal_positive'
    elif bpos_infes:
        samp_methods[met] = 'half_normal_negative'
    else:
        samp_methods[met] = 'normal'

from scipy.stats import norm,halfnorm

dfSamp = pd.DataFrame(index=range(0,size+1), columns=mets_stat)
dfSamp.loc[0,:] = [0]*dfSamp.shape[1]

for met in dfSamp.columns:
    MW = df_mw.MW[met]
    if MW == 0:
        dfSamp.loc[1:, met] = [0]*size
        continue
    else:
        scale_adj = scale*1.80156 * (1000./MW)
    
    if samp_methods[met] == 'zero':
        dfSamp.loc[1:, met] = [0]*size
    elif samp_methods[met] == 'normal':
        dfSamp.loc[1:, met] = scale_adj*norm.rvs(size=size)
    elif samp_methods[met] == 'half_normal_positive':
        dfSamp.loc[1:, met] = scale_adj*halfnorm.rvs(size=size)
    elif samp_methods[met] == 'half_normal_negative':
        dfSamp.loc[1:, met] = -scale_adj*halfnorm.rvs(size=size)
    else:
        print(met, ', need to check something is wrong')

dfSamp = dfSamp.astype(float).round(7)
dfSamp.to_csv(path_A01, sep='\t', index=True, header=True)
