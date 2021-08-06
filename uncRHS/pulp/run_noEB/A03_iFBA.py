import numpy as np
import pandas as pd
import os, cobra, shutil

# Set path and run_id
from path_and_runid import *

# Load adjusted RHS
df_rhs_adj = pd.read_csv(path_A02, sep='\t', index_col=0)

# Rxn list
model = cobra.io.load_json_model(path_model)
rxns = sorted([rxn.id for rxn in model.reactions])

# iFBA simulation and record flux
df_flux = pd.DataFrame(index=df_rhs_adj.index.to_list(), columns=['modelStat'] + rxns)

shutil.copy('../GAMS/iFBA.gms', './iFBA.gms')
shutil.copy('../GAMS/cplex.opt', './cplex.opt')

for i in range(1, df_rhs_adj.shape[0]):
    
    # Write adjusted imbalance for GAMS
    with open(path_GAMS_imbal_adj, 'w') as f:
        text = ['/']
        mets = [col for col in df_rhs_adj.columns if col != 'modelStat']
        for met in mets:
            text.append("'M__" + met + "' " + str(round(df_rhs_adj.loc[i,met],10)))
        text += ['/']
        f.write('\n'.join(text))
        
    # Run iFBA in GAMS
    os.system('module load gams\n' + 'gams ' + path_GAMS_iFBA + ' o=/dev/null')
    #os.system('module load gams\n' + 'gams ' + path_GAMS_iFBA)
    
    # Read GAMS output and record flux to DataFrame
    with open(path_GAMS_ifbaflux) as f:
        x = f.read().split('\n')
    x = [line for line in x if line != '']
    
    df_flux.loc[i, :] = pd.Series({line.split('\t')[0][3:]:line.split('\t')[1] for line in x}).astype(float) / 1e3
    
    # Record modelStat
    with open(path_GAMS_iFBA_modelStat) as f:
        x = f.read().replace(' ', '')
    x = x.replace('\n', '')
    df_flux.loc[i, 'modelStat'] = int(float(x))

df_flux.to_csv(path_A03, sep='\t')
