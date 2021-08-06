import numpy as np
import pandas as pd
import os, shutil

# Set path and run_id
from path_and_runid import *

# Adjust RHS
df_rhs = pd.read_csv(path_A01, sep='\t', index_col=0)
df_rhs_adj = pd.DataFrame(index=df_rhs.index.to_list(), columns=['modelStat'] + df_rhs.columns.to_list())
df_rhs_adj.loc[0,:] = df_rhs.loc[0,:]

shutil.copy('../GAMS/pulp.gms', './pulp.gms')
shutil.copy('../GAMS/cplex.opt', './cplex.opt')

for i in range(1, df_rhs_adj.shape[0]):
    
    # Write initial imbalance for GAMS
    with open(path_GAMS_imbal, 'w') as f:
        text = ['/']
        for met in df_rhs.columns:
            text.append("'M__" + met + "' " + str(round(df_rhs.loc[i,met],7)))
        text += ['/']
        f.write('\n'.join(text))
        
    # Run LP adjustment in GAMS
    os.system('module load gams\n' + 'gams ' + path_GAMS_adjustRHS + ' o=/dev/null')
    #os.system('module load gams\n' + 'gams ' + path_GAMS_adjustRHS)
    
    # Read GAMS output and use it to calculate LP-adjusted RHS
    with open(path_GAMS_sl_vals) as f:
        x = f.read().split('\n')
    x = [line for line in x if line != '']

    mets_in = [line.split('\t')[0] for line in x]
    mets_in = sorted(list(set(mets_in)))
    
    sl_dict = {met:0 for met in mets_in}
    for line in x:
        met,vtype,v = line.split('\t')
        if vtype == 'slL':
            sl_dict[met] += -float(v)
        elif vtype == 'slU':
            sl_dict[met] += float(v)
        else:
            print(i, 'check, something is wrong')
    sl_dict = {met[3:]:v for met,v in sl_dict.items()}

    mets = [met for met in df_rhs.columns if met not in sl_dict.keys()]
    for met in mets:
        sl_dict[met] = 0
    sl_pdseries = pd.core.series.Series(sl_dict)
    
    df_rhs_adj.loc[i,:] = df_rhs.loc[i,:] + sl_pdseries.astype(float) / 1e3
    
    # Record modelStat
    with open(path_GAMS_adjustRHS_modelStat) as f:
        x = f.read().replace(' ', '')
        x = x.replace('\n', '')
        df_rhs_adj.loc[i, 'modelStat'] = int(float(x))

df_rhs_adj.to_csv(path_A02, sep='\t')
