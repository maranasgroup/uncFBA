import numpy as np
import pandas as pd
import os, shutil

# Metabolite list
#path_model = '../input/iML1515_raw.json'
#model = cobra.io.load_json_model(path_model)
#mets_check = sorted([met.id for met in model.metabolites if met.compartment != 'e'])
# List of all intracellular metabolites can be extracted directly from the json model using the above scripts
# Due to an error from uninstalled packages, I cannot load cobrapy package
# I will load the metabolite list from a text file instead
with open('./mets_check.txt') as f:
    mets_check = f.read().split('\n')
imbals_check = sum([[met+'.bpos', met+'.bneg'] for met in mets_check], [])

# Simulation w/ imbalance slackness and record status
cols = ['id', 'modelStat', 'sl_vars']
df_stat = pd.DataFrame(index=imbals_check, columns=cols)
df_stat['id'] = df_stat.index.to_list()

# Initiate
shutil.copy('./GAMS/cmpcheck.gms', './cmpcheck.gms');
shutil.copy('./GAMS/cplex.opt', './cplex.opt');

# Read initial GAMS file
with open('./cmpcheck.gms') as f:
    gams_file = f.read().split('\n')
gams_file = gams_file[:-1]

for i in df_stat.index:
    met,imbal = df_stat.id[i].split('.')
    if imbal == 'bpos':
        gams_file[65] = "slL.fx('M__" + met + "') = 0;"
        gams_file[66] = "slU.fx('M__" + met + "') = 0.01;"
    elif imbal =='bneg':
        gams_file[65] = "slL.fx('M__" + met + "') = 0.01;"
        gams_file[66] = "slU.fx('M__" + met + "') = 0;"
    else:
        print(i, ', check your input')
        continue

    # Write to a GAMS file
    with open('./cmpcheck.gms', 'w') as f:
        f.write('\n'.join(gams_file))
    
    # Run check in GAMS
    os.system('module load gams\n' + 'gams cmpcheck.gms o=/dev/null')
    
    # Record modelStat
    with open('check_infes_imbal.modelStat.txt') as f:
        x = f.read().replace(' ', '')
    x = x.replace('\n', '')
    df_stat.loc[i, 'modelStat'] = int(float(x))
    
    # Record slack_vars
    with open('coupling_imbal.txt') as f:
        x = f.read().split('\n')
    x = [line for line in x if line != '']
    x = [line.split('\t') for line in x]
    
    slack_vars = []
    for entry in x:
        if entry[0][3:] == met:
            continue
        met_coup = entry[0][3:]
        if entry[1] == 'slL':
            imbal_type = 'bneg'
        elif entry[1] == 'slU':
            imbal_type = 'bpos'
        else:
            print(i, ', something is wrong, check')
        slack_vars.append(met_coup + '.' + imbal_type + ':' + entry[2])

    df_stat.loc[i, 'sl_vars'] = ','.join(slack_vars)

df_stat.to_csv('./cmpcheck_results.csv', sep='\t', index=False)
