import pandas as pd
import os, shutil
from copy import deepcopy

def build_integer_cut_from_slvars(slv_text):
    slvs = slv_text.split(', ')

    cut_text = []
    for slv in slvs:
        x,_ = slv.split(':')
        met,imb = x.split('.')
        cut_text.append(yvar_conv[imb] + "('M__" + met + "')")

    size = len(cut_text)
    cut_text = ' + '.join(cut_text)
    cut_text += ' =l= ' + str(size-0.5) + ';'
    return cut_text

# Set metabolite and imbalance
met0 = '2dhp_c'
imb0 = 'bneg'

# Conv dictionaries
imb_conv = {'slL': 'bneg', 'slU': 'bpos'}
imb_conv_rev = {v:k for k,v in imb_conv.items()}
yvar_conv = {'bneg': 'yL', 'bpos': 'yU'}

# Initiate
shutil.copy('../GAMS/cmpfind.gms', './cmpfind.gms');
shutil.copy('../GAMS/cplex.opt', './cplex.opt');

# Initiate file
with open('./cmpfind.gms') as f:
    gms = f.read().split('\n')
    
if imb0 == 'bpos':
    gms[68] = "slL.fx('M__" + met0 + "') = 0;"
    gms[69] = "slU.fx('M__" + met0 + "') = 0.01;"
    gms[70] = "yL.fx('M__" + met0 + "') = 0;"
    gms[71] = "yU.fx('M__" + met0 + "') = 1;"
elif imb0 == 'bneg':
    gms[68] = "slL.fx('M__" + met0 + "') = 0.01;"
    gms[69] = "slU.fx('M__" + met0 + "') = 0;"
    gms[70] = "yL.fx('M__" + met0 + "') = 1;"
    gms[71] = "yU.fx('M__" + met0 + "') = 0;"
    
with open('./cmpfind.gms', 'w') as f:
    f.write('\n'.join(gms))
    
# Run GAMS first time no integer cut
os.system('module load gams\n' + 'gams cmpfind.gms o=/dev/null')
#os.system('module load gams\n' + 'gams cmpfind.gms')

with open('./cmpfind.modelStat.txt') as f:
    modelStat = f.read()
modelStat = modelStat.replace('\n', '')
modelStat = modelStat.replace(' ', '')
modelStat = int(float(modelStat))

if modelStat in [4,10,11,12,13,14,19]:
    with open('./cmpfind_msgOut.txt', 'w') as f:
        f.write('Terminate due to modelStat of ' + str(modelStat))
    quit()
    
cols = ['id', 'modelStat', 'sl_vars', 'cut_text']
df_solns = pd.DataFrame(columns=cols)

count = 0
df_solns.loc[count, ['id', 'modelStat']] = [count, modelStat]

with open('./cmpfind_coupling_imbal.txt') as f:
    coup_imbal = f.read().split('\n')
coup_imbal = [i for i in coup_imbal if i != '']

text = []
for entry in coup_imbal:
    met,imb,v = entry.split('\t')
    v = float(v)
    text.append(met[3:] + '.' + imb_conv[imb] + ':' + str(v))
text = ', '.join(text)

df_solns.loc[count, 'sl_vars'] = text
df_solns.loc[count, 'cut_text'] = build_integer_cut_from_slvars(text)

# Output
with open('./cmpfind_solns.txt', 'w') as f:
    text = [str(i) for i in df_solns.loc[count,:]]
    f.write('\t'.join(text) + '\n')
    
# Run cmpfind, per iteration implemented more and more integer cuts
# Enable integer cuts in GAMS
gms[77] = '$include "cmpfind_integerCuts_declare.txt"'
gms[86] = '$include "cmpfind_integerCuts_eqns.txt"'
gms[92] = '$include "cmpfind_integerCuts_declare.txt"'
with open('./cmpfind.gms', 'w') as f:
    f.write('\n'.join(gms))
    
running = True
slvs_str = df_solns.sl_vars.to_list()[-1]
mets_prev = set([slv.split('.')[0] for slv in slvs_str.split(', ')])

#while running and count < 2:
while running:
    # Write integer cuts to text files
    with open('./cmpfind_integerCuts_declare.txt', 'w') as f:
        f.write('\n'.join(['cut' + str(i) for i in df_solns.index]))

    text = []
    for i in df_solns.index:
        text.append('cut'+str(i) + '.. ' + df_solns.cut_text[i])
    with open('./cmpfind_integerCuts_eqns.txt', 'w') as f:
        f.write('\n'.join(text))

    # Run GAMS
    os.system('module load gams\n' + 'gams cmpfind.gms o=/dev/null')
    #os.system('module load gams\n' + 'gams cmpfind.gms')
    
    # Record
    count += 1
    df_solns.loc[count, 'id'] = count
    
    with open('./cmpfind.modelStat.txt') as f:
        modelStat = f.read()
    modelStat = modelStat.replace('\n', '')
    modelStat = modelStat.replace(' ', '')
    modelStat = int(float(modelStat))
    df_solns.loc[count, 'modelStat'] = modelStat
    
    with open('./cmpfind_coupling_imbal.txt') as f:
        coup_imbal = f.read().split('\n')
    coup_imbal = [i for i in coup_imbal if i != '']

    text = []
    for entry in coup_imbal:
        met,imb,v = entry.split('\t')
        v = float(v)
        text.append(met[3:] + '.' + imb_conv[imb] + ':' + str(v))
    text = ', '.join(text)

    df_solns.loc[count, 'sl_vars'] = text
    df_solns.loc[count, 'cut_text'] = build_integer_cut_from_slvars(text)
    
    # Output
    with open('./cmpfind_solns.txt', 'a') as f:
        text = [str(i) for i in df_solns.loc[count,:]]
        f.write('\t'.join(text) + '\n')
    
    # Check for termination
    if modelStat in [4,10,11,12,13,14,19]:
        with open('./cmpfind_msgOut.txt', 'w') as f:
            f.write('Terminate due to modelStat of ' + str(modelStat))
        running = False
    
    # Check if the run stuck, behave with duplicate solutions being found
    # continuously, probably due to solver tolerance
    slvs_str = df_solns.sl_vars.to_list()[-1]
    mets_now = set([slv.split('.')[0] for slv in slvs_str.split(', ')])
    if mets_now == mets_prev:
        running = False
    else:
        mets_prev = deepcopy(mets_now)
