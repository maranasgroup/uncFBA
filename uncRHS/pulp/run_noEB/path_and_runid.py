import os

# Set path, run_id, and sampling parameters
path = '../'
size = 50
scale = 0.001

path_model = os.path.join(path, 'input/iML1515_raw.json')
path_metMW = os.path.join(path, 'input/met_mw.csv')
path_metStat = os.path.join(path, 'input/cmpcheck_results.csv')
path_A01 = os.path.join('./A01_rhs_init.csv')
path_A02 = os.path.join('./A02_RHS_adjust.csv')
path_A03 = os.path.join('./A03_iFBA_fluxes.csv')


path_GAMS_imbal = os.path.join('./imbal.txt')
path_GAMS_adjustRHS = os.path.join('./pulp_noEB.gms')
path_GAMS_adjustRHS_modelStat = os.path.join('./pulp.modelStat.txt')
path_GAMS_sl_vals = os.path.join('./sl_vals.txt')
path_GAMS_imbal_adj = os.path.join('./imbal_adj.txt')
path_GAMS_iFBA = os.path.join('./iFBA.gms')
path_GAMS_iFBA_modelStat = os.path.join('./iFBA.modelStat.txt')
path_GAMS_ifbaflux = os.path.join('./ifba.txt')
