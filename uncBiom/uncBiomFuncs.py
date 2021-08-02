import pandas as pd

def file_name(coeff_std, gam_std, ftype=None, size=None, remark=None):
    fname_dict = {'ftype':ftype, 'CoeffStd':coeff_std, 'GAMStd':gam_std,
                  'Size':size, 'Remark':remark}
    fname_elem = []
    for elem in ['ftype', 'CoeffStd', 'GAMStd', 'Size', 'Remark']:
        if fname_dict[elem] not in [None, False, '']:
            fname_elem.append(elem + '-' + str(fname_dict[elem]))
    fname = '_'.join(fname_elem)
    return fname

def get_coeff_without_gam(model, biomId, gam_val):
    from collections import OrderedDict 
    import numpy as np

    atpm = ['atp_c', 'adp_c', 'h2o_c', 'pi_c', 'h_c']
    biomrxn = model.reactions.get_by_id(biomId)
    mets = atpm + sorted([met.id for met in biomrxn.metabolites.keys() if met.id not in atpm])

    bMets = OrderedDict()
    for metid in mets:
        met = model.metabolites.get_by_id(metid)
        bMets[met.id] = biomrxn.metabolites[met]

    for metid in atpm:# Exclude ATP maintainance
        bMets[metid] = bMets[metid] - np.sign(bMets[metid])*gam_val
        
    bMets_out = OrderedDict()
    for k,v in bMets.items():
        if v != 0:
            bMets_out[k] = v
            
    return bMets_out
    
def build_reaction_equation_from_metabolites_dict(met_dict, arrow='<=>', floatdecimal=6):
    lhs = []; rhs = [];
    for k,v in met_dict.items():
        v = float(v)
        if v == -1:
            lhs.append(k)
        elif v == 1:
            rhs.append(k)
        elif v < 0 and v != -1 and v.is_integer():
            lhs.append(' '.join([str(-int(v)), k]))
        elif v > 0 and v != 1 and v.is_integer():
            rhs.append(' '.join([str(int(v)), k]))
        elif v < 0 and v != -1:
            lhs.append(' '.join([('{:.' + str(floatdecimal) + 'f}').format(-v), k]))
        elif v > 0 and v != 1:
            rhs.append(' '.join([('{:.' + str(floatdecimal) + 'f}').format(v), k]))
    return ' '.join([ ' + '.join(lhs), arrow, ' + '.join(rhs)])

def make_normal_distributed_coefficients(model, coeff_std, gam_val,
    biomId, mets_select='all', size=10000, bal_adjust=True, mw_adjust=True, zero_tol=1e-9):
    # Index 0 (first line): original coeffs / solution obtained using original coeffs
    # Index 1 ...: index of randomized coeffs / solutions obtained using randomized coeffs
    # This will not renormalize molecular weight
    
    # Get list of metabolites in biomass equation (except growth associated maintenance)
    import numpy as np
    from numpy.random import normal
    
    round_digit = -int(np.log10(zero_tol))
    bMets = get_coeff_without_gam(model, biomId, gam_val)

    # Normally distributed biomass equation coefficients
    mlist = []
    if mets_select == 'all':
        mlist = bMets.keys()
    else:
        mlist = mets_select
        if len(set(mlist) - set(bMets.keys())) > 0:
            raise ValueError('Contain metabolite that is not in the list')

    dfCoeff = pd.DataFrame(index = range(0,size+1), 
        columns=['biomMW'] + list(bMets.keys()))

    i = 0
    biomMW = 0
    for col in dfCoeff.columns[1:]:
        dfCoeff.loc[i, col] = bMets[col]

    seed = -1
    for col in dfCoeff.columns[1:]:
        if col not in mlist:
            dfCoeff.loc[1:, col] = dfCoeff.loc[0, col]
        else:
            if col in ['h2o_c', 'ppi_c']:
                continue
            c0 = dfCoeff.loc[0, col]
            seed += 1
            np.random.seed(seed)
            coeff_norm = c0 + coeff_std*c0*normal(size=size)
            # Only allow negative value for coefficients of biomass constituents
            coeff_norm = [coeff if coeff < -zero_tol else -zero_tol for coeff in coeff_norm]
            dfCoeff.loc[range(1,size+1), col] = coeff_norm

    # Coefficients dependencies on each others in biomass equation (active when bal_adjust=True)
    if bal_adjust:
        aas = ['ala__L','arg__L','asn__L','asp__L','cys__L','gln__L','glu__L','gly','his__L','ile__L',
               'leu__L','lys__L','met__L','phe__L','pro__L','ser__L','thr__L','trp__L','tyr__L','val__L']
        aas = {k + '_c':1 for k in aas}

        nus = ['datp','dctp','dgtp','dttp','ctp','gtp','utp','atp']
        nus = {k + '_c':1 for k in nus}

        for i in range(1, size+1):
            coeff = 0
            for k,v in aas.items():
                coeff += v*dfCoeff.loc[i, k]
            dfCoeff.loc[i, 'h2o_c'] = -coeff

            coeff = 0
            dfCoeff.loc[i, 'ppi_c'] = 0
            for k,v in nus.items():
                coeff += v*dfCoeff.loc[i, k]
            dfCoeff.loc[i, 'ppi_c'] = -coeff

    elif bal_adjust == False:
        seed = 999
        for col in ['h2o_c', 'ppi_c']:
            c0 = dfCoeff.loc[0, col]
            seed += 1
            np.random.seed(seed)
            dfCoeff.loc[range(1,size+1), col] = c0 + coeff_std*c0*normal(size=size)

    # Calculate biomass's MW
    for i in range(1, size+1):
        biomMW = 0
        for col in dfCoeff.columns[1:]:
            met = model.metabolites.get_by_id(col)
            biomMW += dfCoeff.loc[i, col] * met.formula_weight
        dfCoeff.loc[i, 'biomMW'] = -biomMW

    # Normalize to 1 g/mmol
    if mw_adjust:
        cols = [col for col in dfCoeff.columns if col != 'biomMW']
        for i in range(1, size+1):
            biomMW = dfCoeff.loc[i, 'biomMW']
            for col in cols:
                coeff = round(dfCoeff.loc[i, col] * 1000 / biomMW, round_digit)
                dfCoeff.loc[i, col] = coeff

    if bal_adjust:
        aas = ['ala__L','arg__L','asn__L','asp__L','cys__L','gln__L','glu__L','gly','his__L','ile__L',
               'leu__L','lys__L','met__L','phe__L','pro__L','ser__L','thr__L','trp__L','tyr__L','val__L']
        aas = {k + '_c':1 for k in aas}

        nus = ['datp','dctp','dgtp','dttp','ctp','gtp','utp','atp']
        nus = {k + '_c':1 for k in nus}

        for i in range(1, size+1):
            coeff = 0
            for k,v in aas.items():
                coeff += v*dfCoeff.loc[i, k]
            dfCoeff.loc[i, 'h2o_c'] = -coeff

            coeff = 0
            dfCoeff.loc[i, 'ppi_c'] = 0
            for k,v in nus.items():
                coeff += v*dfCoeff.loc[i, k]
            dfCoeff.loc[i, 'ppi_c'] = -coeff

    # Add GAM to the coefficient
    for col in ['atp_c', 'adp_c', 'h2o_c', 'pi_c', 'h_c']:
        if col not in dfCoeff.columns:
            dfCoeff[col] = [0.]*dfCoeff.shape[0] # Create array of zeros
    for metId in ['atp_c', 'h2o_c']:
        dfCoeff.loc[:, metId] -= gam_val
    for metId in ['adp_c', 'pi_c', 'h_c']:
        dfCoeff.loc[:, metId] += gam_val

    atpm = ['atp_c', 'adp_c', 'h2o_c', 'pi_c', 'h_c']
    biomrxn = model.reactions.get_by_id(biomId)
    mets = atpm + sorted([met.id for met in biomrxn.metabolites.keys() if met.id not in atpm])
    dfCoeff = dfCoeff.reindex(['biomMW'] + mets, axis=1)
    
    return dfCoeff

def make_normal_distributed_gam(model, gam_val, gam_std, biomId, size=10000, zero_tol=1e-9):
    # Index 0 (first line): original coeffs / solution obtained using original coeffs
    # Index 1 ...: index of randomized coeffs / solutions obtained using randomized coeffs
    # This will not renormalize molecular weight
    
    # Get list of metabolites in biomass equation (except growth associated maintenance)
    import numpy as np
    from numpy.random import normal
    from scipy.stats import tmean, tstd

    bMets = get_coeff_without_gam(model, biomId, gam_val)
    dfCoeff = pd.DataFrame(index = range(0,size+1), 
        columns=['biomMW'] + list(bMets.keys()))
    vals = [bMets[col] for col in dfCoeff.columns[1:]]
    cols_slice = dfCoeff.columns[1:]
    dfCoeff.loc[:, cols_slice] = vals

    for col in ['atp_c', 'adp_c', 'h2o_c', 'pi_c', 'h_c']:
        if col not in dfCoeff.columns:
            dfCoeff[col] = [0.]*dfCoeff.shape[0] # Create array of zeros

    for metId in ['atp_c', 'h2o_c']:
        dfCoeff.loc[0, metId] -= gam_val
    for metId in ['adp_c', 'pi_c', 'h_c']:
        dfCoeff.loc[0, metId] += gam_val

    np.random.seed(2000)
    gamVal_rand = gam_val + gam_std*gam_val*normal(size=size)
    gamVal_rand = [round(coeff,7) if coeff > zero_tol else zero_tol for coeff in gamVal_rand] # Only allow positive value for GAM
    idx = range(1, size+1)
    for metId in ['atp_c', 'h2o_c']:
        dfCoeff.loc[idx, metId] -= gamVal_rand
    for metId in ['adp_c', 'pi_c', 'h_c']:
        dfCoeff.loc[idx, metId] += gamVal_rand

    atpm = ['atp_c', 'adp_c', 'h2o_c', 'pi_c', 'h_c']
    biomrxn = model.reactions.get_by_id(biomId)
    mets = atpm + sorted([met.id for met in biomrxn.metabolites.keys() if met.id not in atpm])
    dfCoeff = dfCoeff.reindex(['biomMW'] + mets, axis=1) 
        
    return dfCoeff

def make_pfba_csv(model, dfCoeff, dir_pfba, biomId):
    import cobra, sys, os
    
    if not os.path.exists(dir_pfba):
        os.makedirs(dir_pfba)
        
    try:
        model.solver = 'cplex'
    except:
        None
        
    infes_cases = []
    for i in dfCoeff.index:
        metDict = {k:dfCoeff.loc[i,k] for k in dfCoeff.columns[1:]}
        biomEqn = build_reaction_equation_from_metabolites_dict(metDict, arrow='-->')
        model.reactions.get_by_id(biomId).reaction = biomEqn

        try:
            pFba = cobra.flux_analysis.pfba(model, fraction_of_optimum=1)
        except:
            infes_cases.append(i)
            continue
        df = pd.DataFrame({'Rxn': pFba.fluxes.index.tolist(), 'Flux': pFba.fluxes.values.tolist()})
        df = df.loc[:, ['Rxn', 'Flux']]

        fname = 'pFBA' + str(i) + '.csv'
        df.to_csv(os.path.join(dir_pfba, fname), sep=',', index=None)
    try:
        if len(infes_cases) > 0:
            infes_cases_str = [str(case) for case in infes_cases]
            print('List of infeasible cases:', ','.join(infes_cases_str))
        else:
            print('No infeasible cases')
    except:
        print('Error in printing list of infeasible cases, need to inspect manually later')
        
    return None

def make_fluxes_dataframe(index, dir_pfba, zero_tol=1e-9):
    import os
    
    dfFlux = pd.DataFrame(index=index)
    for i in dfFlux.index:
        df = pd.read_csv(os.path.join(dir_pfba, 'pFBA' + str(i) + '.csv'), sep=',')
        df = df[df.Flux.abs() > zero_tol]
        df.index = df.Rxn.tolist()
        for rxn in df.index:
            dfFlux.loc[i, rxn] = df.loc[rxn, 'Flux']
        
    return dfFlux
