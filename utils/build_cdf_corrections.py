#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
This script reads in a DataFrame of paired VS-ST data with corrected T, 
and creates CDF-based RH-correction model.
'''
import numpy as np
import pandas as pd
import os, argparse, logging, csv, re, json
import matplotlib.pyplot as plt

__author__ = "Ting-Shuo Yo"
__copyright__ = "Copyright 2020~2022, DataQualia Lab Co. Ltd."
__credits__ = ["Ting-Shuo Yo"]
__license__ = "Apache License 2.0"
__version__ = "0.2.0"
__maintainer__ = "Ting-Shuo Yo"
__email__ = "tingyo@dataqualia.com"
__status__ = "development"
__date__ = '2023-07-25'


# Create CDF by sorting the values
def ecdf_sort(sample):
    ''' Empirical cumulative distribution function '''
    import numpy as np
    # convert sample to a numpy array, if it isn't already
    sample = np.atleast_1d(sample)
    # Sort the value and increment the CDF by 1/nrec
    quantiles = sample
    cumprob = np.arange(1, len(sample)+1)/len(sample)
    #
    return({'x':quantiles, 'cdf':cumprob})


# Create plot for the CDF
def plot_cdf_t(cdf_vs, cdf_st, dt=None, title=None):
    # Create figures
    fig, ax = plt.subplots(2, 1, figsize=(6, 6), sharex=True, gridspec_kw={'height_ratios': [5, 1]})
    # CDF
    ax[0].plot(cdf_vs['x'], cdf_vs['cdf'], c='k', ls='-', lw=0.8, label='T_vs')
    ax[0].plot(cdf_st['x'], cdf_st['cdf'], c='k', ls='--', lw=0.8, label='T_st')
    #ax[0].xlabel("Temperature (degree C)")
    ax[0].set_ylabel("CDF")
    ax[0].legend()    
    # Show title if specified
    if title is not None:
        ax[0].set_title(title)
    # Adjust values
    if dt is None:
        dt = cdf_vs['x'] - cdf_st['x']
    ax[1].plot(cdf_st['x'], dt, c='k', ls='-', lw=0.8)
    ax[1].set_xlabel("Temperature (degree C)")
    ax[1].set_ylabel("dT")
    # Display
    plt.show()
    return(0)

def plot_cdf_rh(cdf_vs, cdf_st, drh=None, title=None):
    # Create figures
    fig, ax = plt.subplots(2, 1, figsize=(6, 6), sharex=True, gridspec_kw={'height_ratios': [5, 1]})
    # CDF
    ax[0].plot(cdf_vs['x'], cdf_vs['cdf'], c='k', ls='-', lw=0.8, label='RH_vs')
    ax[0].plot(cdf_st['x'], cdf_st['cdf'], c='k', ls='--', lw=0.8, label='RH_st')
    #ax[0].xlabel("Temperature (degree C)")
    ax[0].set_ylabel("Cumulative Probability")
    ax[0].legend()    
    # Show title if specified
    if title is not None:
        ax[0].set_title(title)
    # Adjust values
    if drh is None:
        drh = cdf_vs['x'] - cdf_st['x']
    ax[1].plot(cdf_st['x'], drh, c='k', ls='-', lw=0.8)
    ax[1].set_xlabel("Relative Humidity (%)")
    ax[1].set_ylabel("dRH")
    # Display
    plt.show()
    return(0)

# CDF Model Building functions
def build_cdf_T(df):
    ''' For a given set of T measurement, calculate the CDF within (-80, 45). '''
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    # Get data
    t_vs_sorted = np.sort(df['t_vs'])
    t_st_sorted = np.sort(df['t_st'])
    # Remove upper and lower 2.5% of data
    nobs = df.shape[0]
    ndrop = np.floor(nobs*0.025).astype(int)
    t_vs_95 = t_vs_sorted[ndrop:-ndrop]
    t_st_95 = t_st_sorted[ndrop:-ndrop]
    cdf_vs_95 = ecdf_sort(t_vs_95)
    cdf_st_95 = ecdf_sort(t_st_95)
    # Adjust the range of data to (-80,45)
    xt = np.linspace(-80,45,num=1251)
    # Interpolation to (-80,45) by adding (0., 100.) to both end
    cdf_vs_intp = np.interp(xt, t_vs_95, cdf_vs_95['cdf'])
    cdf_st_intp = np.interp(xt, t_st_95, cdf_st_95['cdf'])
    cdf_vs_ext = {'x':xt,'cdf':cdf_vs_intp}
    cdf_st_ext = {'x':xt,'cdf':cdf_st_intp}
    dt_x = np.append(np.insert(t_st_95,0,-80),45)
    dt_y = np.append(np.insert(t_vs_95-t_st_95,0,0),0)
    dt_intp = np.interp(xt, dt_x, dt_y)
    cdf_dt_table = pd.DataFrame({
        't': np.round(xt,2),
        'dt': np.round(dt_intp,2)
    })
    # Evaluation
    dts = np.round(np.interp(df['t_st'], cdf_dt_table['t'], cdf_dt_table['dt']),2)
    t_st_adj = df['t_st'] + dts
    rmse_adj = np.sqrt(np.nanmean(np.square(df['t_vs'] - t_st_adj)))
    corr_adj = np.corrcoef(df['t_vs'], t_st_adj)[0,1]
    # Derive statistics
    rmse_paired = np.sqrt(np.nanmean(np.square(df['t_vs']-df['t_st'])))
    corr_paired = np.corrcoef(df['t_vs'], df['t_st'])[0,1]
    max_dt = (df['t_vs']-df['t_st']).max()
    min_dt = (df['t_vs']-df['t_st']).min()
    mean_dt = (df['t_vs']-df['t_st']).mean()
    median_dt = (df['t_vs']-df['t_st']).median()    
    rmse_sorted = np.sqrt(np.nanmean(np.square(t_vs_sorted-t_st_sorted)))
    corr_sorted = np.corrcoef(t_vs_sorted, t_st_sorted)[0,1]
    # Prepare output
    summary = {
        'num_obs':nobs,
        'dt_max': max_dt,
        'dt_min': min_dt,
        'dt_mean': mean_dt,
        'dt_med': median_dt,
        'RMSE_raw': rmse_paired,
        'CORR_raw': corr_paired,
        'RMSE_adj': rmse_adj,
        'CORR_adj': corr_adj,
        'STDEV_vs': df['t_vs'].std(),
        'STDEV_st': df['t_st'].std(),
        'STDEV_st_adj': t_st_adj.std(),
        'RMSE_sorted': rmse_sorted,
        'CORR_sorted': corr_sorted,
        'T_vs_95': ((t_vs_95[0],t_vs_95[-1])),
        'T_st_95': ((t_st_95[0],t_st_95[-1]))
    }
    return({'summary':summary, 'model':cdf_dt_table})

def build_cdf_RH(data):
    ''' For a given set of RH measurement, calculate the CDF within (0, 100). '''
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    # Clean up data
    df = data.loc[data['rh_st']<100, ['rh_vs', 'rh_st']].copy().dropna()
    # Get data
    rh_vs_sorted = np.sort(df['rh_vs'])
    rh_st_sorted = np.sort(df['rh_st'])
    # Remove upper and lower 2.5% of data
    nobs = df.shape[0]
    ndrop = np.floor(nobs*0.025).astype(int)
    rh_vs_95 = rh_vs_sorted[ndrop:-ndrop]
    rh_st_95 = rh_st_sorted[ndrop:-ndrop]
    cdf_vs_95 = ecdf_sort(rh_vs_95)
    cdf_st_95 = ecdf_sort(rh_st_95)
    #plot_cdf_rh(cdf_vs_95, cdf_st_95, title='Keep 95% of Data')
    # Adjust the range of data to (0,100)
    x = np.linspace(0,100,num=101)
    # Interpolation to (0,100) by adding (0., 100.) to both end
    cdf_vs_intp = np.interp(x, rh_vs_95, cdf_vs_95['cdf'])
    cdf_st_intp = np.interp(x, rh_st_95, cdf_st_95['cdf'])
    cdf_vs_ext = {'x':x,'cdf':cdf_vs_intp}
    cdf_st_ext = {'x':x,'cdf':cdf_st_intp}
    drh_x = np.append(np.insert(rh_st_95,0,0),100)
    drh_y = np.append(np.insert(rh_vs_95-rh_st_95,0,0),0)
    drh_intp = np.interp(x, drh_x, drh_y)
    #plot_cdf_rh(cdf_vs_ext, cdf_st_ext, drh=drh_intp, title='Extend to (0,100)')
    cdf_drh_table = pd.DataFrame({
        'rh': np.round(x,2),
        'drh': np.round(drh_intp,2)
    })
    # Evaluation
    drhs = np.round(np.interp(df['rh_st'], cdf_drh_table['rh'], cdf_drh_table['drh']),2)
    rh_st_adj = df['rh_st'] + drhs
    rmse_adj = np.sqrt(np.nanmean(np.square(df['rh_vs'] - rh_st_adj)))
    corr_adj = np.corrcoef(df['rh_vs'], rh_st_adj)[0,1]    
    # Derive statistics
    rmse_paired = np.sqrt(np.nanmean(np.square(df['rh_vs']-df['rh_st'])))
    corr_paired = np.corrcoef(df['rh_vs'], df['rh_st'])[0,1]
    max_drh = (df['rh_vs']-df['rh_st']).max()
    min_drh = (df['rh_vs']-df['rh_st']).min()
    mean_drh = (df['rh_vs']-df['rh_st']).mean()
    median_drh = (df['rh_vs']-df['rh_st']).median()    
    rmse_sorted = np.sqrt(np.nanmean(np.square(rh_vs_sorted-rh_st_sorted)))
    corr_sorted = np.corrcoef(rh_vs_sorted, rh_st_sorted)[0,1]
    # Prepare output
    summary = {
        'num_obs':nobs,
        'dt_max': max_drh,
        'dt_min': min_drh,
        'dt_mean': mean_drh,
        'dt_med': median_drh,
        'RMSE_raw': rmse_paired,
        'CORR_raw': corr_paired,
        'RMSE_adj': rmse_adj,
        'CORR_adj': corr_adj,
        'STDEV_vs': df['rh_vs'].std(),
        'STDEV_st': df['rh_st'].std(),
        'STDEV_st_adj': rh_st_adj.std(),
        'RMSE_sorted': rmse_sorted,
        'CORR_sorted': corr_sorted,
        'RH_vs_95': (rh_vs_95[0],rh_vs_95[-1]),
        'RH_st_95': (rh_st_95[0],rh_st_95[-1])
    }
    return({'summary':summary, 'model':cdf_drh_table})

def build_cdf_correction_by_pressure(df, p_mid, p_bins):
    ''' Build CDF-based correction table by pressure bins. '''
    output = []
    model = None
    # Run through bins
    for i in range(len(p_mid)):
        # Select bin
        in_bin = (df['p_st_adj']>p_bins[i][0])&(df['p_st_adj']<=p_bins[i][1])
        bindf = df.loc[in_bin,['t_vs','t_st']].dropna()
        if bindf.shape[0]>1000:
            tmp = build_cdf_T(bindf)
            # Aggregate statistics
            tmp['summary']['group'] = 'P-'+str(p_mid[i])
            output.append(tmp['summary'])
            # Aggregate model
            tmp['model'].columns = ['t', 'P-'+str(p_mid[i])]
            if model is None:
                model = tmp['model']
            else:
                model = model.merge(tmp['model'], left_on='t', right_on='t')
        else:
            logging.warning('Pressure bin '+str(p_mid[i])+' has <1000 records, skipped.')
            output.append({
                'num_obs':bindf.shape[0],
                'dt_max': None,
                'dt_min': None,
                'dt_mean': None,
                'dt_med': None,
                'RMSE_raw': None,
                'CORR_raw': None,
                'RMSE_adj': None,
                'CORR_adj': None,
                'STDEV_vs': None,
                'STDEV_st': None,
                'STDEV_st_adj': None,
                'RMSE_sorted': None,
                'CORR_sorted': None,
                'T_vs_95': None,
                'T_st_95': None,
                'group': 'P-'+str(p_mid[i])
            })

    # Output model and summary
    output = pd.DataFrame(output)
    return({'summary':output, 'model':model})

def build_cdf_correction_by_temperature(df, t_mid, t_bins):
    ''' Build CDF-based correction table by pressure bins. '''
    output = []
    model = None
    # Run through bins
    for i in range(len(t_mid)):
        # Select bin
        in_bin = (df['t_st_adj']>t_bins[i][0])&(df['t_st_adj']<=t_bins[i][1])
        bindf = df.loc[in_bin,['rh_vs','rh_st']].dropna()
        if bindf.shape[0]>1000:
            tmp = build_cdf_RH(bindf)
            # Aggregate statistics
            tmp['summary']['group'] = 'T-'+str(t_mid[i])
            output.append(tmp['summary'])
            # Aggregate model
            tmp['model'].columns = ['rh', 'T-'+str(t_mid[i])]
            if model is None:
                model = tmp['model']
            else:
                model = model.merge(tmp['model'], left_on='rh', right_on='rh')
        else:
            logging.warning('Temperature bin '+str(t_mid[i])+' has <1000 records, skipped.')
            tmpmodel = pd.DataFrame({
                'rh': np.linspace(0,100,num=101),
                'drh': np.zeros(101)              
                })
            # Aggregate model
            tmpmodel.columns = ['rh', 'T-'+str(t_mid[i])]
            if model is None:
                model = tmpmodel
            else:
                model = model.merge(tmpmodel, left_on='rh', right_on='rh')
            # Skip evaluation
            output.append({
                'num_obs':bindf.shape[0],
                'dt_max': None,
                'dt_min': None,
                'dt_mean': None,
                'dt_med': None,
                'RMSE_raw': None,
                'CORR_raw': None,
                'RMSE_adj': None,
                'CORR_adj': None,
                'STDEV_vs': None,
                'STDEV_st': None,
                'STDEV_st_adj': None,
                'RMSE_sorted': None,
                'CORR_sorted': None,
                'RH_vs_95': None,
                'RH_st_95': None,
                'group': 'T-'+str(t_mid[i])
            })

    # Output model and summary
    output = pd.DataFrame(output)
    return({'summary':output, 'model':model})

# Perform CDF-based correction
def correct_T_by_CDF(data, model, p_col='p_st_adj', t_col='t_st', model_surfix='day'):
    ''' Correct T_st with CDF model '''
    import numpy as np
    import pandas as pd
    # Parameters
    p_mid = np.arange(200, 1030, 50)
    p_bins = [(p-25, p+25) for p in p_mid]
    # Copy the P and T from the data source
    df = data.loc[:,['p_st_adj','t_vs','t_st','hour']]
    df['dt'] = np.nan
    # Loop through p_bins
    for i in range(len(p_mid)):
        # Select day bin
        in_bin = (df[p_col]>p_bins[i][0])&(df[p_col]<=p_bins[i][1])&(df['hour']<12)
        model_surfix='day'
        bindf = df.loc[in_bin,:]
        # Select correction table and look for corresponding correction
        bin_col = 'P-'+str(p_mid[i])+'_'+model_surfix
        dt = np.round(np.interp(bindf['t_st'], model['t'], model[bin_col]),2)
        df.loc[in_bin,'dt'] = dt
        # Select night bin
        in_bin = (df[p_col]>p_bins[i][0])&(df[p_col]<=p_bins[i][1])&(df['hour']>=12)
        model_surfix='night'
        bindf = df.loc[in_bin,:]
        # Select correction table and look for corresponding correction
        bin_col = 'P-'+str(p_mid[i])+'_'+model_surfix
        dt = np.round(np.interp(bindf['t_st'], model['t'], model[bin_col]),2)
        df.loc[in_bin,'dt'] = dt    # End loop check
    #
    not_corrected = df['dt'].isnull().sum()
    print('Number of T not corrected: '+str(not_corrected))
    df['t_st_adj'] = df['t_st'] + df['dt']
    # Calculate metrics
    metrics = {
        'rmse_raw': np.sqrt(np.nanmean(np.square(data['t_vs']-data['t_st']))),
        'rmse_adj': np.sqrt(np.nanmean(np.square(data['t_vs']-df['t_st_adj']))),
        'corr_raw': df.dropna().corr().loc['t_vs','t_st'],
        'corr_adj': df.dropna().corr().loc['t_vs','t_st_adj'],
        'sd_vs': data['t_vs'].std(),
        'sd_st': data['t_st'].std(),
        'sd_st_adj': df['t_st_adj'].std()        
    }
    # 
    return({'metrics':metrics, 'data':df})

def correct_RH_by_CDF(data, model, t_col='t_st_adj', rh_col='rh_st'):
    ''' Correct RH_st with CDF model '''
    import numpy as np
    import pandas as pd
    # Parameters
    t_mid = np.arange(-60, 31, 10)
    t_bins = [(t-5, t+5) for t in t_mid]
    # Copy the P and T from the data source
    df = data.loc[:,['t_st_adj','rh_vs','rh_st','hour']]
    df['drh'] = np.nan
    # Loop through p_bins
    for i in range(len(t_mid)):
        # Select day bin
        in_bin = (df['t_st_adj']>t_bins[i][0])&(df['t_st_adj']<=t_bins[i][1])&(df['hour']<12)
        model_surfix = 'day'
        bindf = df.loc[in_bin,:]
        # Select correction table and look for corresponding correction
        bin_col = 'T-'+str(t_mid[i])+'_'+model_surfix
        drh = np.round(np.interp(bindf['rh_st'], model['rh'], model[bin_col]),2)
        df.loc[in_bin,'drh'] = drh
        # Select night bin
        in_bin = (df['t_st_adj']>t_bins[i][0])&(df['t_st_adj']<=t_bins[i][1])&(df['hour']>=12)
        model_surfix = 'night'
        bindf = df.loc[in_bin,:]
        # Select correction table and look for corresponding correction
        bin_col = 'T-'+str(t_mid[i])+'_'+model_surfix
        drh = np.round(np.interp(bindf['rh_st'], model['rh'], model[bin_col]),2)
        df.loc[in_bin,'drh'] = drh
    # End loop check
    not_corrected = df['drh'].isnull().sum()
    print('Number of RH not corrected: '+str(not_corrected))
    df['rh_st_adj'] = df['rh_st'] + df['drh']
    # Calculate metrics
    metrics = {
        'rmse_raw': np.sqrt(np.nanmean(np.square(data['rh_vs']-data['rh_st']))),
        'rmse_adj': np.sqrt(np.nanmean(np.square(data['rh_vs']-df['rh_st_adj']))),
        'corr_raw': df.dropna().corr().loc['rh_vs','rh_st'],
        'corr_adj': df.dropna().corr().loc['rh_vs','rh_st_adj'],
        'sd_vs': data['rh_vs'].std(),
        'sd_st': data['rh_st'].std(),
        'sd_st_adj': df['rh_st_adj'].std()        
    }
    # 
    return({'metrics':metrics, 'data':df})

# Wrapper of T correction
def build_and_evaluate_CDF_T(data):
    # Parameters
    p_mid = np.arange(200, 1030, 50)
    p_bins = [(p-25, p+25) for p in p_mid]
    # Build CDF model for all data
    logging.info('[T] Building CDF-based correction table for the whole dataset.')
    df_t = data.loc[:,['t_vs','t_st']].dropna()
    model_all = build_cdf_T(df_t)
    model_all['summary']['group'] = 'All'
    summary = [pd.DataFrame(model_all['summary'])]
    # Build CDF model of day (hour<12Z)
    logging.info('[T] Building CDF-based correction table for day time by pressure levels.')
    df_day = data.loc[data['hour']<12,:]
    model_day = build_cdf_correction_by_pressure(df_day, p_mid, p_bins)
    model_day['summary']['group'] = [g+'-day' for g in model_day['summary']['group']]
    summary.append(model_day['summary'])
    # Build CDF model of night (hour>=12Z)
    logging.info('[T] Building CDF-based correction table for night time by pressure levels.')
    df_night = data.loc[data['hour']>=12,:]
    model_night = build_cdf_correction_by_pressure(df_night, p_mid, p_bins)
    model_night['summary']['group'] = [g+'-night' for g in model_night['summary']['group']]
    summary.append(model_night['summary'])    
    # Output model and summary
    model = model_day['model'].merge(model_night['model'], left_on='t', right_on='t', suffixes=('_day', '_night'))
    summary = pd.concat(summary)
    return({'model':model, 'summary':summary})

# Wrapper of RH correction
def build_and_evaluate_CDF_RH(data):
    # Parameters
    t_mid = np.arange(-60, 31, 10)
    t_bins = [(t-5, t+5) for t in t_mid]
    # Build CDF model for all data
    logging.info('[RH] Building CDF-based correction table for the whole dataset.')
    df_rh = data.loc[:,['rh_vs','rh_st']].dropna()
    model_all = build_cdf_RH(df_rh)
    model_all['summary']['group'] = 'All'
    summary = [pd.DataFrame(model_all['summary'])]
    # Build CDF model of day (hour<12Z)
    logging.info('[RH] Building CDF-based correction table for day time by temperature levels.')
    df_day = data.loc[data['hour']<12,:]
    model_day = build_cdf_correction_by_temperature(df_day, t_mid, t_bins)
    model_day['summary']['group'] = [g+'-day' for g in model_day['summary']['group']]
    summary.append(model_day['summary'])
    # Build CDF model of night (hour>=12Z)
    logging.info('[RH] Building CDF-based correction table for night time by temperature levels.')
    df_night = data.loc[data['hour']>=12,:]
    model_night = build_cdf_correction_by_temperature(df_night, t_mid, t_bins)
    model_night['summary']['group'] = [g+'-night' for g in model_night['summary']['group']]
    summary.append(model_night['summary'])    
    # Output model and summary
    model = model_day['model'].merge(model_night['model'], left_on='rh', right_on='rh', suffixes=('_day', '_night'))
    summary = pd.concat(summary)
    return({'model':model, 'summary':summary})

### Main Script ###
#-----------------------------------------------------------------------
def main():
    # Configure Argument Parser
    parser = argparse.ArgumentParser(description='Evaluate the performance of GLM with various feature vectors and events.')
    parser.add_argument('--datapath', '-i', help='the paired co-launch data file.')
    parser.add_argument('--output', '-o', default='./', help='the directory for output files.')
    parser.add_argument('--logfile', '-l', default=None, help='the log file.')
    args = parser.parse_args()
    # Set up logging
    if not args.logfile is None:
        logging.basicConfig(level=logging.DEBUG, filename=args.logfile, filemode='w')
    else:
        logging.basicConfig(level=logging.DEBUG)
    logging.debug(args)
    # Load configuration file if specified
    data = pd.read_pickle(args.datapath)
    # Step 1: Build T correction model
    cdf_T = build_and_evaluate_CDF_T(data)
    # Step 2: Correct T
    results_T = correct_T_by_CDF(data, cdf_T['model'])
    data['t_st_adj'] = results_T['data']['t_st_adj']
    # Step 3: Build RH correction model
    cdf_RH = build_and_evaluate_CDF_RH(data)
    # Step 4: Correct RH
    results_RH = correct_RH_by_CDF(data, cdf_RH['model'])
    data['rh_st_adj'] = results_RH['data']['rh_st_adj']
    # Print out results summary
    tmp = pd.DataFrame([results_T['metrics'], results_RH['metrics']])
    tmp.index = ['T','RH']
    print(tmp)
    # Save output
    cdf_T['model'].to_csv(args.output+'/cdf_T_correction.csv', index=False)
    cdf_RH['model'].to_csv(args.output+'/cdf_RH_correction.csv', index=False)
    cdf_T['summary'].to_csv(args.output+'/summary_cdf_T_correction.csv', index=False)
    cdf_RH['summary'].to_csv(args.output+'/summary_cdf_RH_correction.csv', index=False)
    # Done
    return(0)
#==========
# Script
#==========
if __name__=="__main__":
    main()

