#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
This script reads in a DataFrame of paired VS-ST data, and creates 
CDF-based T-correction model.
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
        'dt_medx': median_dt,
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
                'dt_medx': None,
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
    # Parameters
    p_mid = np.arange(200, 1030, 50)
    p_bins = [(p-25, p+25) for p in p_mid]
    # Build CDF model for all data
    logging.info('Building CDF-based correction table for the whole dataset.')
    df_t = data.loc[:,['t_vs','t_st']].dropna()
    model_all = build_cdf_T(df_t)
    model_all['summary']['group'] = 'All'
    summary = [pd.DataFrame(model_all['summary'])]
    # Build CDF model of day (hour<12Z)
    logging.info('Building CDF-based correction table for day time by pressure levels.')
    df_day = data.loc[data['hour']<12,:]
    model_day = build_cdf_correction_by_pressure(df_day, p_mid, p_bins)
    model_day['summary']['group'] = [g+'-day' for g in model_day['summary']['group']]
    summary.append(model_day['summary'])
    # Build CDF model of night (hour>=12Z)
    logging.info('Building CDF-based correction table for night time by pressure levels.')
    df_night = data.loc[data['hour']>=12,:]
    model_night = build_cdf_correction_by_pressure(df_night, p_mid, p_bins)
    model_night['summary']['group'] = [g+'-night' for g in model_night['summary']['group']]
    summary.append(model_night['summary'])    
    # Output model and summary
    model = model_day['model'].merge(model_night['model'], left_on='t', right_on='t', suffixes=('_day', '_night'))
    model.to_csv(args.output+'/cdf_T_correction.csv', index=False)
    summary = pd.concat(summary)
    summary.to_csv(args.output+'/summary_cdf_T_correction.csv', index=False)
    # Done
    return(0)
#==========
# Script
#==========
if __name__=="__main__":
    main()

