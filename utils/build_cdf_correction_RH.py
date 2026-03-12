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
        'dt_medx': median_drh,
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
                'RH_vs_95': None,
                'RH_st_95': None,
                'group': 'T-'+str(t_mid[i])
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
    t_mid = np.arange(-60, 31, 10)
    t_bins = [(t-5, t+5) for t in t_mid]
    # Build CDF model for all data
    logging.info('Building CDF-based correction table for the whole dataset.')
    df_rh = data.loc[:,['rh_vs','rh_st']].dropna()
    model_all = build_cdf_RH(df_rh)
    model_all['summary']['group'] = 'All'
    summary = [pd.DataFrame(model_all['summary'])]
    # Build CDF model of day (hour<12Z)
    logging.info('Building CDF-based correction table for day time by temperature levels.')
    df_day = data.loc[data['hour']<12,:]
    model_day = build_cdf_correction_by_temperature(df_day, t_mid, t_bins)
    model_day['summary']['group'] = [g+'-day' for g in model_day['summary']['group']]
    summary.append(model_day['summary'])
    # Build CDF model of night (hour>=12Z)
    logging.info('Building CDF-based correction table for night time by temperature levels.')
    df_night = data.loc[data['hour']>=12,:]
    model_night = build_cdf_correction_by_temperature(df_night, t_mid, t_bins)
    model_night['summary']['group'] = [g+'-night' for g in model_night['summary']['group']]
    summary.append(model_night['summary'])    
    # Output model and summary
    model = model_day['model'].merge(model_night['model'], left_on='rh', right_on='rh', suffixes=('_day', '_night'))
    model.to_csv(args.output+'/cdf_RH_correction.csv', index=False)
    summary = pd.concat(summary)
    summary.to_csv(args.output+'/summary_cdf_RH_correction.csv', index=False)
    # Done
    return(0)
#==========
# Script
#==========
if __name__=="__main__":
    main()

