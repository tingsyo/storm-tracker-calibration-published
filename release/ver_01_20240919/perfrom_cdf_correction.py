#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
This script takes in a ST sounding profile and correct T and RH
with trained CDF models.
'''
import numpy as np
import pandas as pd
import os, argparse, logging, csv, re, json
from datetime import datetime, timedelta

__author__ = "Ting-Shuo Yo"
__copyright__ = "Copyright 2022~2024, National Taiwan University."
__credits__ = ["Ting-Shuo Yo"]
__license__ = "Apache License 2.0"
__version__ = "0.8.1"
__maintainer__ = "Ting-Shuo Yo"
__email__ = "tsyo@ntu.edu.tw"
__status__ = "development"
__date__ = '2024-09-19'

def correct_T_by_CDF(df, model, model_surfix='day'):
    ''' Correct T_st with CDF model '''
    # Parameters
    p_mid = np.arange(200, 1030, 50)
    p_bins = [(p-25, p+25) for p in p_mid]
    # Copy the P and T from the data source
    df['dt'] = np.nan
    # Loop through p_bins
    for i in range(len(p_mid)):
        # Select bin
        in_bin = (df['p_st_adj']>p_bins[i][0])&(df['p_st_adj']<=p_bins[i][1])
        bindf = df.loc[in_bin,:]
        # Select correction table and look for corresponding correction
        bin_col = 'P-'+str(p_mid[i])+'_'+model_surfix
        dt = np.round(np.interp(bindf['t_st'], model['t'], model[bin_col]),2)
        df.loc[in_bin,'dt'] = dt
    # End loop check
    not_corrected = df['dt'].isnull().sum()
    logging.info('Number of T not corrected: '+str(not_corrected))
    df['t_st_adj'] = df['t_st'] + df['dt']
    # 
    return(df['t_st_adj'])


def correct_RH_by_CDF(df, model, model_surfix='day'):
    ''' Correct RH_st with CDF model '''
    # Parameters
    t_mid = np.arange(-60, 31, 10)
    t_bins = [(t-5, t+5) for t in t_mid]
    # Copy the P and T from the data source
    df['drh'] = np.nan
    # Loop through p_bins
    for i in range(len(t_mid)):
        # Select bin
        in_bin = (df['t_st_adj']>t_bins[i][0])&(df['t_st_adj']<=t_bins[i][1])
        bindf = df.loc[in_bin,:]
        # Select correction table and look for corresponding correction
        bin_col = 'T-'+str(t_mid[i])+'_'+model_surfix
        drh = np.round(np.interp(bindf['rh_st'], model['rh'], model[bin_col]),2)
        df.loc[in_bin,'drh'] = drh
    # End loop check
    not_corrected = df['drh'].isnull().sum()
    logging.info('Number of RH not corrected: '+str(not_corrected))
    df['rh_st_adj'] = df['rh_st'] + df['drh']
    # Done
    return(df['rh_st_adj'])

def evluate_corrections(data):
    # Calculate metrics
    metrics = {
        'mean_T':  data['T'].mean(),
        'sd_T': data['T'].std(),
        'mean_T_adj':  data['T_adj'].mean(),
        'sd_T_adj': data['T_adj'].std(),
        #
        'mean_RH':  data['RH'].mean(),
        'sd_RH': data['RH'].std(),
        'mean_RH_adj':  data['RH_adj'].mean(),
        'sd_RH_adj': data['RH_adj'].std()
    }
    # Done
    return(metrics)

def correct_ST_data(data, hour, dp0, t_model, rh_model):
    # Determine day/night model
    if hour<12:
        model_suffix = 'day'
    else:
        model_suffix = 'night'
    # Correct P
    if dp0 is None:
        p_adj = data['P']
    else:
        p_adj = data['P'] - dp0
    # Correct T
    dft =  pd.DataFrame({'t_st': data['T'].copy(), 'p_st_adj': p_adj})
    t_adj = correct_T_by_CDF(dft, t_model, model_suffix)
    # Correct RH
    dfrh = pd.DataFrame({'rh_st': data['RH'].copy(), 't_st_adj': t_adj})
    rh_adj = correct_RH_by_CDF(dfrh, rh_model, model_suffix)
    # Return corrected data
    data['T_adj'] = t_adj
    data['RH_adj'] = rh_adj
    metrics = evluate_corrections(data)
    # 
    return({'data':data, 'metrics':metrics})


### Main Script ###
#-----------------------------------------------------------------------
def main():
    # Configure Argument Parser
    parser = argparse.ArgumentParser(description='Correct the ST profile (T/RH) with the CDF models.')
    parser.add_argument('--input', '-i', help='the ST observation file in csv format.')
    parser.add_argument('--modelt', default='cdf_T_correction.csv', help='the CDF-T model file.')
    parser.add_argument('--modelrh', default='cdf_RH_correction.csv', help='the CDF-RH model file.')
    parser.add_argument('--output', '-o', default='./', help='the directory for output files.')
    parser.add_argument('--logfile', '-l', default=None, help='the log file.')
    parser.add_argument('--obstime', '-t', default=None, help='the time of observation in YYYYMMDDHH (UTC).')
    parser.add_argument('--dp0', type=float, default=None, help='the ground check reference pressure (hPa).')
    args = parser.parse_args()
    # Set up logging
    if not args.logfile is None:
        logging.basicConfig(level=logging.DEBUG, filename=args.logfile, filemode='w')
    else:
        logging.basicConfig(level=logging.DEBUG)
    logging.debug(args)
    # Load CDF models
    t_model = pd.read_csv(args.modelt)
    rh_model = pd.read_csv(args.modelrh)
    # Parse observation time
    obstime = datetime.strptime(args.obstime, '%Y%m%d%H')
    julianday = int(obstime.strftime('%j'))
    hour = obstime.hour
    # Load data
    tmp = pd.read_csv(args.input)#, dtype='str')
    # correct ST data
    results = correct_ST_data(tmp, hour, args.dp0, t_model, rh_model)
    fname = os.path.basename(args.input)
    output_file = args.output+fname.replace('.csv','_cdf_corrected.csv')
    results['data'].to_csv(output_file, index=False)
    logging.info('Corrected data is stored in '+output_file)
    # Print metrics
    results['metrics']['obs'] = obstime
    results['metrics']['hour'] = hour
    results['metrics']['julianday'] = julianday
    results['metrics']['nrec'] = tmp.shape[0]
    results['metrics']['n_complete_rec'] = tmp.dropna().shape[0]
    logging.info(results['metrics'])
    #metric_file = args.output+fname.replace('.csv','_metrics.csv')
    #pd.DataFrame(results['metrics']).to_csv(metric_file, index=False)
    # Done
    return(0)
#==========
# Script
#==========
if __name__=="__main__":
    main()

