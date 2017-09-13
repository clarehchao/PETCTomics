#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 5/23/17

@author: shuang
@goal: categorize her2 clinical data for data explore

"""

import pandas as pd
import re
import numpy as np
import glob
import os

def str_find_dict(the_str,the_dict):
    for kk in the_dict.keys():
        ss = the_str.lower().replace(' ','')  #lower-case and remove white space
        match_list = the_dict[kk][0]
        for mm in match_list:
            # print 'the_str:{}, match: {}'.format(ss,mm)
            if ss.find(mm.lower()) >= 0:
                return kk
    print '{}: cannot find a key in the dictionary that matches the string'.format(the_str)
    return 0

def recur_func2(x):
    """
    out of the patients who was disease free in his/her life, 0 if no recur, 1 if recur
    :param x: entry in the recurrent type summary data
    :return: return nan if N/A, 0 if disease free/no recur, 1 if some sort of recurrence or never disease free
    """
    if x is np.nan or x.lower() == 'unknown' or x.lower() == 'never disease free':
        return np.nan
    elif x.lower().find('recur') >= 0:
        return 1
    elif x == 'NONE/DISEASE FREE':
        return 0

def recur_func1(x):
    if x is np.nan or x.lower() == 'unknown':
        return np.nan
    if x.lower() == 'never disease free':
        return 2
    elif x == 'NONE/DISEASE FREE':
        return 0
    else:
        return 1

def recur_func0(x):
    if x == 'NEVER DISEASE FREE':
        return 2
    # elif x == 'NONE/DISEASE FREE' or x is np.nan or x.lower() == 'unknown':
    elif x == 'NONE/DISEASE FREE' or x is np.nan:
        return 0
    else:
        return 1

def recur_func3(row):
    if row['Recurrence'] is np.nan or row['Recurrence'] is None:
        cr_recur = row['Recurrence Type Summary']
        if cr_recur is np.nan or cr_recur.lower() == 'unknown':
            print 'use CR data: {},{}'.format(row['Recurrence'], cr_recur)
            return np.nan
        if cr_recur.lower() == 'never disease free':
            return 2
        elif cr_recur == 'NONE/DISEASE FREE':
            return 0
        else:
            return 1
        # print 'use CR data: {},{},{}'.format(row['Recurrence'], cr_recur, val)
    else:
        return int(row['Recurrence'])


def recur_dT(row):
    dT = pd.to_datetime(row['Last contact/death']) - pd.to_datetime(row['Date of Primary Dx'])
    if isinstance(dT, pd.tslib.NaTType):
        dT_yr = -1
    else:
        dT_yr = dT.days / 365.

    return dT_yr

def diseasefree5yr_func(row):
    """
    :param x: a row in a pandas dataframe
    :return: return nan if N/A, 0 if disease free/no recur after 5 years, 1 if some sort of recurrence or never disease free
    """
    dT = pd.to_datetime(row['Last contact/death']) - pd.to_datetime(row['Date of Primary Dx'])
    if isinstance(dT,pd.tslib.NaTType):
        dT_yr = -1
    else:
        dT_yr = dT.days/365.

    if row['Recurrence Type Summary'] is np.nan or row['Recurrence Type Summary'].lower() == 'unknown':
        return np.nan
    elif row['Recurrence Type Summary'] == 'NONE/DISEASE FREE' and dT_yr >= 5.0:
        return 1
    else:
        return 0

# 0: no recur, 1: local recur, 2: distant recur (bone and others),
# nan or 'unknown' or 'never disease free' or 'recurred, type unknow': N/A,
def recur_type_func(x):
    if x is np.nan or x.lower().find('unknown') >= 0 or x.lower() == 'never disease free':
        return np.nan
    elif x == 'NONE/DISEASE FREE':
        return 0
    elif re.search(r'LOCAL RECUR[\w\s.]*',x):
        return 1
    elif re.search(r'DIST RECUR[\w\s.]*',x):
        return 2

def BoneRecurOrNot_func(x):
    if x is np.nan or x.lower().find('unknown') >= 0 or x.lower() == 'never disease free':
        return np.nan
    elif re.search(r'[\w\s.]BONE',x):
        return 1
    else:
        return 0

# def recur_type_func(x):
#     if x is np.nan or x.lower() == 'unknown' or x.lower() == 'never disease free':
#         return np.nan
#     elif x == 'NONE/DISEASE FREE':
#         return 0
#     elif x == 'DIST RECUR, BONE':
#         return 1
#     elif re.search(r'LOCAL RECUR[\w\s.]*',x):
#         return 2
#     elif re.search(r'DIST RECUR[\w\s.]*',x) and x != 'DIST RECUR, BONE':
#         return 3

def Tstage_func(x):
    tstage_dict = {'T0':(['T0','TX','calc','Tis'],0), 'T1':(['T1'],1), 'T2':(['T2'],2),'T3':(['T3'],3),'T4':(['T4'],4)}
    # add the numerical values to each or combined category
    # tstage_k = tstage_dict.keys()
    # for i in range(len(tstage_k)):
    #     tstage_dict[tstage_k[i]] = (tstage_dict[tstage_k[i]],i)

    if x is np.nan:
        # print 'x is nan'
        return x
    else:
        kk = str_find_dict(x,tstage_dict)
        if kk != 0:
            return tstage_dict[kk][1]
        else:
            return np.nan


def Nstage_func(x):
    nstage_dict = {'N0': (['N0','NX','calc'],0), 'N1': (['N1'],1), 'N2': (['N2'],2), 'N3': (['N3'],3)}

    # # add the numerical values to each or combined category
    # nstage_k = nstage_dict.keys()
    # for i in range(len(nstage_k)):
    #     nstage_dict[nstage_k[i]] = (nstage_dict[nstage_k[i]], i)

    if x is np.nan:
        return x
    else:
        kk = str_find_dict(x, nstage_dict)
        if kk != 0:
            return nstage_dict[kk][1]
        else:
            return np.nan

def Overallstage_func(x):
    ostage_dict = {'0': (['0','calc'],0), '1': (['1'],1), '2': (['2'],2), '3': (['3'],3), '4':(['4'],4)}

    # # add the numerical values to each or combined category
    # ostage_k = ostage_dict.keys()
    # for i in range(len(ostage_k)):
    #     ostage_dict[ostage_k[i]] = (ostage_dict[ostage_k[i]], i)

    if x is np.nan:
        return x
    else:
        kk = str_find_dict(x, ostage_dict)
        if kk != 0:
            return ostage_dict[kk][1]
        else:
            return np.nan

def Overallstage_func2(x):
    ostage_dict = {'0': (['0','calc'],0), '1': (['1a','1b','2a'],1), '2': (['2b'],2), '3': (['3'],3), '4':(['4'],4)}

    # # add the numerical values to each or combined category
    # ostage_k = ostage_dict.keys()
    # for i in range(len(ostage_k)):
    #     ostage_dict[ostage_k[i]] = (ostage_dict[ostage_k[i]], i)

    if x is np.nan:
        return x
    else:
        kk = str_find_dict(x, ostage_dict)
        if kk != 0:
            return ostage_dict[kk][1]
        else:
            return np.nan

def BC_subtype_func(row):
# split the tumors into tumor subtypes (suggested by E. Jones)
# 0: HR+/HER2-, 1: HR+/HER2+, 2: HR-/HER2+, 3: triple neg, HR+ is ER+ or PR+
    is_HR_pos = (row['Sjoerd_PR'] == 1) or (row['Sjoerd_ER'] == 1)
    if is_HR_pos and row['Sjoerd_HER2'] == 0:
        return 0
    elif is_HR_pos and row['Sjoerd_HER2'] == 1:
        return 1
    elif ~is_HR_pos and row['Sjoerd_HER2'] == 1:
        return 2
    elif ~is_HR_pos and row['Sjoerd_HER2'] == 0:
        return 3
    else:
        return np.nan

def custom_median(df):
    # inclue NaN in the median calculation
    return df.median(skipna=False)

def Min_Ncluster_CC(cc_dir):
    all_cc = glob.glob('{}/ConsensusClass*'.format(cc_dir))
    dict_Nmincluster = {}
    for ff in all_cc:
        check = re.search('ConsensusClass_kmax(\w+)',os.path.basename(ff))
        if check:
            ncluster = check.group(1)
            df_cs_class = pd.read_csv(ff)
            df_cs_class.columns = ['ptid_side', 'cs_class']

            # find the N of the smallest cluster in each cluster result for future analysis
            N_mincluster = df_cs_class['cs_class'].value_counts().min()
            dict_Nmincluster[ncluster] = N_mincluster
        else:
            print 'cannot determine Ncluster'
            continue

    df_Nmincluster = pd.DataFrame(dict_Nmincluster.items(), columns=['k','N_mincluster'])
    df_Nmincluster['k'] = pd.to_numeric(df_Nmincluster['k'])
    df_Nmincluster['N_mincluster'] = pd.to_numeric(df_Nmincluster['N_mincluster'])
    return df_Nmincluster


def Median_Cluster_Consensus(the_dir, dir_prefix):
    all_dir = [d for d in glob.glob('{}/{}*'.format(the_dir, dir_prefix)) if os.path.isdir(d)]

    df_medCC_all = pd.DataFrame()
    for dd in all_dir:
        # deteremine the cluster method and etc.
        tmp = os.path.basename(dd).split('_')
        if len(tmp) == 3:
            dummy, cm, dm = tmp
        else:
            dummy, cm, cl, dm = tmp
        cc_fname = '{}/ClusterConsensus.csv'.format(dd)
        df_cc = pd.read_csv(cc_fname)
        df_Nmincluster = Min_Ncluster_CC(dd)

        check_nan = df_cc.ix[df_cc['clusterConsensus'].isnull(), 'k']
        if not check_nan.empty:
            check = check_nan.unique().min()
            if check > 5.0:
                mc = df_cc.groupby('k').agg({'clusterConsensus': custom_median})
                mc = mc.reset_index()
                mc['cluster_method'] = cm
                if len(tmp) == 3:
                    mc['cluster_method'] = cm
                    mc['cluster_linkage'] = np.nan
                    mc['dist_method'] = dm
                else:
                    mc['cluster_method'] = cm
                    mc['cluster_linkage'] = cl
                    mc['dist_method'] = dm
                mc.rename(columns={'clusterConsensus': 'medianCC'}, inplace=True)
                mc = pd.merge(mc,df_Nmincluster,how='left',on='k')
                df_medCC_all = pd.concat([df_medCC_all, mc], ignore_index=True)
            else:
                continue
        else:
            mc = df_cc.groupby('k').agg({'clusterConsensus': custom_median})
            mc = mc.reset_index()
            if len(tmp) == 3:
                mc['cluster_method'] = cm
                mc['cluster_linkage'] = np.nan
                mc['dist_method'] = dm
            else:
                mc['cluster_method'] = cm
                mc['cluster_linkage'] = cl
                mc['dist_method'] = dm
            mc.rename(columns={'clusterConsensus': 'medianCC'}, inplace=True)
            mc = pd.merge(mc, df_Nmincluster, how='left', on='k')
            df_medCC_all = pd.concat([df_medCC_all, mc], ignore_index=True)

    # print out the cluster setting for each Ncluster with the highest cluster consensus
    df_tmp = df_medCC_all[df_medCC_all['N_mincluster'] > 5]
    idx = df_tmp.groupby('k').apply(lambda df: df.medianCC.argmax())
    df_oi = df_tmp.ix[idx,['k','N_mincluster','cluster_method','cluster_linkage','dist_method','medianCC']]

    return df_medCC_all, df_oi

def func_feat_rank(the_df, df_info):
    # convert each series to a dataframe and do a join with df_info to do necessary operation to return the values needed
    the_df = pd.DataFrame(the_df)
    the_df.reset_index(inplace=True)
    col = the_df.columns
    col_val = [ss for ss in col if ss is not 'index'][0]
    the_df.rename(columns={'index':'feature', col_val:'values'}, inplace=True)
    jdf = pd.merge(the_df, df_info, on='feature')
    jdf['rank_metric'] = jdf['values'] * jdf['feature_cs_rank']

    # srm = jdf.groupby('feature_cs_class')['rank_metric'].apply(sum)

    srm_dict = jdf.groupby('feature_cs_class')['rank_metric'].apply(sum).to_dict()

    for k,val in srm_dict.items():
        colname = 'featCSclass_{}'.format(k)
        srm_dict[colname] = srm_dict.pop(k)

    return pd.Series(srm_dict)

def ToCorrectYrDate(x):
    # assume the input string is in the format of 'M/D/Y'
    # print x, type(x)
    if isinstance(x, basestring):
        check = re.search('(\d+)/(\d+)/(\d+)',x)
        if check:
            mm, dd, yy = [check.group(ii) for ii in range(1,4)]
            con_str = '/'
            if int(yy) <= 68:
                seq = (mm, dd, '19'+yy) #prefix the year with '19' instead of '20'
                return con_str.join(seq)
            else:
                return x
        else:
            print 'string not in the assumed date format M/D/Y!'
            return np.nan
    else:
        print 'not a string!'
        return np.nan

def To_Age_Dx(the_df):
    dob = the_df['DOB']
    date_dx = the_df['Date of Primary Dx']
    if isinstance(dob, basestring) and isinstance(date_dx, basestring):
        age_dx = ((pd.to_datetime(date_dx) - pd.to_datetime(dob)).days)/365.
        return age_dx
    else:
        return np.nan


def duration_func(x):
    #     print 'Recurrence = {}'.format(x['Recurrence'])
    #     print 'Dx: {}, Recur: {}, last contact: {}'.format(x['Date_Dx'], x['Date_Recur'], x['Date_lastcontactOrDeath'])
    if x['Recurrence'] == 1:
        if isinstance(x['Date_Dx'], pd.tslib.NaTType) is False and isinstance(x['Date_Recur'],
                                                                              pd.tslib.NaTType) is False:
            delta_time = (x['Date_Recur'] - x['Date_Dx']).days
            return delta_time
        else:
            return np.nan
    elif x['Recurrence'] == 0:
        if isinstance(x['Date_Dx'], pd.tslib.NaTType) is False and isinstance(x['Date_lastcontactOrDeath'],
                                                                              pd.tslib.NaTType) is False:
            delta_time = (x['Date_lastcontactOrDeath'] - x['Date_Dx']).days
            return delta_time
        else:
            return np.nan
    else:
        return np.nan


def relapse_free_observed(row):
    rf_dur = row['relapse_free_duration']
    if np.isnan(rf_dur) == False:
        if row['Recurrence'] == 1:
            return 1
        elif row['Recurrence'] == 0 and row['Alive/dead'] == 1:
            return 0
        elif row['Recurrence'] == 0 and row['Alive/dead'] == 0:
            return 0
        elif row['Recurrence'] == 2:
            return np.nan
        else:
            return np.nan
    else:
        return np.nan


def recurrence_free_observed(row):
    rf_dur = row['relapse_free_duration']
    if np.isnan(rf_dur) == False:
        if row['Recurrence'] == 1:
            return 1
        elif row['Recurrence'] == 0 and row['Alive/dead'] == 1:
            return 0
        elif row['Recurrence'] == 0 and row['Alive/dead'] == 0:
            return 1
        elif row['Recurrence'] == 2:
            return np.nan
        else:
            return np.nan
    else:
        return np.nan

def BCHistology_func(row):
    """
    :param x: a dataframe with tumor histology column name as Tumor_Histology and Tumor Histology
    :return: categroized histology described as below

    special case: Adenocarcinoma => look for more info her2_outcome data column titled 'Tumor Histology' for more info, if not, return Nan

    DCIS or LCIS: DCIS or LCIS
    IDC: Invasive ductal
    ILC: Invasive lobular
    Mixed IDC and ILC: Mixed invasive ductal and lobular

    """
    tumor_hist = row['Tumor_Histology'] # tumor histology from Bolouri's studies
    if tumor_hist is np.nan:
        return np.nan
    elif tumor_hist.lower() == 'adenocarcinoma':
        # look for more info
        tumor_hist_more = row['Tumor Histology']
        if tumor_hist_more is np.nan:
            return np.nan
        elif tumor_hist_more.lower() in ['adc', 'adenocarcinoma']:
            return np.nan
        elif re.search('invasive ductal (.+)', tumor_hist_more.lower()) or re.search('idc (.+)', tumor_hist_more.lower()):
            return 'IDC'
        elif re.search('invasive lobular (.+)', tumor_hist_more.lower()) or re.search('ilc (.+)', tumor_hist_more.lower()):
            return 'ILC'
    elif tumor_hist.lower() in ['invasive ductal', 'idc']:
        return 'IDC'
    elif tumor_hist.lower() in ['invasive lobular', 'ilc']:
        return 'ILC'
    elif tumor_hist.lower() == 'mixed invasive ductal and lobular':
        return 'Mixed IDC and ILC'
    elif tumor_hist.lower() in ['dcis', 'lcis']:
        return 'DCIS or LCIS'
