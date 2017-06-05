#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 5/9/17

@author: shuang
@goal:

"""

import pandas as pd
import re
import itertools as itt
import numpy as np
import scipy.stats as stats


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


if __name__ == '__main__':
    # gather outcome data
    rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'
    outcome_fname = '{}/her2_ClinicalData/her2_outcome.csv'.format(rootdir)
    her2_outcome_df = pd.read_csv(outcome_fname,dtype={'MRN': str,'Recurrence Type Summary': str})


    # Recurrence_CR (recurrence info based on Cancer Registry field 'Recurrence Type Summary')
    # her2_outcome_df['Recurrence_CR_0'] = her2_outcome_df['Recurrence Type Summary'].map(recur_func0)
    her2_outcome_df['Recurrence_CR_1'] = her2_outcome_df['Recurrence Type Summary'].map(recur_func1)
    her2_outcome_df['Recurrence_CR_2'] = her2_outcome_df['Recurrence Type Summary'].map(recur_func2)
    her2_outcome_df['Diseasefree_5yr'] = her2_outcome_df.apply(diseasefree5yr_func,axis=1)
    her2_outcome_df['Recurrence_Type'] = her2_outcome_df['Recurrence Type Summary'].map(recur_type_func)

    # catergorize tumor and lymph node stages
    her2_outcome_df['T_stage'] = her2_outcome_df['T-stage at surgery'].map(Tstage_func)
    her2_outcome_df['N_stage'] = her2_outcome_df['N-stage at surgery'].map(Nstage_func)
    her2_outcome_df['Overall_stage'] = her2_outcome_df['Overall stage'].map(Overallstage_func)

    the_outcome_df = her2_outcome_df.loc[:, ['MRN','Diseasefree_5yr','Recurrence_CR_1', 'Recurrence_CR_2', 'Recurrence_Type','N_stage','T_stage','Overall_stage']]
    # the_outcome_df = her2_outcome_df.loc[:,['MRN'] + outcome_name_list]

    # PET image feature data
    pet_dir = '{}/her2_Analysis/PET/IsoVoxel_IMGBIN128_GLCMBIN64'.format(rootdir)
    df_fname = '{}/PETdataAll_glcmNbin64_normNbin128.csv'.format(pet_dir)
    data_df = pd.read_csv(df_fname, dtype={'pt_mrn': str, 'PRIMARY_ID': str, 'MRN': str, 'Anon_Accession_Num': str})
    data_df['BC_subtype'] = data_df.apply(BC_subtype_func,axis=1)
    jdf = pd.merge(data_df,the_outcome_df,on='MRN')
    # print data_df.shape, jdf.shape

    # dist_method_list = ['euclidean','spearman','pearson']
    # cluster_method_dict = {'pam': None, 'km': None, 'kmdist': None,'hc':['centroid','median','average','mcquitty','complete']}
    # theKmax_list = [2,3,4,5]
    # cluster_name = 'cs_class'
    # outcome_name_list = ['Recurrence_Type', 'Diseasefree_5yr', 'Recurrence_CR_1', 'Recurrence_CR_2', 'TripleNeg',
    #                      'Sjoerd_Grade', 'Marjan_Histology', 'T_stage', 'N_stage', 'Overall_stage', 'BC_subtype']

    #for testing
    dist_method_list = ['spearman']
    cluster_method_dict = {'hc': ['mcquitty']}
    theKmax_list = [3]
    cluster_name = 'cs_class'
    outcome_name_list = ['BC_subtype']

    for ncluster in theKmax_list:
        for k, val in cluster_method_dict.items():
            if val:
                for ll,dd in itt.product(val,dist_method_list):
                    cc_dir = '{}/ConsensusCluster_{}_{}_{}'.format(pet_dir,k,ll,dd)
                    cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir,ncluster)
                    cs_class_df = pd.read_csv(cs_class_fname)
                    cs_class_df.columns = ['ptid_side', 'cs_class']

                    # combine the cs_class to the_df
                    the_df = pd.merge(jdf, cs_class_df, how='left', on='ptid_side')
                    fname_save = '{}/dataall_csclass_kmax{}.csv'.format(cc_dir,ncluster)
                    the_df.to_csv(fname_save)
            else:
                for dd in dist_method_list:
                    cc_dir = '{}/ConsensusCluster_{}_{}'.format(pet_dir,k,dd)
                    cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir, ncluster)
                    cs_class_df = pd.read_csv(cs_class_fname)
                    cs_class_df.columns = ['ptid_side', 'cs_class']

                    # combine the cs_class to the_df
                    the_df = pd.merge(jdf, cs_class_df, how='left', on='ptid_side')
                    fname_save = '{}/dataall_csclass_kmax{}.csv'.format(cc_dir, ncluster)
                    the_df.to_csv(fname_save)

