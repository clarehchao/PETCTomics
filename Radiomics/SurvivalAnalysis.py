#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 8/8/17

@author: shuang
@goal: perform survival analysis with her2 dataset

"""

import pandas as pd
import numpy as np
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
import seaborn as sns
import DataHelper as dh
from lifelines.statistics import logrank_test

rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'
outcome_fname = '{}/her2_ClinicalData/her2_outcome_NJupdated.csv'.format(rootdir)
her2_outcome_df = pd.read_csv(outcome_fname,dtype={'MRN': str,'Recurrence Type Summary': str,'Brith Date': str, 'Date at Dx': str, 'Date Recurrence': str, 'Vital status': str})
her2_outcome_df = her2_outcome_df.replace('#VALUE!',np.nan)

her2_outcome_df['Date_Recur'] = pd.to_datetime(her2_outcome_df['Date Recurrence'])
her2_outcome_df['Date_Dx'] = pd.to_datetime(her2_outcome_df['Date of Primary Dx'])
her2_outcome_df['Date_lastcontactOrDeath'] = pd.to_datetime(her2_outcome_df['Last contact/death'])
her2_outcome_df['relapse_free_duration'] = her2_outcome_df.apply(dh.duration_func,axis=1)
# her2_outcome_df['relapse_observed'] = her2_outcome_df.apply(dh.relapse_free_observed,axis=1)
her2_outcome_df['recurrence_observed'] = her2_outcome_df.apply(dh.recurrence_free_observed,axis=1)

# combine all outcome data with radiomics data together
# PET and MRI image setting
the_mri_tp = 2
the_mri_glcmbin = 128
the_pet_imgnorm_bin = 128
the_pet_glcmbin = 64

# get all the consensus clustering result and determine the BEST cluster result
im_dir = '{}/her2_Analysis/PETMRI/PETgbin{}_imgbin{}_MRItp{}_gbin{}'.format(rootdir,the_pet_glcmbin,the_pet_imgnorm_bin,the_mri_tp,the_mri_glcmbin)
df_data_all = pd.read_csv('{}/data_all.csv'.format(im_dir),dtype={'MRN': str})
df_data_all = df_data_all.rename(columns={'Sjoerd_Grade': 'Tumor_Grade', 'Marjan_Histology': 'Tumor_Histology'})
# print df_data_all.columns.tolist()
outcome_list = ['Recurrence_Type', 'Diseasefree_5yr', 'BoneMetsOrNot','Recurrence_3', 'Recurrence_4','TripleNeg',
                'Tumor_Grade','Tumor_Histology', 'T_stage', 'N_stage', 'Overall_stage','BC_subtype']
df_outcome = df_data_all.ix[:,['MRN','ptid_side'] + outcome_list]

# get all cluster result
cc_prefix = 'patientConsensusCluster'
the_N_cluster = 3
df_medCC_all, df_oi = dh.Median_Cluster_Consensus(im_dir, cc_prefix)
df_tmp = df_medCC_all.ix[(df_medCC_all['k'] == the_N_cluster) & (df_medCC_all['N_mincluster'] > 5),:]

for ii in range(df_tmp.shape[0]):
    the_cm, the_cl, the_dm = df_tmp.ix[df_tmp.index[ii], ['cluster_method', 'cluster_linkage', 'dist_method']].values
    if isinstance(the_cl, str):
        cc_dir = '{}/{}_{}_{}_{}'.format(im_dir, cc_prefix, the_cm, the_cl, the_dm)
    else:
        cc_dir = '{}/{}_{}_{}'.format(im_dir, cc_prefix, the_cm, the_dm)

    cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir,the_N_cluster)
    cs_class_df = pd.read_csv(cs_class_fname)
    cs_class_df.columns = ['ptid_side', 'cs_class']

    # combine the cs_class to the_df
    jdf1 = pd.merge(df_outcome, cs_class_df, on='ptid_side')

    # join with her2 recurrence data
    jdf2 = pd.merge(jdf1,her2_outcome_df, on='MRN')

    kmf = KaplanMeierFitter()

    T = jdf2['relapse_free_duration'].dropna()
    E = jdf2['recurrence_observed'].dropna()
    kmf.fit(T, event_observed=E)
    kmf.plot()
    plt.xlabel('Recurrence-Free Survival (days)')
    if isinstance(the_cl, str):
        plt.title('{}, {}, {}'.format(the_cm, the_cl, the_dm))
    else:
        plt.title('{}, {}'.format(the_cm, the_dm))
    plt.show()

    ax = plt.subplot(111)
    the_subgroup_outcome = 'cs_class'
    subgroup_list = jdf2[the_subgroup_outcome].unique().tolist()
    for ii in subgroup_list:
        cc = (jdf2[the_subgroup_outcome] == ii)
        T2 = T[cc]
        E2 = E[cc]
        kmf.fit(T2, event_observed=E2, label='subgroup {}'.format(ii))
        kmf.plot(ax=ax)
        print 'class {}, median survival time: {}'.format(ii, kmf.median_)

    plt.xlabel('Recurrence-Free Survival (days)')
    plt.ylim(0, 1)
    if isinstance(the_cl, str):
        the_title = '{}, {}, {}'.format(the_cm, the_cl, the_dm)
    else:
        the_title = '{}, {}'.format(the_cm, the_dm)
    plt.title(the_title)
    plt.show()

    cc = (jdf2[the_subgroup_outcome] == 2)
    results = logrank_test(T[cc],T[~cc],E[cc],E[~cc])
    print the_title
    results.print_summary()