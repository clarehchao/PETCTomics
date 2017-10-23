#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 10/10/17

@author: shuang
@goal: extract info from the nested cross-validation for reports

"""

import pandas as pd
import numpy as np


clf_name = 'LogReg'
df_yr = range(1,6,1)
# feat_name = 'pet_mr_radiomics'
# feat_name = 'Tumor_Grade'
feat_name = 'BC_subtype'
outcome_names = ['DF_{}yr'.format(yy) for yy in df_yr]
n_trial = 10
k_fold = 3

rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'
the_mri_tp = 2
the_mri_bin_width = 5
the_pet_bin_width = 0.1
data_dir = '{}/her2_Analysis/PETMRI/PETbinwidth{:.1f}_MRItp{}_binwidth{}/Learner'.format(rootdir,the_pet_bin_width, the_mri_tp, the_mri_bin_width)

feat_importance_thresh = 0.6


for oc in outcome_names:
    print('outcome: {}'.format(oc))
    fname = '{}/{}_IDV{}_DV{}_Trial{}_{}folds.json'.format(data_dir, clf_name, feat_name, oc, n_trial, k_fold)
    print(fname)
    # get mean and stdev AUC among all the trials
    the_df = pd.read_json(fname)
    auc_mean = the_df.groupby('trial_n').aggregate(np.mean)['auc'].mean()
    auc_std = the_df.groupby('trial_n').aggregate(np.mean)['auc'].std()
    print('mean AUC: {}, stdev: {}'.format(auc_mean, auc_std))

    # get a sense of key features
    key_feat_list = []
    for nn in range(n_trial):
        # print('nn = {}'.format(nn))
        np_feat_importance = np.abs(np.vstack(the_df.ix[the_df['trial_n'] == nn, 'feat_importance'].values))
        np_feat_name = np.hstack(the_df.ix[(the_df['trial_n'] == nn) & (the_df['fold_n'] == 0), 'feat_name'])
        avg_keyfeat = np.mean(np_feat_importance, axis=0)

        tmp = np_feat_name[avg_keyfeat > feat_importance_thresh].tolist()
        key_feat_list = tmp + key_feat_list
        # print(avg_keyfeat)
        # print(np.std(np_feat_importance, axis=0))

    print('key features: {}'.format(list(set(key_feat_list))))




