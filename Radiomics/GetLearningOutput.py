#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 10/10/17

@author: shuang
@goal: extract info from the nested cross-validation for reports

"""

import pandas as pd
import numpy as np


feat_name = 'pet_mr_radiomics'
# feat_name = 'Tumor_Grade'
# feat_name = 'BC_subtype'

df_yr = range(1,6,1)
outcome_names = ['DF_{}yr'.format(yy) for yy in df_yr]
# outcome_names = ['DF_2yr']

n_trial = 10
k_fold = 3

rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'
the_mri_tp = 2
the_mri_bin_width = 5
the_pet_bin_width = 0.1
data_dir = '{}/her2_Analysis/PETMRI/PETbinwidth{:.1f}_MRItp{}_binwidth{}/Learner'.format(rootdir,the_pet_bin_width, the_mri_tp, the_mri_bin_width)


# featimp_thresh_lut = {'LogReg': 0.6, 'SGD': 0.6, 'SVM': 0.6, 'RandomForest': 0.08}

clf_names = ['ElasticNet', 'L1LiblinearLogReg', 'L1SagaLogReg', 'L2LbfgsLogReg', 'L2SagLogReg', 'L2NewtoncgLogReg', 'L2LiblinearLogReg', 'RandomForest', 'SVM']

lst_data_all = []
data_col_names = ['clf_name', 'n_trial', 'k_fold', 'Indep_var', 'Dep_var', 'AUC_mean', 'CI_lo_mean', 'CI_hi_mean']

for oc in outcome_names:
    print('outcome: {}'.format(oc))

    # for clf_name, feat_thresh in featimp_thresh_lut.items():
    for clf_name in clf_names:
        print('clf: {}'.format(clf_name))
        fname = '{}/{}_IDV{}_DV{}_Trial{}_{}folds.json'.format(data_dir, clf_name, feat_name, oc, n_trial, k_fold)

        # get mean and stdev AUC among all the trials
        the_df = pd.read_json(fname)
        auc_mean = the_df.groupby('trial_n').aggregate(np.mean)['auc']
        auc_std = the_df.groupby('trial_n').aggregate(np.std)['auc']
        if clf_name in ['RandomForest', 'SVM']:
            Nsample_mean = expl_Nsamplemean
        else:
            Nsample_mean = the_df.groupby('trial_n').aggregate(np.mean)['N_sample']

        if clf_name == 'ElasticNet':
            expl_Nsamplemean = Nsample_mean
        CI_lo = auc_mean - 1.96*auc_std/np.sqrt(Nsample_mean)
        CI_hi = auc_mean + 1.96*auc_std/np.sqrt(Nsample_mean)
        print('auc_mean: {}, 95% CI: {}, {}'.format(auc_mean.mean(), CI_lo.mean(), CI_hi.mean()))

        tmp = dict(zip(data_col_names, [clf_name, n_trial, k_fold, feat_name, oc, auc_mean.mean(), CI_lo.mean(), CI_hi.mean()]))
        lst_data_all.append(tmp)

        # auc_mean = the_df.groupby('trial_n').aggregate(np.mean)['auc'].mean()
        # auc_std = the_df.groupby('trial_n').aggregate(np.mean)['auc'].std()
        # print('{}, mean AUC: {}, stdev: {}'.format(clf_name, auc_mean, auc_std))

        # # get a sense of key features
        # key_feat_list = []
        # for nn in range(n_trial):
        #     # print('nn = {}'.format(nn))
        #     np_feat_importance = np.abs(np.vstack(the_df.ix[the_df['trial_n'] == nn, 'feat_importance'].values))
        #     # print(np_feat_importance)
        #     np_feat_name = np.hstack(the_df.ix[(the_df['trial_n'] == nn) & (the_df['fold_n'] == 0), 'feat_name'])
        #     avg_keyfeat = np.mean(np_feat_importance, axis=0)
        #
        #     tmp = np_feat_name[avg_keyfeat > feat_thresh].tolist()
        #     key_feat_list = tmp + key_feat_list
        #     # print(avg_keyfeat)
        #     # print(np.std(np_feat_importance, axis=0))
        #
        # print('{}, key features: {}'.format(clf_name, list(set(key_feat_list))))

df_data_all = pd.DataFrame(lst_data_all)
df_data_all = df_data_all.ix[:, data_col_names]
fname = '{}/CLF_output_all.csv'.format(data_dir)
df_data_all.to_csv(fname, index=False)


