#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 1/30/18

@author: shuang
@goal: obtain the learning output and confidence interval by boostrap approach

"""

import pandas as pd
import numpy as np
import LearningTool as lt


feat_name = 'pet_mr_radiomics'
# feat_name = 'Tumor_Grade'
# feat_name = 'BC_subtype'

df_yr = range(1,6,1)
outcome_names = ['DF_{}yr'.format(yy) for yy in df_yr]
# outcome_names = ['DF_2yr']

k_fold = 3

rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'
the_mri_tp = 2
the_mri_bin_width = 5
the_pet_bin_width = 0.1
data_dir = '{}/her2_Analysis/PETMRI/PETbinwidth{:.1f}_MRItp{}_binwidth{}/Learner'.format(rootdir,the_pet_bin_width, the_mri_tp, the_mri_bin_width)


# featimp_thresh_lut = {'LogReg': 0.6, 'SGD': 0.6, 'SVM': 0.6, 'RandomForest': 0.08}

clf_names = ['ElasticNet', 'L1LiblinearLogReg', 'L1SagaLogReg', 'L2LbfgsLogReg', 'L2SagLogReg', 'L2NewtoncgLogReg', 'L2LiblinearLogReg', 'RandomForest', 'SVM']
n_trials = [100] + (len(clf_names)-1)*[1000]

lst_data_all = []
data_col_names = ['clf_name', 'n_trial', 'k_fold', 'Indep_var', 'Dep_var', 'AUC_hat', 'CI_lo', 'CI_hi']

for oc in outcome_names:
    print('outcome: {}'.format(oc))

    # for clf_name, feat_thresh in featimp_thresh_lut.items():
    for clf_name, n_trial in zip(clf_names, n_trials):
        print('clf: {}'.format(clf_name))

        if clf_name == 'ElasticNet':
            the_df = lt.CombineFiles(data_dir, clf_name, oc, feat_name, n_trial, k_fold)
            # save the combined dataframe
            fname = '{}/{}_IDV{}_DV{}_Trial{}_{}folds.json'.format(data_dir, clf_name, feat_name, oc, len(the_df['trial_n'].unique()), k_fold)
            the_df.to_json(fname)
        else:
            fname = '{}/{}_IDV{}_DV{}_Trial{}_{}folds.json'.format(data_dir, clf_name, feat_name, oc, n_trial, k_fold)
            print(fname)
            the_df = pd.read_json(fname)

        AUC_star = the_df.groupby('trial_n').aggregate(np.mean)['auc'].as_matrix()
        AUC_hat = np.mean(AUC_star)
        delta = AUC_star - AUC_hat

        # determine 95% confidence interval (boostrap brute force approach)
        CI_lo = AUC_hat - np.percentile(delta, 97.5)
        CI_hi = AUC_hat - np.percentile(delta, 2.5)

        print('auc_mean: {}, 95% CI: {}, {}'.format(AUC_hat, CI_lo, CI_hi))

        tmp = dict(zip(data_col_names, [clf_name, len(delta), k_fold, feat_name, oc, AUC_hat, CI_lo, CI_hi]))
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
fname = '{}/CLF_output_all_boostrapCI.csv'.format(data_dir)
df_data_all.to_csv(fname, index=False)


