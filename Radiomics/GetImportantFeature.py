#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 12/28/17

@author: shuang
@goal: do one last training to deteremine feature importance

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

# read the learning output data and more
fname = '{}/CLF_output_all.csv'.format(data_dir)
df_clf_all = pd.read_csv(fname)

# drop the one with Random Forest since we would like to use logistic regression only
df_clf_all.drop(df_clf_all[(df_clf_all.Dep_var == 'DF_1yr') & (df_clf_all.clf_name == 'RandomForest')].index, inplace=True)
idx = df_clf_all.groupby('Dep_var').apply(lambda df: df.AUC_mean.argmax())
test = [idx[0]]

for ii in test:
    clf_name, oc, n_trial, k_fold = df_clf_all.ix[ii, ['clf_name', 'Dep_var', 'n_trial', 'k_fold']].tolist()
    json_fname = '{}/{}_IDV{}_DV{}_Trial{}_{}folds.json'.format(data_dir, clf_name, feat_name, oc, n_trial, k_fold)
    df_learning_output = pd.read_json(json_fname)
    print(df_learning_output.best_params[2])
