#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 9/14/17

@author: shuang
@goal: test out different classification algorithms

"""

import pandas as pd
import re
from scipy.stats import zscore
import LearningTool as lt
import sys
from sklearn.linear_model import SGDClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn import tree

if __name__ == '__main__':
    finput = sys.argv[1]
    if len(sys.argv) > 2:
        nRun = sys.argv[2]
    else:
        nRun = 0
    param_dict = lt.LoadInputParameter(finput)

    # load up nested CV settings
    n_trials = param_dict['n_trials']
    rootdir = param_dict['rootdir']
    n_fold = param_dict['n_fold']

    the_mri_tp = 2
    the_mri_bin_width = 5
    the_pet_bin_width = 0.1
    im_dir = '{}/her2_Analysis/PETMRI/PETbinwidth{:.1f}_MRItp{}_binwidth{}'.format(rootdir,the_pet_bin_width, the_mri_tp, the_mri_bin_width)
    fname = '{}/data_all.csv'.format(im_dir)
    df_data = pd.read_csv(fname)
    # df_data.rename(columns={'Sjoerd_Grade':'Tumor_Grade', 'Marjan_Histology': 'Tumor_Histology'},inplace=True)

    # define X and y
    # print(df_data.columns.tolist())
    pat = re.compile('texture_|FOstats_|ShapeSize_')
    feat_names = [ss for ss in df_data.columns.tolist() if pat.match(ss)]
    feat_tag = 'pet_mr_radiomics'

    # scale the features to z-score
    df_data[feat_names] = df_data[feat_names].apply(zscore)

    # feat_names = ['BC_subtype']
    # feat_tag = 'BC_subtype'

    # feat_names = ['Tumor_Grade']
    # feat_tag = 'Tumor_Grade'

    # df_yr = range(1,6,1)
    # outcome_name_lst = ['DF_{}yr'.format(yy) for yy in df_yr]
    # outcome_name_lst = ['TripleNeg']
    # outcome_name_lst = ['DF_2yr']
    # outcome_name_lst = ['Tumor_Grade']
    outcome_name_lst = param_dict['outcome_names']

    coef_thresh = 0.5

    # determine which classifier to use
    clf_name = param_dict['clf_name']
    clf_dict = {'ElasticNet': SGDClassifier, 'LogReg': LogisticRegression, 'SVM': SVC, 'RandomForest': RandomForestClassifier, 'DecTree': tree, 'extrDecTree': ExtraTreesClassifier}
    clf = clf_dict[clf_name]()
    hyper_param_dict = param_dict['hyper_params']

    for ocn in outcome_name_lst:
        df_tmp = df_data.ix[:,feat_names + [ocn]]
        df_tmp = df_tmp.dropna()
        X, y = df_tmp.ix[:,feat_names].as_matrix(), df_tmp.ix[:,ocn].astype('int').as_matrix()
        lt.nested_CV(X, y, clf, hyper_param_dict, n_fold, n_trials, feat_names, feat_tag, coef_thresh, im_dir, param_dict['clf_tag'], ocn, n_Run=nRun)
