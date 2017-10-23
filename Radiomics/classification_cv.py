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
import numpy as np

rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'

the_mri_tp = 2
the_mri_bin_width = 5
the_pet_bin_width = 0.1
im_dir = '{}/her2_Analysis/PETMRI/PETbinwidth{:.1f}_MRItp{}_binwidth{}'.format(rootdir,the_pet_bin_width, the_mri_tp, the_mri_bin_width)
fname = '{}/data_all.csv'.format(im_dir)
df_data = pd.read_csv(fname)
# df_data.rename(columns={'Sjoerd_Grade':'Tumor_Grade', 'Marjan_Histology': 'Tumor_Histology'},inplace=True)

# # define X and y
# # print(df_data.columns.tolist())
# pat = re.compile('texture_|FOstats_|ShapeSize_')
# feat_names = [ss for ss in df_data.columns.tolist() if pat.match(ss)]
# feat_tag = 'pet_mr_radiomics'
#
# # scale the features to z-score
# df_data[feat_names] = df_data[feat_names].apply(zscore)

feat_names = ['BC_subtype']
feat_tag = 'BC_subtype'

df_yr = range(1,6,1)
outcome_name_lst = ['DF_{}yr'.format(yy) for yy in df_yr]
# outcome_name_lst = ['TripleNeg']
# outcome_name_lst = ['DF_2yr']
# outcome_name_lst = ['Tumor_Grade']


# cross-validation setting
n_fold = 3
# n_trials = 100
n_trials = 10
coef_thresh = 0.5

# from sklearn import tree
# clf = tree.DecisionTreeClassifier()
# clf_name = 'DecTree'
# params_dict = {'criterion': ['gini','entropy'], 'splitter':['best','random'], 'max_depth': range(1,20,3), 'min_samples_split': range(2,6,1)}

# from sklearn.linear_model import LogisticRegression
# clf = LogisticRegression()
# clf_name = 'LogReg'
# clf_name = 'LogReg_{}'.format(feat_names[0])
# param for logistic regression
# params_dict = {'solver':['liblinear','sag','newton-cg','lbfgs','saga'], 'max_iter': [1000,5000,10000,20000], 'random_state': [random_state]}
# params_dict = {'solver': ['sag', 'newton-cg', 'lbfgs', 'saga'], 'multi_class': ['multinomial'], 'max_iter': [5000, 10000, 20000, 30000]}
# params_dict = {'solver':['liblinear','sag','newton-cg','lbfgs','saga'], 'max_iter': [5000,10000,20000,30000]}

# from sklearn.linear_model import SGDClassifier
# clf = SGDClassifier()
# clf_name = 'SGD'
# params_dict = {'loss':['log'],
#                'penalty':['l2','l1','elasticnet'],'l1_ratio': np.arange(0.0, 1.1,0.1),'max_iter':[5000, 100000,20000,30000]}

# from sklearn.svm import SVC
# clf = SVC()
# clf_name = 'SVM'
# params_dict = {'C': np.arange(0.00001,1.00001,0.1), 'kernel':['linear'], 'max_iter':[5000,10000,20000,30000],'probability':[True]}


from sklearn.ensemble import RandomForestClassifier
clf = RandomForestClassifier()
clf_name = 'RandomForest'
params_dict = {'n_estimators': range(1,15,2), 'criterion':['gini','entropy'],
               'max_features':['auto','sqrt','log2'], 'min_samples_split': range(2,12,2)}

for ocn in outcome_name_lst:
    df_tmp = df_data.ix[:,feat_names + [ocn]]
    df_tmp = df_tmp.dropna()
    X, y = df_tmp.ix[:,feat_names].as_matrix(), df_tmp.ix[:,ocn].astype('int').as_matrix()
    lt.nested_CV(X, y, clf, params_dict, n_fold, n_trials, feat_names, feat_tag, coef_thresh, im_dir, clf_name, ocn)