#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 6/19/18

@author: shuang
@goal: determine the final feature importance from nested CV result

"""

import LearningTool as lt

ID_var_name = 'pet_mr_radiomics'

n_trial = 1000
k_fold = 3

rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'
the_mri_tp = 2
the_mri_bin_width = 5
the_pet_bin_width = 0.1
im_dir = '{}/her2_Analysis/PETMRI/PETbinwidth{:.1f}_MRItp{}_binwidth{}'.format(rootdir,the_pet_bin_width, the_mri_tp, the_mri_bin_width)

# the selected model to use for all DF outcome
# df_yr = range(1,6,1)
# outcome_names = ['DF_{}yr'.format(yy) for yy in df_yr]
# the_clf_name = 'ElasticNet'
# lt.Feature_Importance(im_dir, the_clf_name, ID_var_name, n_trial, k_fold, outcome_names)

outcome_names = ['Tumor_Grade_Binary']
the_clf_name = 'L2lbfgsLogReg'
coef_thresh = 0.01
lt.Feature_Importance(im_dir, the_clf_name, ID_var_name, n_trial, k_fold, outcome_names, coef_thresh=coef_thresh)


