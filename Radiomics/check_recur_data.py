#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 7/5/17

@author: shuang
@goal: check the recurrence data consistency

"""

import pandas as pd
import numpy as np
import DataHelper as dh


rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'
outcome_fname = '{}/her2_ClinicalData/her2_outcome.csv'.format(rootdir)
her2_outcome_df = pd.read_csv(outcome_fname,dtype={'MRN': str,'Recurrence Type Summary': str})
her2_outcome_df['recur1'] = pd.to_numeric(her2_outcome_df['Recurrence Type Summary'].map(dh.recur_func1)).fillna(1000)
her2_outcome_df['Recurrence'] = pd.to_numeric(her2_outcome_df['Recurrence'].map(lambda x: np.nan if x == '#VALUE!' else x)).fillna(1000)
her2_outcome_df['recur_diff'] = abs(her2_outcome_df['recur1'] - her2_outcome_df['Recurrence'])
oh = her2_outcome_df.ix[(her2_outcome_df['recur_diff'] == 1.0) | (her2_outcome_df['recur_diff'] == 2.0), :]
oh.to_csv('{}/her2_ClinicalData/Recur_diff_MRN_2.csv'.format(rootdir), index=False, index_label=False)