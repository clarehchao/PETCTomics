#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 9/19/17

@author: shuang
@goal: plot clustermap with multiple variables

"""

import pandas as pd
import VizPlot as vp
import DataHelper as dh
import StatsTool as st
import seaborn as sns
import re
import numpy as np


rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'

the_mri_tp = 2
the_pet_binwidth = 0.1
the_mri_binwidth = 5

im_dir = '{}/her2_Analysis/PETMRI/PETbinwidth{}_MRItp{}_binwidth{}'.format(rootdir,the_pet_binwidth, the_mri_tp, the_mri_binwidth)
cc_prefix = 'patientConsensusCluster'

fname = '{}/data_all.csv'.format(im_dir)
df_data = pd.read_csv(fname)

imf_pat = re.compile('texture_|ShapeSize_|FOstats_')
imf_petmri_names = [ss for ss in df_data.columns.tolist() if imf_pat.match(ss)] + ['MRI_SERROI_SER_MEAN','MRI_PE2ROI_PE2_PEAK']

cm = 'hc'
cl = 'complete'
dm = 'spearman'
the_N_cluster = 3

if isinstance(cl, str):
    cc_dir = '{}/{}_{}_{}_{}'.format(im_dir, cc_prefix, cm,cl,dm)
else:
    cc_dir = '{}/{}_{}_{}'.format(im_dir, cc_prefix, cm, dm)

cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir, the_N_cluster)
cs_class_df = pd.read_csv(cs_class_fname)
cs_class_df.columns = ['ptid_side', 'cs_class']

# combine the cs_class to the_df
the_df = pd.merge(df_data, cs_class_df, on='ptid_side')

the_outcome_of_interest = ['Tumor_Grade', 'T_stage', 'BC_subtype']

fig_fname = '{}/clustermap_Combo_kmax{}.pdf'.format(cc_dir, the_N_cluster)

light_pal = sns.light_palette((210, 90, 60), input="husl")
grade_pal = light_pal[1:len(light_pal):2]
tstage_pal = sns.color_palette('Blues', len(the_df['T_stage'].unique()))
bc_pal = sns.color_palette('YlOrRd', len(the_df['BC_subtype'].unique()))
vp.ClustermapPlot(the_df, the_outcome_of_interest, fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr, idvar_color_pal_mode=2,var_color_pal=[grade_pal, tstage_pal, bc_pal])