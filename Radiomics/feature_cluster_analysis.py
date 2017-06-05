#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 5/31/17

@author: shuang
@goal: determine the best cluster setting for feature cluster

"""

import pandas as pd
import DataHelper as dh
import VizPlot as vp

rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'

# PET and MRI image setting
the_mri_tp = 2
the_mri_glcmbin = 128
the_pet_imgnorm_bin = 128
the_pet_glcmbin = 64

# get all the consensus clustering result and determine the BEST cluster result
im_dir = '{}/her2_Analysis/PETMRI/PETgbin{}_imgbin{}_MRItp{}_gbin{}'.format(rootdir,the_pet_glcmbin,the_pet_imgnorm_bin,the_mri_tp,the_mri_glcmbin)

df_medCC_all, df_oi = dh.Median_Cluster_Consensus(im_dir, 'featureConsensusCluster')

# plot the heatmap of assoc p-value and corr coeff
assoc_fname = '{}/assoc_corr_test_result_all.csv'.format(im_dir)
df_assoc = pd.read_csv(assoc_fname)

the_Ncluster = 3
the_cm, the_cl, the_dm = df_oi.ix[df_oi['k'] == the_Ncluster,['cluster_method','cluster_linkage','dist_method']].as_matrix()[0]

# read the feature cluster result
if isinstance(the_cl, str):
    cc_dir = '{}/featureConsensusCluster_{}_{}_{}'.format(im_dir, the_cm, the_cl, the_dm)
else:
    cc_dir = '{}/featureConsensusCluster_{}_{}'.format(im_dir, the_cm, the_dm)
cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir, the_Ncluster)
df_cs = pd.read_csv(cs_class_fname)
df_cs.columns = ['feature', 'cs_class']


df_cs = df_cs.sort_values('cs_class')
idx = df_cs['feature'].tolist()
the_tab = df_assoc.pivot('feature', 'outcome', 'pval')
the_tab = the_tab.ix[idx]
fig_name = '{}/featureVSoutcome_pval.pdf'.format(cc_dir)

# define the row label
df_cs.index = df_cs['feature']
row_labels = df_cs.ix[idx,'cs_class'].tolist()
print row_labels
vp.heatmap_plot(the_tab, row_labels, fig_name = fig_name)

