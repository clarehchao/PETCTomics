#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 5/31/17

@author: shuang
@goal: determine the best cluster setting for feature cluster and plot heatmaps of the association test p-value and corr coefficient

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
# assoc_fname = '{}/assoc_corr_test_result_all_wrt_pt_cs.csv'.format(im_dir)
# assoc_fname = '{}/assoc_spearmancorr_all.csv'.format(im_dir)
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

# define the row label
df_cs.index = df_cs['feature']
row_labels = df_cs.ix[idx,'cs_class'].tolist()
col_label_lut = {'BC_subtype':'Breast cancer subtype','T_stage':'T-stage', 'N_stage': 'N-stage',
                 'Overall_stage': 'Overall-stage', 'TripleNeg': 'Triple Neg', 'Sjoerd_Grade': 'Tumor Grade',
                 'Marjan_Histology': 'Tumor Histology', 'Recurrence_CR_1':'Recurrence_3grp', 'Recurrence_CR_2':'Recurrence_2grp',
                 'Diseasefree_5yr': '5-yr Disease Free', 'Recurrence_Type': 'Recurrence site'}


# col_label_lut = {'BC_subtype':'Breast cancer subtype','T_stage':'T-stage', 'N_stage': 'N-stage',
#                  'Overall_stage': 'Overall-stage', 'TripleNeg': 'Triple Neg', 'Tumor_Grade': 'Tumor Grade',
#                  'Tumor_Histology': 'Tumor Histology', 'Recurrence_CR_1':'Recurrence_3grp', 'Recurrence_CR_2':'Recurrence_2grp',
#                  'Diseasefree_5yr': '5-yr Disease Free', 'Recurrence_Type': 'Recurrence site'}

# imf_order = ['Tumor Grade', 'Breast cancer subtype', 'Tumor Histology', 'T-stage', 'N-stage', 'Overall-stage', 'Recurrence_3grp', 'Recurrence site']
# imf_order = df_assoc['outcome'].unique().tolist()
# imf_order = ['Tumor Grade','T-stage','N-stage','Overall-stage']
imf_order = ['Tumor Grade','Breast cancer subtype','T-stage','N-stage','Overall-stage']
imf_order.reverse()

# heatmap of KW H test p-val
the_tab1 = df_assoc.pivot('feature', 'outcome', 'pval')
the_tab1 = the_tab1.ix[idx]
the_tab1.columns = [col_label_lut[ss] for ss in the_tab1.columns]
the_tab1_oi = the_tab1.ix[:,imf_order]

# the mask: true, will mask data, false: will NOT mask data
row_label_title = 'feature CS class'
fig_name = '{}/featureVSoutcome_pval_Ncluster{}_nomask.pdf'.format(cc_dir, the_Ncluster)
# fig_name = '{}/featureVSoutcome_pval_Ncluster{}_nomask_spc.pdf'.format(cc_dir, the_Ncluster)
vp.clustermap_plot_simple(the_tab1_oi, row_labels, row_label_title=row_label_title, fig_name=fig_name, value_title='P-value')

fig_name = '{}/featureVSoutcome_pval_Ncluster{}_withmask.pdf'.format(cc_dir, the_Ncluster)
# fig_name = '{}/featureVSoutcome_pval_Ncluster{}_withmask_spc.pdf'.format(cc_dir, the_Ncluster)
the_mask = the_tab1_oi > 0.05
vminmax = [0, 0.05]
vp.clustermap_plot_simple(the_tab1_oi, row_labels, row_label_title=row_label_title, vminmax=vminmax,mask=the_mask, fig_name=fig_name, value_title='P-value')

# heatmap of corr coeff
the_tab2 = df_assoc.pivot('feature', 'outcome', 'corr_coeff')
the_tab2 = the_tab2.ix[idx]
the_tab2.columns = [col_label_lut[ss] for ss in the_tab2.columns]
the_tab2_oi = the_tab2.ix[:,imf_order]

fig_name = '{}/featureVSoutcome_corrcoeff_Ncluster{}_nomask.pdf'.format(cc_dir, the_Ncluster)
# fig_name = '{}/featureVSoutcome_spearmancorr_Ncluster{}_nomask.pdf'.format(cc_dir, the_Ncluster)
vp.clustermap_plot_simple(the_tab2_oi, row_labels, row_label_title=row_label_title, fig_name=fig_name, value_title='Multiple\nCorr. Coeff.')

# the mask: true, will mask data, false: will NOT mask data
fig_name = '{}/featureVSoutcome_corrcoeff_Ncluster{}_withmask.pdf'.format(cc_dir, the_Ncluster)
# fig_name = '{}/featureVSoutcome_spearmancorr_Ncluster{}_withmask.pdf'.format(cc_dir, the_Ncluster)
# the_mask = the_tab2_oi < 0.3
# the_mask = the_tab1_oi > 0.05
cor_coeff_cutoff = 0.3
the_mask = (the_tab1_oi > 0.05) | (abs(the_tab2_oi) < cor_coeff_cutoff)
# vminmax = [df_assoc['corr_coeff'].min(), df_assoc['corr_coeff'].max()]
vminmax = [cor_coeff_cutoff, df_assoc['corr_coeff'].max()]
vp.clustermap_plot_simple(the_tab2_oi, row_labels, row_label_title=row_label_title, vminmax=vminmax, mask=the_mask, fig_name=fig_name, value_title='Multiple\nCorr. Coeff.')




