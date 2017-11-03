#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 10/30/17

@author: shuang
@goal: plot cluster plot with the features ranked by correlation coefficients

"""


import pandas as pd
import VizPlot as vp
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'

the_mri_tp = 2
the_pet_binwidth = 0.1
the_mri_binwidth = 5
im_dir = '{}/her2_Analysis/PETMRI/PETbinwidth{}_MRItp{}_binwidth{}'.format(rootdir,the_pet_binwidth, the_mri_tp, the_mri_binwidth)


fname = '{}/data_all.csv'.format(im_dir)
df_data_all = pd.read_csv(fname)

# look at the data for clustering algm settings
the_N_cluster = 3
# outcome_name_list = ['Recurrence_Type','Recurrence','Tumor_Grade','Tumor_Histology', 'T_stage', 'N_stage', 'Overall_stage','BC_subtype']
outcome_name_list = ['Recurrence','Tumor_Grade','BC_subtype']

# the_N_cluster = 2
# outcome_name_list = ['BoneMetsOrNot','TripleNeg'] + more_outcome

the_N_cluster = 4
outcome_name_list

chi2_outcome_fname = '{}/BestChi2_outcome_Ncluster{}.csv'.format(im_dir, the_N_cluster)
df_chi2_outcome_all = pd.read_csv(chi2_outcome_fname)

# get all the correlation coeff info
corr_tag = ['spearmancorr','corr']
for k in corr_tag:
    corr_fname = '{}/assoc_{}_all_v2_pvalcorr.csv'.format(im_dir,k)
    if k == 'spearmancorr':
        df_assoc_rs = pd.read_csv(corr_fname)
        rs_outcome_list = df_assoc_rs['outcome'].tolist()
    elif k == 'corr':
        df_assoc_rmrg = pd.read_csv(corr_fname)
        rmrg_outcome_list = df_assoc_rmrg['outcome'].tolist()

cluster_name = 'cs_class'
light_pal = sns.light_palette((210, 90, 60), input="husl")
the_color_pals = [light_pal[1:len(light_pal):2], sns.color_palette('Blues', len(df_data_all['T_stage'].unique())), sns.color_palette('YlOrRd', len(df_data_all['BC_subtype'].unique()))]
cc_prefix = 'patientConsensusCluster'
corr_thresh = 0.20
drop_feat_name = ['MRI_VOLUME_BLUE','MRI_VOLUME_GREEN','MRI_VOLUME_PURPLE','MRI_VOLUME_RED',
                  'MRI_VOLUME_TUMOR', 'MRI_VOLUME_TUM_BLU','MRI_VOLUME_WHITE','PET_SUV_max',
                  'MRI_PE2ROI_SER_MEAN','MRI_SERROI_PE2_PEAK']

for ov, pal in zip(outcome_name_list, the_color_pals):
    cm, cl, dm = df_chi2_outcome_all.ix[df_chi2_outcome_all['outcome_var'] == ov, ['cluster_method','cluster_linkage','dist_method']].values.flatten()
    print(ov)
    if isinstance(cl, str):
        cc_dir = '{}/{}_{}_{}_{}'.format(im_dir, cc_prefix, cm,cl,dm)
    else:
        cc_dir = '{}/{}_{}_{}'.format(im_dir, cc_prefix, cm, dm)

    # find the Ncluster that gives the highest median cluster consensus
    cc_fname = '{}/ClusterConsensus.csv'.format(cc_dir)
    cc_df = pd.read_csv(cc_fname)
    print(cc_df.groupby('k')['clusterConsensus'].median())

    # cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir, the_N_cluster)
    # cs_class_df = pd.read_csv(cs_class_fname)
    # cs_class_df.columns = ['ptid_side', 'cs_class']
    #
    # # combine the cs_class to the_df
    # the_df = pd.merge(df_data_all, cs_class_df, on='ptid_side')
    #
    # if ov in rs_outcome_list:
    #     df_corr = df_assoc_rs.ix[df_assoc_rs['outcome'] == ov, ['feature','corr_coeff']]
    # elif ov in rmrg_outcome_list:
    #     df_corr = df_assoc_rmrg.ix[df_assoc_rmrg['outcome'] == ov, ['feature', 'corr_coeff']]
    #
    # # drop the un-necessary features
    # df_corr = df_corr[~df_corr['feature'].isin(drop_feat_name)]
    # df_corr['abs_corr_coeff'] = df_corr['corr_coeff'].map(lambda x: np.abs(x))
    # df_corr.sort_values('abs_corr_coeff', ascending=False, inplace=True)
    #
    # # the_feat_order = df_corr.ix[df_corr['abs_corr_coeff'] >= corr_thresh, 'feature'].tolist()
    # the_feat_order = df_corr.feature.tolist()
    #
    # if ov == 'Recurrence':
    #     fig_fname = '{}/clustermap_Recur_kmax{}.pdf'.format(cc_dir, the_N_cluster)
    #     # fig_fname = '{}/clustermap_Recur_kmax{}_KeyFeat.pdf'.format(cc_dir, the_N_cluster)
    #     print(fig_fname)
    #     # recur_pal1 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3], sns.color_palette('Set2', 10)[7]]
    #     recur_pal1 = sns.color_palette('YlOrRd', len(the_df[ov].unique()))
    #     vp.ClustermapPlot(the_df, [ov], fig_fname, the_feat_order, vp.featurename_redef_petmr,
    #                       vp.featurelabel_redef_petmr, idvar_color_pal_mode=2, var_color_pal=[recur_pal1],
    #                       var_title=['Recurrence Status'])
    # elif ov == 'BC_subtype':
    #     fig_fname = '{}/clustermap_BCsubtype_kmax{}.pdf'.format(cc_dir, the_N_cluster)
    #     # fig_fname = '{}/clustermap_BCsubtype_kmax{}_KeyFeat.pdf'.format(cc_dir, the_N_cluster)
    #     print(fig_fname)
    #     vp.ClustermapPlot(the_df, [ov], fig_fname, the_feat_order, vp.featurename_redef_petmr,
    #                       vp.featurelabel_redef_petmr, idvar_color_pal_mode=2, var_title=['Breast cancer subtype'])
    #
    # elif ov == 'Tumor_Grade':
    #     fig_fname = '{}/clustermap_Grade_kmax{}.pdf'.format(cc_dir, the_N_cluster)
    #     # fig_fname = '{}/clustermap_Grade_kmax{}_KeyFeat.pdf'.format(cc_dir, the_N_cluster)
    #     print(fig_fname)
    #     light_pal = sns.light_palette((210, 90, 60), input="husl")
    #     grade_pal = light_pal[1:len(light_pal):2]
    #     vp.ClustermapPlot(the_df, [ov], fig_fname, the_feat_order, vp.featurename_redef_petmr,
    #                       vp.featurelabel_redef_petmr, idvar_color_pal_mode=2, var_title=['Tumor Grade'],
    #                       var_color_pal=[grade_pal])