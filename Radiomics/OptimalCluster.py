#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 5/17/17

@author: shuang
@goal: determine the cluster setting that is optimal with the best cluster consensus

"""

import pandas as pd
import glob
import os
import numpy as np
import VizPlot as vp
import DataHelper as dh
import StatsTool as st
import seaborn as sns


rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'

# PET and MRI image setting
the_mri_tp = 2
the_mri_glcmbin = 128
the_pet_imgnorm_bin = 128
the_pet_glcmbin = 64

# get all the consensus clustering result and determine the BEST cluster result
im_dir = '{}/her2_Analysis/PETMRI/PETgbin{}_imgbin{}_MRItp{}_gbin{}'.format(rootdir,the_pet_glcmbin,the_pet_imgnorm_bin,the_mri_tp,the_mri_glcmbin)

df_medCC_all = dh.Median_Cluster_Consensus(im_dir, 'ConsensusCluster')

#display the maximum consensus clustering algorithm setting for Ncluster = 2,3,4, etc.
# idx = df_medCC_all.groupby('k')['medianCC'].max()

# get all the MRI and PET data in one space
# define the image feature of interest
texture_cols = ['texture_autocorrelation', 'texture_cluster_prominence', 'texture_cluster_shade',
                'texture_cluster_tendency', 'texture_contrast', 'texture_correlation', 'texture_diff_entropy',
                'texture_dissimilarity',
                'texture_energy', 'texture_entropy', 'texture_homogeneity1', 'texture_homogeneity2', 'texture_idmn',
                'texture_idn', 'texture_inv_var', 'texture_maxprob', 'texture_sum_avg', 'texture_sum_entropy',
                'texture_sum_var']

# re-format the radiomic features into columnes of feature name and feature values
img_feature_names = ['FOstats_energy', 'FOstats_entropy', 'FOstats_kurtosis', 'FOstats_mean', 'FOstats_min',
                     'FOstats_max', 'FOstats_skewness', 'FOstats_uniformity', 'FOstats_variance',
                     'ShapeSize_compactness1', 'ShapeSize_compactness2', 'ShapeSize_max_euc_dis',
                     'ShapeSize_spherical_disproportion', 'ShapeSize_sphericity', 'ShapeSize_surf_area_cm2',
                     'ShapeSize_surface2volratio', 'ShapeSize_vol_cm3'] + [ss + '_avg' for ss in texture_cols]

mri_fname = '{}/her2_Analysis/MRI/IsoVoxel_TP{}_GLCMBIN{}/MRIdataAll_tp{}_Nbin{}.csv'.format(rootdir,the_mri_tp,the_mri_glcmbin,the_mri_tp,the_mri_glcmbin)
df_mri = pd.read_csv(mri_fname,dtype={'PRIMARY_ID': str, 'MRN': str})
df_mri['BC_subtype'] = df_mri.apply(dh.BC_subtype_func, axis=1)

pet_fname = '{}/her2_Analysis/PET/IsoVoxel_IMGBIN{}_GLCMBIN{}/PETdataAll_glcmNbin{}_normNbin{}.csv'.format(rootdir,the_pet_imgnorm_bin,the_pet_glcmbin,the_pet_glcmbin,the_pet_imgnorm_bin)
df_pet = pd.read_csv(pet_fname, dtype={'pt_mrn': str, 'PRIMARY_ID': str, 'MRN': str, 'Anon_Accession_Num': str})
df_pet = df_pet.ix[:,img_feature_names + ['MRN','breast_side']]

# merge PET and MRI data to evaluate
df_petmr = pd.merge(df_mri, df_pet, left_on=['MRN','Laterality'], right_on=['MRN','breast_side'],suffixes=['_mri','_pet'])

# get other clinical data
outcome_fname = '{}/her2_ClinicalData/her2_outcome.csv'.format(rootdir)
df_her2_outcome = pd.read_csv(outcome_fname,dtype={'MRN': str,'Recurrence Type Summary': str})
df_her2_outcome['Recurrence_CR_1'] = df_her2_outcome['Recurrence Type Summary'].map(dh.recur_func1)
df_her2_outcome['Recurrence_CR_2'] = df_her2_outcome['Recurrence Type Summary'].map(dh.recur_func2)
df_her2_outcome['Diseasefree_5yr'] = df_her2_outcome.apply(dh.diseasefree5yr_func,axis=1)
df_her2_outcome['Recurrence_Type'] = df_her2_outcome['Recurrence Type Summary'].map(dh.recur_type_func)

# catergorize tumor and lymph node stages
df_her2_outcome['T_stage'] = df_her2_outcome['T-stage at surgery'].map(dh.Tstage_func)
df_her2_outcome['N_stage'] = df_her2_outcome['N-stage at surgery'].map(dh.Nstage_func)
df_her2_outcome['Overall_stage'] = df_her2_outcome['Overall stage'].map(dh.Overallstage_func)

# define the clinical outcome variables of interest
df_outcome = df_her2_outcome.loc[:, ['MRN','Diseasefree_5yr','Recurrence_CR_1', 'Recurrence_CR_2', 'Recurrence_Type','N_stage','T_stage','Overall_stage']]

jdf = pd.merge(df_petmr, df_outcome, on='MRN')

# save the joint data for further analysis for feature-to-outcome
fname = '{}/data_all.csv'.format(im_dir)
jdf.to_csv(fname)

# look at the data for clustering algm settings
outcome_name_list = ['Recurrence_Type', 'Diseasefree_5yr', 'Recurrence_CR_1', 'Recurrence_CR_2', 'TripleNeg',
                     'Sjoerd_Grade','Marjan_Histology', 'T_stage', 'N_stage', 'Overall_stage','BC_subtype']
cluster_name = 'cs_class'

chi2_all_df = pd.DataFrame(columns=('cluster_method','cluster_linkage','dist_method','N_mincluster','N_sample','v1','v2','chi2','pval','cramersV'))
count = 0
the_N_cluster = 3
df_tmp = df_medCC_all[df_medCC_all['k'] == the_N_cluster]

for ii in range(df_tmp.shape[0]):
    the_cm, the_cl, the_dm = df_tmp.ix[df_tmp.index[ii], ['cluster_method', 'cluster_linkage', 'dist_method']].values
    if isinstance(the_cl, str):
        cc_dir = '{}/ConsensusCluster_{}_{}_{}'.format(im_dir, the_cm, the_cl, the_dm)
    else:
        cc_dir = '{}/ConsensusCluster_{}_{}'.format(im_dir, the_cm, the_dm)

    cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir,the_N_cluster)
    cs_class_df = pd.read_csv(cs_class_fname)
    cs_class_df.columns = ['ptid_side', 'cs_class']

    # combine the cs_class to the_df
    the_df = pd.merge(jdf, cs_class_df, on='ptid_side')
    # print the_df.shape

    N_mincluster = the_df['cs_class'].value_counts().min()
    for v2 in outcome_name_list:
        # print 'v2 is {}'.format(v2)
        chi2,pval,dum1,dum2,N_sample,cv = st.chi2test(the_df,cluster_name,v2)
        chi2_all_df.loc[count] = [the_cm, the_cl, the_dm, N_mincluster, N_sample, cluster_name, v2, chi2, pval,cv]
        count = count + 1

# the_ov_interest = 'T_stage'
# print chi2_all_df[chi2_all_df['v2'] == the_ov_interest]

# get the optimal cluster setting for each outcome variable and plot the clustermap
fdf = chi2_all_df[chi2_all_df['N_mincluster'] > 5]
idx = fdf.groupby('v2').apply(lambda df: df.pval.argmin())

# determine the number of outcome variables
the_outcome_vars_list = chi2_all_df['v2'].unique().tolist()

imf_petmri_names = [ss + '_pet' for ss in img_feature_names] + [ss + '_mri' for ss in img_feature_names]
df_chi2_outcome_all = pd.DataFrame(columns=('outcome_var','cluster_method','cluster_linkage','dist_method','N_mincluster','N_sample','pval','cramersV','medianCC'))
count = 0
for ov in the_outcome_vars_list:
    ifdf = idx[ov]
    cm, cl, dm, N_mincluster, Nsample, pval, cv = fdf.ix[ifdf,['cluster_method','cluster_linkage','dist_method','N_mincluster','N_sample','pval','cramersV']].tolist()
    if isinstance(cl, str):
        cc_dir = '{}/ConsensusCluster_{}_{}_{}'.format(im_dir,cm,cl,dm)
    else:
        cc_dir = '{}/ConsensusCluster_{}_{}'.format(im_dir,cm,dm)

    cc_fname = '{}/ClusterConsensus.csv'.format(cc_dir)
    cc_df = pd.read_csv(cc_fname)
    medCC = cc_df.ix[cc_df['k'] == the_N_cluster,'clusterConsensus'].median()
    df_chi2_outcome_all.loc[count] = [ov, cm, cl, dm, N_mincluster, Nsample, pval, cv, medCC]
    count = count + 1

    if ov == 'BC_subtype':
        fig_fname = '{}/clustermap_BCsubtype_kmax{}.pdf'.format(cc_dir,the_N_cluster)
        vp.ClustermapPlot(the_df,['BC_subtype'],fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Breast cancer subtype'])
    elif ov == 'T_stage':
        fig_fname = '{}/clustermap_Tstage_kmax{}.pdf'.format(cc_dir,the_N_cluster)
        vp.ClustermapPlot(the_df,['T_stage'],fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['T-stage'])
    elif ov == 'N_stage':
        fig_fname = '{}/clustermap_Nstage_kmax{}.pdf'.format(cc_dir, the_N_cluster)
        vp.ClustermapPlot(the_df, ['N_stage'], fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['N-stage'])
    elif ov == 'Overall_stage':
        fig_fname = '{}/clustermap_Overallstage_kmax{}.pdf'.format(cc_dir, the_N_cluster)
        vp.ClustermapPlot(the_df, ['Overall_stage'], fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Overall stage'])
    elif ov == 'TripleNeg':
        fig_fname = '{}/clustermap_TripleNeg_kmax{}.pdf'.format(cc_dir, the_N_cluster)
        tripneg_pal = [(1., 1., 1.),sns.color_palette('Set2', 10)[2]]
        vp.ClustermapPlot(the_df, ['TripleNeg'], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Triple Neg'], var_color_pal=[tripneg_pal])
    elif ov == 'Sjoerd_Grade':
        fig_fname = '{}/clustermap_Grade_kmax{}.pdf'.format(cc_dir, the_N_cluster)
        light_pal = sns.light_palette((210, 90, 60), input="husl")
        grade_pal = light_pal[1:len(light_pal):2]
        vp.ClustermapPlot(the_df, ['Sjoerd_Grade'], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Grade'], var_color_pal=[grade_pal])
    elif ov == 'Marjan_Histology':
        fig_fname = '{}/clustermap_Histology_kmax{}.pdf'.format(cc_dir, the_N_cluster)
        colors = ['windows blue', 'amber', 'greyish', 'brick']
        histology_pal = sns.xkcd_palette(colors)
        vp.ClustermapPlot(the_df, ['Marjan_Histology'], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Histology'],var_color_pal=[histology_pal])
    elif ov == 'Recurrence_CR_1':
        fig_fname = '{}/clustermap_recurstatus1_kmax{}.pdf'.format(cc_dir, the_N_cluster)
        recur_pal1 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3],sns.color_palette('Set2', 10)[7]]
        vp.ClustermapPlot(the_df, ['Recurrence_CR_1'], fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Recurrence'],var_color_pal=[recur_pal1])
    elif ov == 'Recurrence_CR_2':
        fig_fname = '{}/clustermap_recurstatus2_kmax{}.pdf'.format(cc_dir, the_N_cluster)
        recur_pal2 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3]]
        vp.ClustermapPlot(the_df, ['Recurrence_CR_2'], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Recurrence'],var_color_pal=[recur_pal2])
    elif ov == 'Diseasefree_5yr':
        fig_fname = '{}/clustermap_df_5yr_kmax{}.pdf'.format(cc_dir, the_N_cluster)
        df5yr_pal = [(1., 1., 1.), sns.color_palette('Set2', 10)[3]]
        vp.ClustermapPlot(the_df, ['Diseasefree_5yr'], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['5-yr Disease Free'],var_color_pal=[df5yr_pal])
    elif ov == 'Recurrence_Type':
        fig_fname = '{}/clustermap_recurtype_kmax{}.pdf'.format(cc_dir, the_N_cluster)
        recurtype_pal = sns.light_palette('navy', n_colors=len(the_df['Recurrence_Type'].unique()), reverse=True)
        vp.ClustermapPlot(the_df, ['Recurrence_Type'], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Recurrence Type'],var_color_pal=[recurtype_pal])

print df_chi2_outcome_all


# # combine the cs_class to jdf (combined PET and MRI feature data)
# # determine the optimal clustering setting for a given N_cluster (determined by maximum median cluster consensus)
# idx = df_medCC_all.groupby('k').apply(lambda df: df.medianCC.argmax())
# df_tmp = df_medCC_all.loc[idx,:]
# the_N_cluster = 3
# the_cm, the_cl, the_dm, the_medCC = df_tmp.ix[df_tmp['k'] == the_N_cluster,['cluster_method','cluster_linkage','dist_method','medianCC']].values[0]
# # print the_cm, the_cl, the_dm,the_medCC
#
#
# # obtain the cluster class result
# the_cm = 'pam'
# the_dm = 'spearman'
# the_cl = np.nan
# if isinstance(the_cl, str):
#     cc_dir = '{}/ConsensusCluster_{}_{}_{}'.format(im_dir, the_cm, the_cl, the_dm)
# else:
#     cc_dir = '{}/ConsensusCluster_{}_{}'.format(im_dir, the_cm, the_dm)
#
# cc_file = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir,the_N_cluster)
# print 'read the file {}'.format(cc_file)
# df_cc = pd.read_csv(cc_file)
# df_cc.columns = ['ptid_side', 'cs_class']
# the_df = pd.merge(jdf, df_cc, on='ptid_side')

# # save dataframe for R testing
# fname = '{}/data_cs_all.csv'.format(cc_dir)
# the_df.to_csv(fname)

# # # Plot clustermap of radiomics vs. clinical outcome for samples
# imf_petmri_names = [ss + '_pet' for ss in img_feature_names] + [ss + '_mri' for ss in img_feature_names]
#
# # fig_fname = '{}/clustermap_BCsubtype_kmax{}.pdf'.format(cc_dir,the_N_cluster)
# # vp.ClustermapPlot(the_df,['BC_subtype'],fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Breast cancer subtype'])
#
# fig_fname = '{}/clustermap_Tstage_kmax{}.pdf'.format(cc_dir,the_N_cluster)
# vp.ClustermapPlot(the_df,['T_stage'],fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['T-stage'])

# fig_fname = '{}/clustermap_Nstage_kmax{}.pdf'.format(cc_dir, the_N_cluster)
# vp.ClustermapPlot(the_df, ['N_stage'], fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title='N-stage')
#
# fig_fname = '{}/clustermap_Overallstage_kmax{}.pdf'.format(cc_dir, the_N_cluster)
# vp.ClustermapPlot(the_df, ['Overall_stage'], fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title='Overall stage')
#
# fig_fname = '{}/clustermap_TripleNeg_kmax{}.pdf'.format(cc_dir, the_N_cluster)
# tripneg_pal = [(1., 1., 1.),sns.color_palette('Set2', 10)[2]]
# vp.ClustermapPlot(the_df, ['TripleNeg'], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Triple Neg'], var_color_pal=[tripneg_pal])
#
# fig_fname = '{}/clustermap_Grade_kmax{}.pdf'.format(cc_dir, the_N_cluster)
# light_pal = sns.light_palette((210, 90, 60), input="husl")
# grade_pal = light_pal[1:len(light_pal):2]
# vp.ClustermapPlot(the_df, ['Sjoerd_Grade'], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Grade'], var_color_pal=[grade_pal])
#
# fig_fname = '{}/clustermap_Histology_kmax{}.pdf'.format(cc_dir, the_N_cluster)
# colors = ['windows blue', 'amber', 'greyish', 'brick']
# histology_pal = sns.xkcd_palette(colors)
# vp.ClustermapPlot(the_df, ['Marjan_Histology'], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Histology'],var_color_pal=[histology_pal])
#
# fig_fname = '{}/clustermap_recurstatus1_kmax{}.pdf'.format(cc_dir, the_N_cluster)
# recur_pal1 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3],sns.color_palette('Set2', 10)[7]]
# vp.ClustermapPlot(the_df, ['Recurrence_CR_1'], fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Recurrence'],var_color_pal=[recur_pal1])
#
# fig_fname = '{}/clustermap_recurstatus2_kmax{}.pdf'.format(cc_dir, the_N_cluster)
# recur_pal2 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3]]
# vp.ClustermapPlot(the_df, ['Recurrence_CR_2'], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Recurrence'],var_color_pal=[recur_pal2])
#
# fig_fname = '{}/clustermap_df_5yr_kmax{}.pdf'.format(cc_dir, the_N_cluster)
# df5yr_pal = [(1., 1., 1.), sns.color_palette('Set2', 10)[3]]
# vp.ClustermapPlot(the_df, ['Diseasefree_5yr'], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['5-yr Disease Free'],
#                var_color_pal=[df5yr_pal])
#
# fig_fname = '{}/clustermap_recurtype_kmax{}.pdf'.format(cc_dir, the_N_cluster)
# recurtype_pal = sns.light_palette('navy', n_colors=len(the_df['Recurrence_Type'].unique()), reverse=True)
# vp.ClustermapPlot(the_df, ['Recurrence_Type'], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Recurrence Type'],var_color_pal=[recurtype_pal])
#
# fig_fname = '{}/clustermap_df5yr_TN_kmax{}.pdf'.format(cc_dir, the_N_cluster)
# varname = ['Diseasefree_5yr','TripleNeg']
# vartitle = ['5-yr Disease free','Triple Negative']
# pal1  = [(1., 1., 1.), sns.color_palette('Set2', 10)[3]]
# pal2 = [(1., 1., 1.), sns.color_palette('Set2', 10)[8]]
# var_pal_list = [pal1,pal2]
# vp.ClustermapPlot(the_df,varname, fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2,var_title=vartitle,var_color_pal=var_pal_list)


# # perform chi2 test for CC clusters vs clinical outcome
# chi2_all_df = pd.DataFrame(columns=('N_sample','v1','v2','chi2','pval'))
# count = 0
# for v2 in outcome_name_list:
#     # print 'v2 is {}'.format(v2)
#     chi2,pval,dum1,dum2,N_sample = st.chi2test(the_df,cluster_name,v2)
#     # print chi2,pval, N_sample
#     chi2_all_df.loc[count] = [N_sample,cluster_name,v2,chi2,pval]
#     count = count + 1
#
# print chi2_all_df

