#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 5/17/17

@author: shuang
@goal: determine the cluster setting that is optimal with the best cluster consensus

"""

import pandas as pd
import VizPlot as vp
import DataHelper as dh
import StatsTool as st
import seaborn as sns
import re
import numpy as np


rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'

# PET and MRI image setting
# the_mri_tp = 2
# the_mri_glcmbin = 128
# the_pet_imgnorm_bin = 128
# the_pet_glcmbin = 128

the_mri_tp = 2
the_pet_binwidth = 0.1
the_mri_binwidth = 5

# get all the consensus clustering result and determine the BEST cluster result
# im_dir = '{}/her2_Analysis/PETMRI/PETgbin{}_imgbin{}_MRItp{}_gbin{}'.format(rootdir,the_pet_glcmbin,the_pet_imgnorm_bin,the_mri_tp,the_mri_glcmbin)
im_dir = '{}/her2_Analysis/PETMRI/PETbinwidth{}_MRItp{}_binwidth{}'.format(rootdir,the_pet_binwidth, the_mri_tp, the_mri_binwidth)
Nrep = 10000
# cc_prefix = 'patientConsensusCluster'
cc_prefix = 'patientConsensusCluster_Rep{}'.format(Nrep)

#display the maximum consensus clustering algorithm setting for Ncluster = 2,3,4, etc.
# idx = df_medCC_all.groupby('k')['medianCC'].max()

# get all the MRI and PET data in one space
# define the image feature of interest
# texture_cols = ['texture_autocorrelation', 'texture_cluster_prominence', 'texture_cluster_shade',
#                 'texture_cluster_tendency', 'texture_contrast', 'texture_correlation', 'texture_diff_entropy',
#                 'texture_dissimilarity',
#                 'texture_energy', 'texture_entropy', 'texture_homogeneity1', 'texture_homogeneity2', 'texture_idmn',
#                 'texture_idn', 'texture_inv_var', 'texture_maxprob', 'texture_sum_avg', 'texture_sum_entropy',
#                 'texture_sum_var']
# # re-format the radiomic features into columnes of feature name and feature values
# img_feature_names = ['FOstats_energy', 'FOstats_entropy', 'FOstats_kurtosis', 'FOstats_mean', 'FOstats_min',
#                      'FOstats_max', 'FOstats_skewness', 'FOstats_uniformity', 'FOstats_variance',
#                      'ShapeSize_compactness1', 'ShapeSize_compactness2', 'ShapeSize_max_euc_dis',
#                      'ShapeSize_spherical_disproportion', 'ShapeSize_sphericity', 'ShapeSize_surf_area_cm2',
#                      'ShapeSize_surface2volratio', 'ShapeSize_vol_cm3'] + [ss + '_avg' for ss in texture_cols]


feature_list = ['autocorrelation', 'cluster_prominence', 'cluster_shade', 'cluster_tendency', 'contrast',
                    'correlation','diff_entropy', 'dissimilarity', 'energy', 'entropy', 'homogeneity1', 'homogeneity2',
                    'idmn','idn', 'inv_var', 'maxprob', 'sum_avg', 'sum_entropy', 'sum_var', 'imc1', 'imc2', 'diff_avg',
                    'diff_var', 'avg_intensity', 'sum_squares']

# re-format the radiomic features into columnes of feature name and feature values
img_feature_names = ['FOstats_energy', 'FOstats_entropy', 'FOstats_kurtosis', 'FOstats_mean', 'FOstats_min',
                     'FOstats_max', 'FOstats_skewness', 'FOstats_uniformity', 'FOstats_variance',
                     'ShapeSize_compactness1', 'ShapeSize_compactness2', 'ShapeSize_max_euc_dis',
                     'ShapeSize_spherical_disproportion', 'ShapeSize_sphericity', 'ShapeSize_surf_area_cm2',
                     'ShapeSize_surface2volratio', 'ShapeSize_vol_cm3'] + ['texture_{}_avg'.format(ss) for ss in feature_list]

# mri_fname = '{}/her2_Analysis/MRI/IsoVoxel_TP{}_GLCMBIN{}/MRIdataAll_tp{}_Nbin{}.csv'.format(rootdir,the_mri_tp,the_mri_glcmbin,the_mri_tp,the_mri_glcmbin)
mri_fname = '{}/her2_Analysis/MRI/IsoVoxel_binWidth{}/MRIdataAll_tp{}_binWidth{}.csv'.format(rootdir,the_mri_binwidth,the_mri_tp,the_mri_binwidth)
df_mri = pd.read_csv(mri_fname,dtype={'PRIMARY_ID': str, 'MRN': str})
df_mri['BC_subtype'] = df_mri.apply(dh.BC_subtype_func, axis=1)
df_mri_colnames = df_mri.columns.tolist()
df_mri_colnames[df_mri_colnames.index('SUV max')] = 'PET_SUV_max'
df_mri.columns = df_mri_colnames

# pet_fname = '{}/her2_Analysis/PET/IsoVoxel_IMGBIN{}_GLCMBIN{}/PETdataAll_glcmNbin{}_normNbin{}.csv'.format(rootdir,the_pet_imgnorm_bin,the_pet_glcmbin,the_pet_glcmbin,the_pet_imgnorm_bin)
pet_fname = '{}/her2_Analysis/PET/IsoVoxel_binWidth{:.1f}/PETdataAll_binWidth{:.1f}.csv'.format(rootdir,the_pet_binwidth, the_pet_binwidth)
df_pet = pd.read_csv(pet_fname, dtype={'pt_mrn': str, 'PRIMARY_ID': str, 'MRN': str, 'Anon_Accession_Num': str})
df_pet = df_pet.ix[:,img_feature_names + ['MRN','breast_side']]

# merge PET and MRI data to evaluate
df_petmr = pd.merge(df_mri, df_pet, left_on=['MRN','Laterality'], right_on=['MRN','breast_side'],suffixes=['_mri','_pet'])
#TODO: why was the shape of df_petmr different from the shape of tmp3 (see her2_data_img_org.py in /data/francgrp1/...)

# get other clinical data
# outcome_fname = '{}/her2_ClinicalData/her2_outcome.csv'.format(rootdir)
outcome_fname = '{}/her2_ClinicalData/her2_outcome_NJupdated.csv'.format(rootdir)
df_her2_outcome = pd.read_csv(outcome_fname,dtype={'MRN': str,'Recurrence Type Summary': str,'Brith Date': str, 'Date at Dx': str})

# MAKE SURE MRN are in the same format as MRN from df_pet and df_mri (saved and written by ClusterAnalysis_PET.py and ClusterAnalysis_MRI.py on petct:/data/francgrp1/clare_work/PETCTomics)
# or else, some MRN were missed!
df_her2_outcome['MRN'] = df_her2_outcome['MRN'].apply(lambda x: '{0:0>8}'.format(x))
df_her2_outcome = df_her2_outcome.replace('#VALUE!',np.nan)
df_her2_outcome['DOB'] = df_her2_outcome['Brith Date'].map(dh.ToCorrectYrDate)
df_her2_outcome['Age_Dx'] = df_her2_outcome.apply(dh.To_Age_Dx, axis=1)
print(df_her2_outcome['Age_Dx'].median(), df_her2_outcome['Age_Dx'].min(), df_her2_outcome['Age_Dx'].max())

# df_her2_outcome['Recurrence_3'] = df_her2_outcome.apply(dh.recur_func3, axis=1)
# df_her2_outcome['Recurrence_4'] = df_her2_outcome['Recurrence_3'].map(lambda x: np.nan if x == 2 else x)
# df_her2_outcome['Recurrence_CR_1'] = df_her2_outcome['Recurrence Type Summary'].map(dh.recur_func1)
# df_her2_outcome['Recurrence_CR_2'] = df_her2_outcome['Recurrence Type Summary'].map(dh.recur_func2)
# df_her2_outcome['Diseasefree_5yr'] = df_her2_outcome.apply(dh.diseasefree5yr_func,axis=1)
df_her2_outcome['Recurrence_Type'] = df_her2_outcome['Recurrence Type Summary'].map(dh.recur_type_func)
df_her2_outcome['BoneMetsOrNot'] = df_her2_outcome['Recurrence Type Summary'].map(dh.BoneRecurOrNot_func)

# catergorize tumor and lymph node stages
df_her2_outcome['T_stage'] = df_her2_outcome['T-stage at surgery'].map(dh.Tstage_func)
df_her2_outcome['N_stage'] = df_her2_outcome['N-stage at surgery'].map(dh.Nstage_func)
df_her2_outcome['Overall_stage'] = df_her2_outcome['Overall stage'].map(dh.Overallstage_func)
df_her2_outcome['Overall_stage_2'] = df_her2_outcome['Overall stage'].map(dh.Overallstage_func2)

# recurrence free and survival free duration info
df_her2_outcome['Date_Recur'] = pd.to_datetime(df_her2_outcome['Date Recurrence'])
df_her2_outcome['Date_Dx'] = pd.to_datetime(df_her2_outcome['Date of Primary Dx'])
df_her2_outcome['Date_lastcontactOrDeath'] = pd.to_datetime(df_her2_outcome['Last contact/death'])

n_yr_list = [1., 2., 3., 4., 5.]
more_outcome = []
for yy in n_yr_list:
    DF_colname = 'DF_{:.0f}yr'.format(yy)
    OS_colname = 'OS_{:.0f}yr'.format(yy)
    more_outcome.extend([DF_colname, OS_colname]) # don't use append when adding with a list
    df_her2_outcome[DF_colname] = df_her2_outcome.apply(dh.Diseasefree_byYear_func, n_year = yy, axis=1)
    df_her2_outcome[OS_colname] = df_her2_outcome.apply(dh.OverallSurvival_byYear_func, n_year=yy, axis=1)


# define the clinical outcome variables of interest
# the_clinical_vars = ['MRN','Diseasefree_5yr','BoneMetsOrNot','Recurrence','Recurrence_Type','N_stage','T_stage', 'Tumor Histology','Overall_stage','Overall_stage_2']
the_clinical_vars = ['MRN','BoneMetsOrNot','Recurrence','Recurrence_Type','N_stage','T_stage', 'Tumor Histology','Overall_stage','Overall_stage_2'] + more_outcome

df_outcome = df_her2_outcome.loc[:, the_clinical_vars]

jdf = pd.merge(df_petmr, df_outcome, on='MRN')

# save the joint data for further analysis for feature-to-outcome

# change clinical outcome name
jdf = jdf.rename(columns={'Sjoerd_Grade': 'Tumor_Grade', 'Marjan_Histology': 'Tumor_Histology'})
jdf['Tumor_Histology_updated'] = jdf.apply(dh.BCHistology_func, axis=1)

print(jdf.columns.tolist())

# make tumor grade into binary category: 0: G1, G2, 1: G3, G4
tg_mapping = {1:0, 2:0, 3:1, 4:1}
jdf['Tumor_Grade_Binary'] = jdf['Tumor_Grade'].map(tg_mapping)

fname = '{}/data_all.csv'.format(im_dir)
jdf.to_csv(fname)

# look at the data for clustering algm settings
# the_N_cluster = 3
the_N_cluster = [3]
# outcome_name_list = ['Recurrence_Type','Recurrence','Tumor_Grade','Tumor_Histology', 'T_stage', 'N_stage', 'Overall_stage','BC_subtype']
outcome_name_list = ['Recurrence_Type','Recurrence','Tumor_Grade','Tumor_Histology_updated',
                     'T_stage', 'N_stage', 'Overall_stage','BC_subtype','Overall_stage_2']

# the_N_cluster = 2
# outcome_name_list = ['BoneMetsOrNot','TripleNeg'] + more_outcome

cluster_name = 'cs_class'

# compute median cluster consensus for all cluster result
df_medCC_all, df_oi = dh.Median_Cluster_Consensus(im_dir, cc_prefix, the_N_cluster)

#TODO: find the CC setting with the highest CC (not lowest p-value..)!!!

chi2_all_df = pd.DataFrame(columns=('cluster_method','cluster_linkage','dist_method','N_mincluster','N_sample','N_cluster','v1','v2','chi2','pval','cramersV','medianCC'))
count = 0
# df_tmp = df_medCC_all[df_medCC_all['k'] == the_N_cluster]
df_tmp = df_medCC_all[df_medCC_all['k'].isin(the_N_cluster)]
# df_tmp = df_oi

for ii in range(df_tmp.shape[0]):
    k, the_cm, the_cl, the_dm, medCC = df_tmp.ix[df_tmp.index[ii], ['k','cluster_method', 'cluster_linkage', 'dist_method', 'medianCC']].values
    if isinstance(the_cl, str):
        cc_dir = '{}/{}_{}_{}_{}'.format(im_dir, cc_prefix, the_cm, the_cl, the_dm)
    else:
        cc_dir = '{}/{}_{}_{}'.format(im_dir, cc_prefix, the_cm, the_dm)

    cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir,k)
    cs_class_df = pd.read_csv(cs_class_fname)
    cs_class_df.columns = ['ptid_side', 'cs_class']

    # combine the cs_class to the_df
    the_df = pd.merge(jdf, cs_class_df, on='ptid_side')

    N_mincluster = the_df['cs_class'].value_counts().min()
    for v2 in outcome_name_list:
        # print 'v2 is {}'.format(v2)
        chi2,pval,dum1,dum2,N_sample,cv,tab = st.chi2test(the_df,cluster_name,v2)
        chi2_all_df.loc[count] = [the_cm, the_cl, the_dm, N_mincluster, N_sample, k, cluster_name, v2, chi2, pval,cv, medCC]
        count = count + 1

# the_ov_interest = 'T_stage'
# print chi2_all_df[chi2_all_df['v2'] == the_ov_interest]

# save the chi2 result for all data
chi2_fname = '{}/chi2_Ncluster{}_Rep{}.csv'.format(im_dir, '_'.join(str(x) for x in the_N_cluster), Nrep)
chi2_all_df.to_csv(chi2_fname)

# get the optimal cluster setting for each outcome variable and plot the clustermap
fdf = chi2_all_df[chi2_all_df['N_mincluster'] > 5]
# idx = fdf.groupby('v2').apply(lambda df: df.pval.argmin())
idx = fdf.groupby('v2').apply(lambda df: df.cramersV.argmax())
# idx = fdf.groupby('v2').apply(lambda df: df.medianCC.argmax())

# determine the number of outcome variables
the_outcome_vars_list = chi2_all_df['v2'].unique().tolist()
# the_outcome_vars_list = ['BC_subtype']

imf_pat = re.compile('texture_|ShapeSize_|FOstats_')
# imf_petmri_names = [ss for ss in jdf.columns.tolist() if imf_pat.match(ss)] + ['MRI_SERROI_SER_MEAN','MRI_PE2ROI_PE2_PEAK','PET_SUV_max']
imf_petmri_names = [ss for ss in jdf.columns.tolist() if imf_pat.match(ss)] + ['MRI_SERROI_SER_MEAN','MRI_PE2ROI_PE2_PEAK']
# imf_petmri_names = [ss + '_pet' for ss in img_feature_names] + [ss + '_mri' for ss in img_feature_names]

df_chi2_outcome_all = pd.DataFrame(columns=('outcome_var','cluster_method','cluster_linkage','dist_method','N_mincluster','N_sample','N_cluster','pval','cramersV','medianCC'))
count = 0
for ov in the_outcome_vars_list:
    ifdf = idx[ov]
    # print 'outcome var: {}'.format(ov)
    cm, cl, dm, N_mincluster, Nsample, N_cluster, pval, cv, medCC = fdf.ix[ifdf,['cluster_method','cluster_linkage','dist_method','N_mincluster','N_sample','N_cluster','pval','cramersV','medianCC']].tolist()
    if isinstance(cl, str):
        cc_dir = '{}/{}_{}_{}_{}'.format(im_dir, cc_prefix, cm,cl,dm)
    else:
        cc_dir = '{}/{}_{}_{}'.format(im_dir, cc_prefix, cm, dm)
    print(cm,cl,dm)

    cs_class_fname = '{}/ConsensusClass_kmax{:.0f}.csv'.format(cc_dir, N_cluster)
    cs_class_df = pd.read_csv(cs_class_fname)
    cs_class_df.columns = ['ptid_side', 'cs_class']

    # combine the cs_class to the_df
    the_df = pd.merge(jdf, cs_class_df, on='ptid_side')

    df_chi2_outcome_all.loc[count] = [ov, cm, cl, dm, N_mincluster, Nsample, N_cluster, pval, cv, medCC]
    count = count + 1

    # print proportion table for chi2 test
    _,_,_ = st.proportion_table(the_df, 'cs_class', ov, table_margin=True)

    if ov == 'Recurrence':
        fig_fname = '{}/clustermap_Recur_kmax{:.0f}.pdf'.format(cc_dir, N_cluster)
        # recur_pal1 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3], sns.color_palette('Set2', 10)[7]]
        recur_pal1 = sns.color_palette('YlOrRd', len(the_df[ov].unique()))
        vp.ClustermapPlot(the_df, [ov], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr, idvar_color_pal_mode=2,var_color_pal=[recur_pal1], var_title=['Recurrence Status'])
    elif ov == 'BoneMetsOrNot':
        fig_fname = '{}/clustermap_BoneMetsOrNot_kmax{:.0f}.pdf'.format(cc_dir, N_cluster)
        the_pal = [(1., 1., 1.), sns.color_palette('Set2', 10)[3]]
        vp.ClustermapPlot(the_df, [ov], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr, idvar_color_pal_mode=2, var_title=['Bone Mets'],var_color_pal=[the_pal])
    elif re.search('OS_(.+)', ov) or re.search('DF_(.+)',ov):
        fig_fname = '{}/clustermap_{}_kmax{:.0f}.pdf'.format(cc_dir, ov, N_cluster)
        the_pal = [(1., 1., 1.), sns.color_palette('Set2', 10)[3]]
        print(ov, fig_fname)
        vp.ClustermapPlot(the_df, [ov], fig_fname, imf_petmri_names, vp.featurename_redef_petmr, vp.featurelabel_redef_petmr, idvar_color_pal_mode=2,var_color_pal=[the_pal])
    elif ov == 'BC_subtype':
        fig_fname = '{}/clustermap_BCsubtype_kmax{:.0f}.pdf'.format(cc_dir, N_cluster)
        vp.ClustermapPlot(the_df,[ov],fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Breast cancer subtype'])
    elif ov == 'T_stage':
        fig_fname = '{}/clustermap_Tstage_kmax{:.0f}.pdf'.format(cc_dir, N_cluster)
        vp.ClustermapPlot(the_df,[ov],fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['T-stage'])
    elif ov == 'N_stage':
        fig_fname = '{}/clustermap_Nstage_kmax{:.0f}.pdf'.format(cc_dir, N_cluster)
        vp.ClustermapPlot(the_df, [ov], fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['N-stage'])
    elif ov == 'Overall_stage':
        fig_fname = '{}/clustermap_Overallstage_kmax{:.0f}.pdf'.format(cc_dir, N_cluster)
        the_pal = sns.color_palette('Purples', len(the_df[ov].unique()))
        vp.ClustermapPlot(the_df, [ov], fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Overall stage'], var_color_pal=[the_pal])
    elif ov == 'Overall_stage_2':
        fig_fname = '{}/clustermap_Overallstage2_kmax{:.0f}.pdf'.format(cc_dir, N_cluster)
        vp.ClustermapPlot(the_df, [ov], fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Overall stage'])
    elif ov == 'TripleNeg':
        fig_fname = '{}/clustermap_TripleNeg_kmax{:.0f}.pdf'.format(cc_dir, N_cluster)
        tripneg_pal = [(1., 1., 1.),sns.color_palette('Set2', 10)[2]]
        vp.ClustermapPlot(the_df, [ov], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Triple Neg'], var_color_pal=[tripneg_pal])
    elif ov == 'Tumor_Grade':
        fig_fname = '{}/clustermap_Grade_kmax{:.0f}.pdf'.format(cc_dir, N_cluster)
        light_pal = sns.light_palette((210, 90, 60), input="husl")
        grade_pal = light_pal[1:len(light_pal):2]
        vp.ClustermapPlot(the_df, [ov], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Tumor Grade'], var_color_pal=[grade_pal])
    elif ov == 'Tumor_Histology_updated':
        fig_fname = '{}/clustermap_Histology_kmax{:.0f}.pdf'.format(cc_dir, N_cluster)
        colors = ['windows blue', 'amber', 'greyish', 'brick']
        histology_pal = sns.xkcd_palette(colors)
        vp.ClustermapPlot(the_df, [ov], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Tumor Histology'], var_color_pal=[histology_pal])
    elif ov == 'Recurrence_Type':
        fig_fname = '{}/clustermap_recurtype_kmax{:.0f}.pdf'.format(cc_dir, N_cluster)
        recurtype_pal = sns.light_palette('navy', n_colors=len(the_df['Recurrence_Type'].unique()), reverse=True)
        vp.ClustermapPlot(the_df, [ov], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Recurrence Type'],var_color_pal=[recurtype_pal])

    # elif ov == 'Recurrence_CR_1':
    #     fig_fname = '{}/clustermap_recurstatus1_kmax{}.pdf'.format(cc_dir, the_N_cluster)
    #     recur_pal1 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3], sns.color_palette('Set2', 10)[7]]
    #     vp.ClustermapPlot(the_df, [ov], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,
    #                       vp.featurelabel_redef_petmr, idvar_color_pal_mode=2, var_title=['Recurrence'],
    #                       var_color_pal=[recur_pal1])
    # elif ov == 'Recurrence_CR_2':
    #     fig_fname = '{}/clustermap_recurstatus2_kmax{}.pdf'.format(cc_dir, the_N_cluster)
    #     recur_pal2 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3]]
    #     vp.ClustermapPlot(the_df, [ov], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,
    #                       vp.featurelabel_redef_petmr, idvar_color_pal_mode=2, var_title=['Recurrence'],
    #                       var_color_pal=[recur_pal2])
    # elif ov == 'Recurrence_3':
    #     fig_fname = '{}/clustermap_recur3_kmax{}.pdf'.format(cc_dir, the_N_cluster)
    #     recur_pal1 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3], sns.color_palette('Set2', 10)[7]]
    #     vp.ClustermapPlot(the_df, [ov], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,
    #                       vp.featurelabel_redef_petmr, idvar_color_pal_mode=2, var_title=['Recurrence'],
    #                       var_color_pal=[recur_pal1])
    # if ov == 'Recurrence_4':
    #     fig_fname = '{}/clustermap_recur4_kmax{}.pdf'.format(cc_dir, the_N_cluster)
    #     recur_pal1 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3]]
    #     vp.ClustermapPlot(the_df, [ov], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,
    #                       vp.featurelabel_redef_petmr, idvar_color_pal_mode=2, var_title=['Recurrence'],
    #                       var_color_pal=[recur_pal1])

# print df_chi2_outcome_all
# chi2_outcome_fname = '{}/BestChi2_outcome_Ncluster{}.csv'.format(im_dir, '_'.join(str(x) for x in the_N_cluster))
chi2_outcome_fname = '{}/BestCVChi2_outcome_Ncluster{}_Rep{}.csv'.format(im_dir, '_'.join(str(x) for x in the_N_cluster), Nrep)
# chi2_outcome_fname = '{}/BestMedCCChi2_outcome_Ncluster{}_3.csv'.format(im_dir, '_'.join(str(x) for x in the_N_cluster))
df_chi2_outcome_all.to_csv(chi2_outcome_fname)


# # print the number break down for DF_xxxx
# more_outcome = ['DF_{:.0f}yr'.format(x) for x in n_yr_list]
# for ov in more_outcome:
#     print('outcome: {}'.format(ov))
#     print(jdf[ov].value_counts())
#
# # # combine the cs_class to jdf (combined PET and MRI feature data)
# # # determine the optimal clustering setting for a given N_cluster (determined by maximum median cluster consensus)
# # idx = df_medCC_all.groupby('k').apply(lambda df: df.medianCC.argmax())
# # df_tmp = df_medCC_all.loc[idx,:]
# # the_N_cluster = 3
# # the_cm, the_cl, the_dm, the_medCC = df_tmp.ix[df_tmp['k'] == the_N_cluster,['cluster_method','cluster_linkage','dist_method','medianCC']].values[0]
# # # print the_cm, the_cl, the_dm,the_medCC
# #
# #
# # # obtain the cluster class result
# # the_cm = 'pam'
# # the_dm = 'spearman'
# # the_cl = np.nan
# # if isinstance(the_cl, str):
# #     cc_dir = '{}/ConsensusCluster_{}_{}_{}'.format(im_dir, the_cm, the_cl, the_dm)
# # else:
# #     cc_dir = '{}/ConsensusCluster_{}_{}'.format(im_dir, the_cm, the_dm)
# #
# # cc_file = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir,the_N_cluster)
# # print 'read the file {}'.format(cc_file)
# # df_cc = pd.read_csv(cc_file)
# # df_cc.columns = ['ptid_side', 'cs_class']
# # the_df = pd.merge(jdf, df_cc, on='ptid_side')
#
# # # save dataframe for R testing
# # fname = '{}/data_cs_all.csv'.format(cc_dir)
# # the_df.to_csv(fname)
#
# # # # Plot clustermap of radiomics vs. clinical outcome for samples
# # imf_petmri_names = [ss + '_pet' for ss in img_feature_names] + [ss + '_mri' for ss in img_feature_names]
# #
# # # fig_fname = '{}/clustermap_BCsubtype_kmax{}.pdf'.format(cc_dir,the_N_cluster)
# # # vp.ClustermapPlot(the_df,['BC_subtype'],fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Breast cancer subtype'])
# #
# # fig_fname = '{}/clustermap_Tstage_kmax{}.pdf'.format(cc_dir,the_N_cluster)
# # vp.ClustermapPlot(the_df,['T_stage'],fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['T-stage'])
#
# # fig_fname = '{}/clustermap_Nstage_kmax{}.pdf'.format(cc_dir, the_N_cluster)
# # vp.ClustermapPlot(the_df, ['N_stage'], fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title='N-stage')
# #
# # fig_fname = '{}/clustermap_Overallstage_kmax{}.pdf'.format(cc_dir, the_N_cluster)
# # vp.ClustermapPlot(the_df, ['Overall_stage'], fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title='Overall stage')
# #
# # fig_fname = '{}/clustermap_TripleNeg_kmax{}.pdf'.format(cc_dir, the_N_cluster)
# # tripneg_pal = [(1., 1., 1.),sns.color_palette('Set2', 10)[2]]
# # vp.ClustermapPlot(the_df, ['TripleNeg'], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Triple Neg'], var_color_pal=[tripneg_pal])
# #
# # fig_fname = '{}/clustermap_Grade_kmax{}.pdf'.format(cc_dir, the_N_cluster)
# # light_pal = sns.light_palette((210, 90, 60), input="husl")
# # grade_pal = light_pal[1:len(light_pal):2]
# # vp.ClustermapPlot(the_df, ['Sjoerd_Grade'], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Grade'], var_color_pal=[grade_pal])
# #
# # fig_fname = '{}/clustermap_Histology_kmax{}.pdf'.format(cc_dir, the_N_cluster)
# # colors = ['windows blue', 'amber', 'greyish', 'brick']
# # histology_pal = sns.xkcd_palette(colors)
# # vp.ClustermapPlot(the_df, ['Marjan_Histology'], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Histology'],var_color_pal=[histology_pal])
# #
# # fig_fname = '{}/clustermap_recurstatus1_kmax{}.pdf'.format(cc_dir, the_N_cluster)
# # recur_pal1 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3],sns.color_palette('Set2', 10)[7]]
# # vp.ClustermapPlot(the_df, ['Recurrence_CR_1'], fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Recurrence'],var_color_pal=[recur_pal1])
# #
# # fig_fname = '{}/clustermap_recurstatus2_kmax{}.pdf'.format(cc_dir, the_N_cluster)
# # recur_pal2 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3]]
# # vp.ClustermapPlot(the_df, ['Recurrence_CR_2'], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Recurrence'],var_color_pal=[recur_pal2])
# #
# # fig_fname = '{}/clustermap_df_5yr_kmax{}.pdf'.format(cc_dir, the_N_cluster)
# # df5yr_pal = [(1., 1., 1.), sns.color_palette('Set2', 10)[3]]
# # vp.ClustermapPlot(the_df, ['Diseasefree_5yr'], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['5-yr Disease Free'],
# #                var_color_pal=[df5yr_pal])
# #
# # fig_fname = '{}/clustermap_recurtype_kmax{}.pdf'.format(cc_dir, the_N_cluster)
# # recurtype_pal = sns.light_palette('navy', n_colors=len(the_df['Recurrence_Type'].unique()), reverse=True)
# # vp.ClustermapPlot(the_df, ['Recurrence_Type'], fig_fname, imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2, var_title=['Recurrence Type'],var_color_pal=[recurtype_pal])
# #
# # fig_fname = '{}/clustermap_df5yr_TN_kmax{}.pdf'.format(cc_dir, the_N_cluster)
# # varname = ['Diseasefree_5yr','TripleNeg']
# # vartitle = ['5-yr Disease free','Triple Negative']
# # pal1  = [(1., 1., 1.), sns.color_palette('Set2', 10)[3]]
# # pal2 = [(1., 1., 1.), sns.color_palette('Set2', 10)[8]]
# # var_pal_list = [pal1,pal2]
# # vp.ClustermapPlot(the_df,varname, fig_fname,imf_petmri_names, vp.featurename_redef_petmr,vp.featurelabel_redef_petmr,idvar_color_pal_mode=2,var_title=vartitle,var_color_pal=var_pal_list)
#
#
# # # perform chi2 test for CC clusters vs clinical outcome
# # chi2_all_df = pd.DataFrame(columns=('N_sample','v1','v2','chi2','pval'))
# # count = 0
# # for v2 in outcome_name_list:
# #     # print 'v2 is {}'.format(v2)
# #     chi2,pval,dum1,dum2,N_sample = st.chi2test(the_df,cluster_name,v2)
# #     # print chi2,pval, N_sample
# #     chi2_all_df.loc[count] = [N_sample,cluster_name,v2,chi2,pval]
# #     count = count + 1
# #
# # print chi2_all_df

