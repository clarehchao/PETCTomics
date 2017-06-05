#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 8/23/16 1:59 PM

@author: shuang Shih-ying Huang
@goal: plot a cluster plot image features vs. PRIMARY_ID to get an idea of the cluster distribution

"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import re
import numpy as np
import os
import errno

def makedirs(file_path):
    if os.path.exists(file_path):
        if not os.path.isdir(file_path):
            raise IOError('O_O: {} is not a directory!'.format(file_path))
    else:
        os.makedirs(file_path)

rootdir = '/data/francgrp1/clare_work/Data'
feature_data_dir = '{}/her2_ImageFeatures/IsoVoxelSize'.format(rootdir)
pt_data_dir = '{}/her2_ClinicalData'.format(rootdir)

# MRI_PET_ImgFeature json
fname = '{}/MRI_PET_imgfeat_data.json'.format(pt_data_dir)
df_mri_pet_imgfeat = pd.read_json(fname,dtype={'pt_mrn':str,'pt_id':str,'PRIMARY_ID':str,'MRN':str})
df_mri_pet_imgfeat['MRN'] = df_mri_pet_imgfeat['MRN'].apply(lambda x: '{0:0>8}'.format(x))

# fill in the clinical data like ER, HER2 HR status with nan for '8888' or '9999'
df_mri_pet_imgfeat['Sjoerd_ER'] = df_mri_pet_imgfeat['Sjoerd_ER'].replace([8888,9999],np.nan)
df_mri_pet_imgfeat['Sjoerd_HER2'] = df_mri_pet_imgfeat['Sjoerd_HER2'].replace([8888,9999],np.nan)
df_mri_pet_imgfeat['Sjoerd_PR'] = df_mri_pet_imgfeat['Sjoerd_PR'].replace([8888,9999],np.nan)
df_mri_pet_imgfeat['Sjoerd_Grade'] = df_mri_pet_imgfeat['Sjoerd_Grade'].replace([8888,9999],np.nan)

# pick the needed columns for cluster analysis
clinical_colname = ['Laterality','PRIMARY_ID','MRN','CR_AccessionSeq','MRI_VOLUME_TUM_BLU','MRI_VOLUME_TUMOR','MRI_VOLUME_WHITE','MRI_VOLUME_RED','MRI_VOLUME_GREEN',
                    'MRI_VOLUME_PURPLE','MRI_VOLUME_BLUE','MRI_SERROI_SER_MEAN','MRI_PE2ROI_PE2_PEAK','MRI_SERROI_PE2_PEAK','MRI_PE2ROI_SER_MEAN','SUV max',
                    'Marjan_Histology','Marjan_Size (mm)','Sjoerd_HER2','Sjoerd_ER','Sjoerd_PR','Sjoerd_Grade']
df_clinicaldata = df_mri_pet_imgfeat.ix[:,clinical_colname]

# read all the Image Feature data
# all_jsons = glob.glob('{}/*.json'.format(feature_data_dir))
all_jsons = glob.glob('{}/MRI*.json'.format(feature_data_dir))
df_img_features = pd.DataFrame()
for fj in all_jsons:
    df = pd.read_json(fj,dtype={'pt_id':str,'pt_mrn':str}) #make sure pt_id and MRN are read in as string
    df['dce_MRI_time_point'] = df['dce_series_dmi_fname'].apply(lambda x: re.search(r'(.+)tp(\d+).dmi',x).group(2))
    df['Laterality'] = df['dce_series_dmi_fname'].apply(lambda x: re.search(r'(.+)sag(\w+)_tp(\d+).dmi',x).group(2))
    df_img_features = df_img_features.append(df,ignore_index=True)

df_img_clinical_data = pd.merge(df_img_features,df_clinicaldata,how='left',left_on=['pt_mrn','Laterality'],right_on=['MRN','Laterality'])

# select the data of interest
# clustermap for each time point feature data
# # the_tp = '1'
# the_tp = '2'
# the_glcm_Nbin = 128

# drop un-needed columns
features_drop_list = ['CR_AccessionSeq','MRI_VOLUME_TUM_BLU','MRI_VOLUME_TUMOR','MRI_VOLUME_WHITE','MRI_VOLUME_RED',
                     'MRI_VOLUME_GREEN','MRI_VOLUME_PURPLE','MRI_VOLUME_BLUE','MRI_SERROI_SER_MEAN','MRI_PE2ROI_PE2_PEAK',
                     'MRI_SERROI_PE2_PEAK','MRI_PE2ROI_SER_MEAN','Marjan_Size (mm)','organ_mask', 'process_name', 'process_version',
                     'texture_glcm_offset', 'dce_MRI_time_point','dce_series_dmi_fname', 'glcm_Nbin', 'img_norm_Nbin', 'pt_id',
                     'pt_accession_num', 'pt_mrn','SUV max']

# data munging: get the average of all the features values (over all offsets) and split minmax into separate variables
texture_cols = ['texture_autocorrelation', 'texture_cluster_prominence', 'texture_cluster_shade',
                'texture_cluster_tendency', 'texture_contrast', 'texture_correlation', 'texture_diff_entropy',
                'texture_dissimilarity','texture_energy', 'texture_entropy', 'texture_homogeneity1', 'texture_homogeneity2',
                'texture_idmn','texture_idn', 'texture_inv_var', 'texture_maxprob', 'texture_sum_avg', 'texture_sum_entropy',
                'texture_sum_var']

the_tp_list = df_img_features['dce_MRI_time_point'].unique().tolist()
the_glcmNbin_list = df_img_features['glcm_Nbin'].unique().tolist()

for the_tp in the_tp_list:
    for the_glcm_Nbin in the_glcmNbin_list:
        df_tmp = df_img_clinical_data[(df_img_clinical_data['dce_MRI_time_point'] == the_tp) & (df_img_clinical_data['glcm_Nbin'] == the_glcm_Nbin)]

        # combine primary_id with breast_side
        # Need to make a copy of df_tmp so SettingWithCopyException doesn't appear
        df_tp_glcmNbin = df_tmp.copy()

        # NOTE: must use 'pt_id' not 'PRIMARY_ID' since 'PRIMARY_ID' numbering is different from the pt_id
        df_tp_glcmNbin['ptid_side'] = df_tmp.apply(lambda x: '{}_{}'.format(x['pt_id'],x['Laterality']),axis=1)
        # print df_tp_glcmNbin.loc[df_tp_glcmNbin['pt_id'] == '96','ptid_side']

        # drop un-needed features for optimal N for future analysis
        df_tp_glcmNbin = df_tp_glcmNbin.drop(features_drop_list,axis=1)

        for tc in texture_cols:
            df_tp_glcmNbin[tc +'_avg'] = df_tp_glcmNbin[tc].apply(np.mean)
            df_tp_glcmNbin = df_tp_glcmNbin.drop(tc,axis=1)

        df_tp_glcmNbin['FOstats_min'] = df_tp_glcmNbin['FOstats_minmax'].apply(lambda x: x[0])
        df_tp_glcmNbin['FOstats_max'] = df_tp_glcmNbin['FOstats_minmax'].apply(lambda x: x[1])
        df_tp_glcmNbin = df_tp_glcmNbin.drop('FOstats_minmax',axis=1)

        # # drop any NaN
        # print df_tp_glcmNbin.shape
        # test = df_tp_glcmNbin[df_tp_glcmNbin.isnull().any(axis=1)]
        # print test.ix[31,:].tolist()

        df_tp_glcmNbin = df_tp_glcmNbin.dropna()
        print 'the_tp: {}, glcm_Nbin: {}, df.shape: {}'.format(the_tp, the_glcm_Nbin, df_tp_glcmNbin.shape)

        # separate the patient clinical data into Triple-Neg vs Non Triple-Neg
        df_tp_glcmNbin['TripleNeg'] = df_tp_glcmNbin.apply(lambda x: 1.0 if x['Sjoerd_ER'] == 0.0 and x['Sjoerd_HER2'] == 0.0 and x['Sjoerd_PR'] == 0.0 else 0.0, axis=1)
        fdir = '{}/her2_Analysis/MRI/IsoVoxel_TP{}_GLCMBIN{}'.format(rootdir,the_tp,the_glcm_Nbin)
        makedirs(fdir)
        tp_glcmNbin_fname = '{}/MRIdataAll_tp{}_Nbin{}.csv'.format(fdir,the_tp,the_glcm_Nbin)
        df_tp_glcmNbin.to_csv(tp_glcmNbin_fname)


# df_tp_glcmNbin.sort_values('TripleNeg',inplace=True)
# the_TripleNegsorted_ptidside = df_tp_glcmNbin['ptid_side'].tolist()
#
# # re-format the radiomic features into columnes of feature name and feature values
# img_feature_names = ['MRI_PE2ROI_PE2_PEAK','MRI_PE2ROI_SER_MEAN','MRI_SERROI_PE2_PEAK','MRI_SERROI_SER_MEAN','MRI_VOLUME_BLUE','MRI_VOLUME_GREEN','MRI_VOLUME_PURPLE','MRI_VOLUME_RED','MRI_VOLUME_TUMOR','MRI_VOLUME_TUM_BLU','MRI_VOLUME_WHITE','FOstats_energy',
#                      'FOstats_entropy','FOstats_kurtosis','FOstats_mean','FOstats_min','FOstats_max','FOstats_skewness','FOstats_uniformity','FOstats_variance','ShapeSize_compactness1','ShapeSize_compactness2','ShapeSize_max_euc_dis',
#                      'ShapeSize_spherical_disproportion','ShapeSize_sphericity','ShapeSize_surf_area_cm2','ShapeSize_surface2volratio','ShapeSize_vol_cm3'] + [ss + '_avg' for ss in texture_cols]
#
#
# df_flatten_data = pd.melt(df_tp_glcmNbin,id_vars='ptid_side',value_vars=img_feature_names,var_name='Image Feature')
# df_pivot = df_flatten_data.pivot(index='Image Feature',columns='ptid_side',values='value')
# df_pivot_sorted = df_pivot[the_TripleNegsorted_ptidside]
#
# TripleNeg_status = df_tp_glcmNbin['TripleNeg'].tolist()
# tumor_grade = df_tp_glcmNbin['Sjoerd_Grade'].tolist()
# histology = df_tp_glcmNbin['Marjan_Histology'].tolist()
#
# # build in tumor status for each primary_ID
# mcol = pd.MultiIndex.from_tuples(zip(histology,tumor_grade,TripleNeg_status),names=['Histology','Grade','TripleNeg'])
# df_pivot_sorted.columns = mcol
#
# # cluster colorbar
# img_feature_vars = df_pivot_sorted.index.tolist()
#
# # categorical label for the image feature vars
# img_feature_labels = []
# for ss in img_feature_vars:
#     if re.search(r'FOstats_(.+)',ss):
#         img_feature_labels.append('First-order Stats')
#     elif re.search(r'MRI_(.+)',ss):
#         img_feature_labels.append('Marjan MRI features')
#     elif re.search(r'ShapeSize_(.+)',ss):
#         img_feature_labels.append('Shape and Size')
#     elif re.search(r'texture_(.+)',ss):
#         img_feature_labels.append('GLCM texture features')
#
# img_feature_pal = sns.light_palette('navy',n_colors=len(set(img_feature_labels)),reverse=True)
# img_feature_lut = dict(zip(map(str, list(set(img_feature_labels))), img_feature_pal))
# img_feature_colors = pd.Series(img_feature_labels,index=df_pivot.index).map(img_feature_lut)
# df_img_feature_colors = pd.DataFrame(dict(ImageFeature=img_feature_colors))
# df_img_feature_colors.columns = ['']
#
#  # categorize the clinical data into sub-group for visualization (set col_colors for the clustermap
# clinical_data_colors = pd.DataFrame()
# the_pal = sns.color_palette('Set2',10)
# df_pivot_col_unique = df_pivot_sorted.columns.unique()
# TripleNeg_colors_labels = [df_pivot_col_unique[i][2] for i in range(len(df_pivot_col_unique))]
# clinical_data_colors['TripleNeg'] = pd.Series(TripleNeg_colors_labels,index=df_pivot_col_unique).map({0:(1.,1.,1.),1:the_pal[2]})
#
# light_pal = sns.light_palette((210, 90, 60), input="husl")
# grade_pal = light_pal[1:len(light_pal):2]
# grade_labels = [df_pivot_col_unique[i][1] for i in range(len(df_pivot_col_unique))]
# unique_grade_labels = list(set(grade_labels))
# grade_lut = dict(zip(unique_grade_labels, grade_pal))
# clinical_data_colors['Grade'] = pd.Series(grade_labels,index=df_pivot_col_unique).map(grade_lut)
#
# colors = ['windows blue', 'amber', 'greyish','brick']
# histology_pal = sns.xkcd_palette(colors)
# histology_labels = [df_pivot_col_unique[i][0] for i in range(len(df_pivot_col_unique))]
# unique_histology_labels = list(set(histology_labels))
# histology_lut = dict(zip(unique_histology_labels, histology_pal))
# clinical_data_colors['Histology'] = pd.Series(histology_labels,index=df_pivot_col_unique).map(histology_lut)
#
# # Create a custom colormap for the heatmap values
# cmap = sns.diverging_palette(h_neg=10, h_pos=220, s=80, n=7, as_cmap=True)
#
# # go through different methods and metrics for plotting patients separated into triple neg and non triple-negative
# sns.set()
# sns.set(font="monospace")
# cluster_method = ['average','complete','weighted']  # try 'average, 'complete' (max of distance), 'weighted'; method='centroid','ward', 'median' didn't work
# cluster_metric = ['euclidean','minkowski','cityblock','seuclidean','sqeuclidean','cosine','correlation','chebyshev','canberra','braycurtis']
#
# for mt in cluster_method:
#     for mm in cluster_metric:
#         try:
#             g = sns.clustermap(df_pivot_sorted, method=mt, metric=mm, col_cluster=False,
#                                col_colors=clinical_data_colors, row_colors=df_img_feature_colors,
#                                z_score=0, cmap=cmap, linewidths=0, xticklabels=False, yticklabels=False,
#                                figsize=(15, 15))
#             # display the legend for the image feature colormap by adding a empty bar plot but display the legend
#             for label in list(set(img_feature_labels)):
#                 g.ax_row_dendrogram.bar(0, 0, color=img_feature_lut[label], label=label, linewidth=0)
#             g.ax_row_dendrogram.legend(bbox_to_anchor=(1.4, 1.2), loc='best', title='Image Features')
#
#             for label in list(set(histology_labels)):
#                 g.ax_col_dendrogram.bar(0, 0, color=histology_lut[label], label=label, linewidth=0)
#             g.ax_col_dendrogram.legend(bbox_to_anchor=(0.95, 1.62), loc='best', title='Tumor histology')
#
#             for label in list(set(grade_labels)):
#                 g.ax_heatmap.bar(0, 0, color=grade_lut[label], label=label, linewidth=0)
#             g.ax_heatmap.legend(bbox_to_anchor=(0.5, 1.34), loc='best', title='Tumor grade')
#
#             # position the heatmap colorbar appropriately
#             g.cax.set_position([0.92, .2, .03, .45])
#             g.ax_heatmap.set(xlabel='', ylabel='')
#
#             g.fig.suptitle('cluster method = {}, metric = {}'.format(mt, mm), fontweight='bold')
#
#             # create the file-save directories
#             fw_dir = '{}/her2_Analysis/MRI/IsoVoxel_TP{}_GLCMBIN{}'.format(rootdir,the_tp,the_glcm_Nbin)
#             try:
#                 os.makedirs(fw_dir)
#             except OSError as exception:
#                 if not os.path.isdir(fw_dir):
#                     raise
#             clusterplotfname = '{}/clustermap_TripleNeg_{}_{}.pdf'.format(fw_dir,mt,mm)
#             g.savefig(clusterplotfname)
#
#             plt.show()
#
#         except:
#             print 'error in cluster method: {}, metric: {}'.format(mt,mm)
#             continue

# Plot a clustermap plot for all patients with all clinical info
# # re-format the radiomic features into columnes of feature name and feature values
# pt_info_names = ['PRIMARY_ID','MRN','CR_AccessionSeq']
# clinical_feature_names = ['Marjan_Histology','Marjan_Size (mm)','Sjoerd_ER','Sjoerd_HER2','Sjoerd_PR','Sjoerd_Grade']
# img_feature_names = ['MRI_PE2ROI_PE2_PEAK','MRI_PE2ROI_SER_MEAN','MRI_SERROI_PE2_PEAK','MRI_SERROI_SER_MEAN','MRI_VOLUME_BLUE','MRI_VOLUME_GREEN','MRI_VOLUME_PURPLE','MRI_VOLUME_RED','MRI_VOLUME_TUMOR','MRI_VOLUME_TUM_BLU','MRI_VOLUME_WHITE','FOstats_energy',
#                      'FOstats_entropy','FOstats_kurtosis','FOstats_mean','FOstats_min','FOstats_max','FOstats_skewness','FOstats_uniformity','FOstats_variance','ShapeSize_compactness1','ShapeSize_compactness2','ShapeSize_max_euc_dis',
#                      'ShapeSize_spherical_disproportion','ShapeSize_sphericity','ShapeSize_surf_area_cm2','ShapeSize_surface2volratio','ShapeSize_vol_cm3'] + [ss + '_avg' for ss in texture_cols]
#
#
# df_flatten_data = pd.melt(df_tp_glcmNbin,id_vars='PRIMARY_ID',value_vars=img_feature_names,var_name='Image Feature')
# df_pivot = df_flatten_data.pivot(index='Image Feature',columns='PRIMARY_ID',values='value')
#
# # index df by 'PRIMARY_ID' for ease of data selection
# df_tp_glcmNbin.index = df_tp_glcmNbin['PRIMARY_ID'].tolist()
#
# # numerate the histology column for multi-indexing purpose
# the_histology_labels = df_tp_glcmNbin['Marjan_Histology'].unique()
# histology_map = dict(zip(list(set(the_histology_labels)), range(len(the_histology_labels))))
# the_grade_labels = df_tp_glcmNbin['Sjoerd_Grade'].unique()
# grade_map = dict(zip(list(set(the_grade_labels)), range(len(the_grade_labels))))
# df_tp_glcmNbin['Marjan_Histology_category'] = df_tp_glcmNbin['Marjan_Histology'].apply(lambda x: histology_map[x])
# df_tp_glcmNbin['Sjoerd_Grade_category'] = df_tp_glcmNbin['Sjoerd_Grade'].apply(lambda x: grade_map[x])
#
# the_id_list = df_pivot.columns.tolist()
# ER_status = df_tp_glcmNbin.loc[the_id_list,'Sjoerd_ER'].tolist()
# HER2_status = df_tp_glcmNbin.loc[the_id_list,'Sjoerd_HER2'].tolist()
# PR_status = df_tp_glcmNbin.loc[the_id_list,'Sjoerd_PR'].tolist()
# # tumor_grade = df_tp_glcmNbin.loc[the_id_list,'Sjoerd_Grade_category'].tolist()
# # histology = df_tp_glcmNbin.loc[the_id_list,'Marjan_Histology_category'].tolist()
# tumor_grade = df_tp_glcmNbin.loc[the_id_list,'Sjoerd_Grade'].tolist()
# histology = df_tp_glcmNbin.loc[the_id_list,'Marjan_Histology'].tolist()
#
# # build in tumor status for each primary_ID
# # mcol = pd.MultiIndex(levels=[the_histology_labels,the_grade_labels,[0,1],[0,1],[0,1]],labels=[histology,tumor_grade,ER_status,HER2_status,PR_status],names=['Histology','Grade','ER','HER2','PR'])
# mcol = pd.MultiIndex.from_tuples(zip(histology,tumor_grade,ER_status,HER2_status,PR_status),names=['Histology','Grade','ER','HER2','PR'])
# df_pivot.columns = mcol
#
# # cluster colorbar
# img_feature_vars = df_pivot.index.tolist()
#
# # categorical label for the image feature vars
# img_feature_labels = []
# for ss in img_feature_vars:
#     if re.search(r'FOstats_(.+)',ss):
#         img_feature_labels.append('First-order Stats')
#     elif re.search(r'MRI_(.+)',ss):
#         img_feature_labels.append('Marjan MRI features')
#     elif re.search(r'ShapeSize_(.+)',ss):
#         img_feature_labels.append('Shape and Size')
#     elif re.search(r'texture_(.+)',ss):
#         img_feature_labels.append('GLCM texture features')
#
# img_feature_pal = sns.light_palette('navy',n_colors=len(set(img_feature_labels)),reverse=True)
# img_feature_lut = dict(zip(map(str, list(set(img_feature_labels))), img_feature_pal))
# img_feature_colors = pd.Series(img_feature_labels,index=df_pivot.index).map(img_feature_lut)
# df_img_feature_colors = pd.DataFrame(dict(ImageFeature=img_feature_colors))
# df_img_feature_colors.columns = ['']
#
#  # categorize the clinical data into sub-group for visualization (set col_colors for the clustermap
# clinical_data_colors = pd.DataFrame()
# the_pal = sns.color_palette('Set2',10)
# df_pivot_col_unique = df_pivot.columns.unique()
# ER_colors_labels = [df_pivot_col_unique[i][2] for i in range(len(df_pivot_col_unique))]
# clinical_data_colors['ER'] = pd.Series(ER_colors_labels,index=df_pivot_col_unique).map({0:(1.,1.,1.),1:the_pal[2]})
#
# HER2_colors_labels = [df_pivot_col_unique[i][3] for i in range(len(df_pivot_col_unique))]
# clinical_data_colors['HER2'] = pd.Series(HER2_colors_labels,index=df_pivot_col_unique).map({0:(1.,1.,1.),1:the_pal[3]})
#
# PR_colors_labels = [df_pivot_col_unique[i][4] for i in range(len(df_pivot_col_unique))]
# clinical_data_colors['PR'] = pd.Series(PR_colors_labels,index=df_pivot_col_unique).map({0:(1.,1.,1.),1:the_pal[6]})
#
# light_pal = sns.light_palette((210, 90, 60), input="husl")
# grade_pal = light_pal[1:len(light_pal):2]
# grade_labels = [df_pivot_col_unique[i][1] for i in range(len(df_pivot_col_unique))]
# unique_grade_labels = list(set(grade_labels))
# grade_lut = dict(zip(unique_grade_labels, grade_pal))
# clinical_data_colors['Grade'] = pd.Series(grade_labels,index=df_pivot_col_unique).map(grade_lut)
#
# colors = ['windows blue', 'amber', 'greyish','brick']
# histology_pal = sns.xkcd_palette(colors)
# histology_labels = [df_pivot_col_unique[i][0] for i in range(len(df_pivot_col_unique))]
# unique_histology_labels = list(set(histology_labels))
# histology_lut = dict(zip(unique_histology_labels, histology_pal))
# clinical_data_colors['Histology'] = pd.Series(histology_labels,index=df_pivot_col_unique).map(histology_lut)

# # Create a custom colormap for the heatmap values
# cmap = sns.diverging_palette(h_neg=10, h_pos=220, s=80, n=7, as_cmap=True)

# # plot clustermap via seaborn (hierarchical clustering, unsupervised clustering)
# sns.set()
# sns.set(font="monospace")
# # cluster_method = ['average','complete','weighted']  # try 'average, 'complete' (max of distance), 'weighted'; method='centroid','ward', 'median' didn't work
# # cluster_metric = ['euclidean','minkowski','cityblock','seuclidean','sqeuclidean','cosine','correlation','chebyshev','canberra','braycurtis']
#
# cluster_method = ['average']
# cluster_metric = ['euclidean']
#
# for mt in cluster_method:
#     for mm in cluster_metric:
#         try:
#             g = sns.clustermap(df_pivot, method=mt, metric=mm,col_colors=clinical_data_colors,row_colors=df_img_feature_colors,
#                                z_score=0, cmap=cmap, linewidths=0, xticklabels=False, yticklabels=False,
#                                figsize=(15, 15))
#             # display the legend for the image feature colormap by adding a empty bar plot but display the legend
#             for label in list(set(img_feature_labels)):
#                 g.ax_row_dendrogram.bar(0, 0, color=img_feature_lut[label], label=label, linewidth=0)
#             g.ax_row_dendrogram.legend(bbox_to_anchor=(1.4, 1.2), loc='best', title='Image Features')
#
#             for label in list(set(histology_labels)):
#                 g.ax_col_dendrogram.bar(0, 0, color=histology_lut[label], label=label, linewidth=0)
#             g.ax_col_dendrogram.legend(bbox_to_anchor=(0.95, 1.62), loc='best', title='Tumor histology')
#
#             for label in list(set(grade_labels)):
#                 g.ax_heatmap.bar(0, 0, color=grade_lut[label], label=label, linewidth=0)
#             g.ax_heatmap.legend(bbox_to_anchor=(0.5, 1.34), loc='best', title='Tumor grade')
#
#             # position the heatmap colorbar appropriately
#             g.cax.set_position([0.92, .2, .03, .45])
#             g.ax_heatmap.set(xlabel='', ylabel='')
#             # print df_pivot.index
#             # print df_pivot.columns
#             # print g.dendrogram_col.linkage
#             # print g.dendrogram_row.linkage
#             # print g.dendrogram_row.linkage.size
#
#             # g.fig.suptitle('cluster method = {}, metric = {}'.format(mt, mm), fontweight='bold')
#
#             # create the file-save directories
#             fw_dir = '{}/her2_Analysis/TP{}_GLCMBIN{}'.format(rootdir,the_tp,the_glcm_Nbin)
#             try:
#                 os.makedirs(fw_dir)
#             except OSError as exception:
#                 if not os.path.isdir(fw_dir):
#                     raise
#             clusterplotfname = '{}/clustermap_alldata_{}_{}.pdf'.format(fw_dir,mt,mm)
#             g.savefig(clusterplotfname)
#
#             plt.show()
#
#         except:
#             print 'error in cluster method: {}, metric: {}'.format(mt,mm)
#             continue
