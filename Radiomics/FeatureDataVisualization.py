# -*- coding: utf-8 -*-
"""
Created on 11/15/16 11:13 AM

@author: shuang Shih-ying Huang
@goal: visualize the image feature data to understand the data distribution

"""
import glob
import pandas as pd
import re
from itertools import product
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


rootdir = '/data/francgrp1/clare_work/Data'
feature_data_dir = '{}/her2_ImageFeatures'.format(rootdir)

the_pt_id = '16'
the_feature = ['texture_cluster_prominence','texture_cluster_shade','texture_cluster_tendency','texture_contrast','texture_correlation',
               'texture_diff_entropy','texture_dissimilarity','texture_energy','texture_entropy','texture_homogeneity1',
               'texture_homogeneity2','texture_idmn','texture_idn','texture_inv_var','texture_maxprob','texture_sum_avg','texture_sum_entropy',
               'texture_sum_var']


# # visualize MRI feature data
# all_jsons = glob.glob('{}/IsoVoxelSize/MRI*.json'.format(feature_data_dir))
# df_img_features = pd.DataFrame()
# for fj in all_jsons:
#     df = pd.read_json(fj,dtype={'pt_id':str,'pt_mrn':str}) #make sure pt_id and MRN are read in as string
#     df['dce_MRI_time_point'] = df['dce_series_dmi_fname'].apply(lambda x: re.search(r'(.+)tp(\d+).dmi',x).group(2))
#     df['Laterality'] = df['dce_series_dmi_fname'].apply(lambda x: re.search(r'(.+)sag(\w+)_tp(\d+).dmi',x).group(2))
#     df_img_features = df_img_features.append(df,ignore_index=True)
#
# for ff in the_feature:
#     print ff
#     df_test = df_img_features.ix[df_img_features['pt_id'] == the_pt_id,['dce_MRI_time_point','glcm_Nbin','img_norm_Nbin',ff]]
#
#     # reformat the dataframe to flatten the data in the numpy array into a dataframe entry
#     tp_list = df_test['dce_MRI_time_point'].unique()
#     glcmNbin_list = df_test['glcm_Nbin'].unique()
#     combo_list = list(product(tp_list,glcmNbin_list))
#     df_all = pd.DataFrame()
#     for ii in range(len(combo_list)):
#         tp = combo_list[ii][0]
#         gbin = combo_list[ii][1]
#         tmp_list = df_test.ix[(df_test['dce_MRI_time_point'] == tp) & (df_test['glcm_Nbin'] == gbin),ff].values[0]
#         if ii == 0:
#             col_dir = ['dir_{}'.format(dd) for dd in range(len(tmp_list))]
#             data_all = np.array([['','mri_time_point','glcm_Nbin'] + col_dir])
#         data_all = np.append(data_all,[[ii,tp,gbin] + tmp_list],axis=0)
#     df_all = pd.DataFrame(data=data_all[1:,1:],index=data_all[1:,0],columns=data_all[0,1:])
#     df_all = df_all.convert_objects(convert_numeric=True)
#
#     # flatten the dataframe
#     df2 = pd.melt(df_all, id_vars=['mri_time_point', 'glcm_Nbin'], value_vars=col_dir, value_name=ff)
#
#     # Initialize a grid of plots with an Axes for each walk
#     grid = sns.FacetGrid(df2, col='mri_time_point', row='glcm_Nbin', col_order=[1, 2, 3], row_order=[64, 128, 256])
#
#     # Draw a line plot to show the trajectory of each random walk
#     grid.map(sns.pointplot, 'variable', ff, scale=0.7)
#     grid.set(xticks=np.arange(13))
#     plt.show()

# # PET image features
# all_jsons = glob.glob('{}/IsoVoxelSize/PET*.json'.format(feature_data_dir))
# df_img_features = pd.DataFrame()
# for fj in all_jsons:
#     df = pd.read_json(fj,dtype={'pt_id':str,'pt_mrn':str,'pt_accession_num':str}) #make sure pt_id and MRN are read in as string
#     df_img_features = df_img_features.append(df,ignore_index=True)
#
# # visualize PET feature data
# for ff in the_feature:
#     print ff
#     df_test = df_img_features.ix[df_img_features['pt_id'] == the_pt_id, ['glcm_Nbin', 'img_norm_Nbin', ff]]
#
#     # reformat the dataframe to flatten the data in the numpy array into a dataframe entry
#     imgNbin_list = df_test['img_norm_Nbin'].unique()
#     glcmNbin_list = df_test['glcm_Nbin'].unique()
#     print imgNbin_list
#     print glcmNbin_list
#     combo_list = list(product(imgNbin_list, glcmNbin_list))
#     df_all = pd.DataFrame()
#     for ii in range(len(combo_list)):
#         ibin = combo_list[ii][0]
#         gbin = combo_list[ii][1]
#         tmp_list = df_test.ix[(df_test['img_norm_Nbin'] == ibin) & (df_test['glcm_Nbin'] == gbin), ff].values[0]
#         if ii == 0:
#             col_dir = ['dir_{}'.format(dd) for dd in range(len(tmp_list))]
#             data_all = np.array([['', 'img_norm_Nbin', 'glcm_Nbin'] + col_dir])
#         data_all = np.append(data_all, [[ii, ibin, gbin] + tmp_list], axis=0)
#     df_all = pd.DataFrame(data=data_all[1:, 1:], index=data_all[1:, 0], columns=data_all[0, 1:])
#     df_all = df_all.convert_objects(convert_numeric=True)
#
#     # flatten the dataframe
#     df2 = pd.melt(df_all,id_vars=['img_norm_Nbin','glcm_Nbin'],value_vars=col_dir,value_name=ff)
#
#     # Initialize a grid of plots with an Axes for each walk
#     grid = sns.FacetGrid(df2, col='img_norm_Nbin',row='glcm_Nbin',col_order=imgNbin_list,row_order=glcmNbin_list)
#
#     # Draw a line plot to show the trajectory of each random walk
#     grid.map(sns.pointplot,'variable',ff, scale=0.7)
#     # grid.set(xticks=np.arange(13))
#     plt.show()


# # visualize all texture features across all tumors
# # MRI feature data
# all_jsons = glob.glob('{}/IsoVoxelSize/MRI*.json'.format(feature_data_dir))
# df_img_features = pd.DataFrame()
# for fj in all_jsons:
#     df = pd.read_json(fj,dtype={'pt_id':str,'pt_mrn':str}) #make sure pt_id and MRN are read in as string
#     df['dce_MRI_time_point'] = df['dce_series_dmi_fname'].apply(lambda x: re.search(r'(.+)tp(\d+).dmi',x).group(2))
#     df['Laterality'] = df['dce_series_dmi_fname'].apply(lambda x: re.search(r'(.+)sag(\w+)_tp(\d+).dmi',x).group(2))
#     df_img_features = df_img_features.append(df,ignore_index=True)
#
# for tc in the_feature:
#     df_img_features[tc + '_avg'] = df_img_features[tc].apply(np.mean)
#     df_img_features = df_img_features.drop(tc, axis=1)
#
# col_name = df_img_features['dce_MRI_time_point'].unique()
# row_name = df_img_features['glcm_Nbin'].unique()
#
# the_avg_feature = [ss + '_avg' for ss in the_feature]
# save_fdir = '/data/francgrp1/clare_work/Data/her2_Analysis/DataViz'
# for ff in the_avg_feature:
#     print ff
#
#     # Initialize a grid of plots with an Axes for each walk
#     grid = sns.FacetGrid(df_img_features,col='dce_MRI_time_point',row='glcm_Nbin',col_order=col_name,row_order=row_name)
#
#     # Draw a line plot to show the trajectory of each random walk
#     grid.map(plt.hist,ff)
#
#     # find min and max
#     # grid.set(ylim=(0.5*df_img_features[ff].min(),1.5*df_img_features[ff].max()))
#
#     fig_fname = '{}/MRI_facetgrid_hist_{}.pdf'.format(save_fdir, ff)
#     plt.savefig(fig_fname)
#
#     plt.show()

# PET feature data
all_jsons = glob.glob('{}/IsoVoxelSize/PET*.json'.format(feature_data_dir))
df_img_features = pd.DataFrame()
for fj in all_jsons:
    df = pd.read_json(fj,dtype={'pt_id':str,'pt_mrn':str,'pt_accession_num':str}) #make sure pt_id and MRN are read in as string
    df_img_features = df_img_features.append(df,ignore_index=True)

for tc in the_feature:
    df_img_features[tc] = df_img_features[tc].apply(lambda x: x if all(x) else np.nan)
    if tc == 'texture_correlation':
        print df_img_features[tc]
    df_img_features[tc + '_avg'] = df_img_features[tc].apply(np.nanmean)
    df_img_features = df_img_features.drop(tc, axis=1)

the_avg_feature = [ss + '_avg' for ss in the_feature]
save_fdir = '/data/francgrp1/clare_work/Data/her2_Analysis/DataViz'
imgNbin_list = df_img_features['img_norm_Nbin'].unique()
glcmNbin_list = df_img_features['glcm_Nbin'].unique()
for ff in the_avg_feature:
    print ff

    # Initialize a grid of plots with an Axes for each walk
    grid = sns.FacetGrid(df_img_features,col='img_norm_Nbin',row='glcm_Nbin',col_order=imgNbin_list,row_order=glcmNbin_list)

    # Draw a line plot to show the trajectory of each random walk
    grid.map(plt.hist,ff)

    # find min and max
    # grid.set(ylim=(0.5*df_img_features[ff].min(),1.5*df_img_features[ff].max()))

    fig_fname = '{}/PET_facetgrid_hist_{}.pdf'.format(save_fdir, ff)
    plt.savefig(fig_fname)

    plt.show()



