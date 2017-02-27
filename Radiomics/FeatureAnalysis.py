# -*- coding: utf-8 -*-
"""
Created on 7/25/16 1:34 PM

@author: shuang Shih-ying Huang
@goal: initial effort for the image features vs clinical data analysis

"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import re
import numpy as np

rootdir = '/data/francgrp1/clare_work/Data'
feature_data_dir = '{}/her2_ImageFeatures'.format(rootdir)
pt_data_dir = '{}/her2_ClinicalData'.format(rootdir)

# # MRI mrn vs tumor clinical data
# df_clinical = pd.read_excel('{}/TumorMRIDCEvsPETSUV_Marjan.xlsx'.format(pt_data_dir),sheetname='Data n=126',converters={'PRIMARY_ID':str,'MRN':str,'CR_AccessionSeq':str})

# MRI_PET_data.json inclues all the clinical data and MRI volumn, SER, and PE and PET SUV, MRI and PET acc # from Marjan's study
mri_pet_jname = '{}/MRI_PET_data.json'.format(pt_data_dir)
df_mri_pet = pd.read_json(mri_pet_jname,dtype={'PRIMARY_ID':str,'pt_id':str,'Anon_Accession_num_PET':str,'Accession_num_PET':str,'anon_acc_num':str,'exam_dir_num':str,'pt_id':str,'pt_mrn':str,'MRI_pt_accession_num':str,'MRN_MRI':str,'MRN_PET':str})

def get_pet_breast_side(row):
    side = row['Side']
    lesion_loc = row['Lesion locations']
    if not pd.isnull(side) and 'l' in side.lower():
        return 'L'
    elif not pd.isnull(side) and 'r' in side.lower():
        return 'R'
    elif not pd.isnull(side) and 'b' in side.lower():
        if lesion_loc and 'right' in lesion_loc.lower():
            return 'R'
        elif lesion_loc and 'left' in lesion_loc.lower():
            return 'L'
    else:
        return 'NaN'

# print df_mri_pet.loc[df_mri_pet['pt_id'] == '109',['Side','pt_id','Lesion locations']]
df_mri_pet['PET_breast_side'] = df_mri_pet.apply(get_pet_breast_side,axis=1)
print df_mri_pet.columns
oh = df_mri_pet[df_mri_pet['pt_id'] == '90'].to_dict()
for kk,val in oh.items():
    print kk,val


# fill in the clinical data like ER, HER2 HR status with nan for '8888' or '9999'
df_mri_pet['Sjoerd_ER'] = df_mri_pet['Sjoerd_ER'].replace([8888,9999],np.nan)
df_mri_pet['Sjoerd_HER2'] = df_mri_pet['Sjoerd_HER2'].replace([8888,9999],np.nan)
df_mri_pet['Sjoerd_PR'] = df_mri_pet['Sjoerd_PR'].replace([8888,9999],np.nan)
df_mri_pet['Sjoerd_Grade'] = df_mri_pet['Sjoerd_Grade'].replace([8888,9999],np.nan)


# read all the Image Feature data
all_jsons = glob.glob('{}/*.json'.format(feature_data_dir))
df_img_features = pd.DataFrame()
for fj in all_jsons:
    df = pd.read_json(fj,dtype={'pt_id':str,'pt_mrn':str}) #make sure pt_id and MRN are read in as string
    df['dce_MRI_time_point'] = df['dce_series_dmi_fname'].apply(lambda x: re.search(r'(.+)tp(\d+).dmi',x).group(2))
    df['MRI_breast_side'] = df['dce_series_dmi_fname'].apply(lambda x: re.search(r'(.+)sag(\w+)_tp(\d+).dmi',x).group(2))
    df_img_features = df_img_features.append(df,ignore_index=True)

print df_img_features.ix[df_img_features['pt_id'] == '90',['dce_series_dmi_fname','MRI_breast_side','dce_MRI_time_point','glcm_Nbin','texture_autocorrelation']]


# merge the feature vs. all MRI/PET pt info
df_img_clinical_data = pd.merge(df_mri_pet,df_img_features,how='right',left_on=['pt_mrn','PET_breast_side'],right_on=['pt_mrn','MRI_breast_side'])
# print df_img_clinical_data.ix[1330,['pt_mrn','PET_breast_side','MRI_breast_side','pt_id_x','pt_accession_num','texture_autocorrelation']]
# print df_img_clinical_data.loc[(df_img_clinical_data['pt_id_x'] == '89') & (df_img_features['dce_MRI_time_point'] == '1') & (df_img_features['glcm_Nbin'] == 128),:]

# # print df_img_clinical_data[(df_img_clinical_data['pt_id_x'] == '109') & (df_img_clinical_data['dce_MRI_time_point'] == '1') & (df_img_clinical_data['glcm_Nbin'] == 128)]
#
# select the data of interest
# clustermap for each time point feature data
the_tp = '1'
the_glcm_Nbin =128

df_tp_glcmNbin = df_img_clinical_data[(df_img_clinical_data['dce_MRI_time_point'] == the_tp) & (df_img_clinical_data['glcm_Nbin'] == the_glcm_Nbin)]
df_tp_glcmNbin['MRN_side'] = df_tp_glcmNbin.apply(lambda row: '{}_{}'.format(row['pt_mrn'],row['MRI_breast_side']),axis=1)
print df_tp_glcmNbin.columns
# print df_tp_glcmNbin.loc[df_tp_glcmNbin.duplicated('MRN_side',keep=False)]

oh1 = df_tp_glcmNbin.ix[770,['pt_id_x','dce_series_dmi_fname','MRI_dce_fname','MRI_breast_side','PET_breast_side','texture_autocorrelation']].tolist()
oh2 = df_tp_glcmNbin.ix[771,['pt_id_x','dce_series_dmi_fname','MRI_dce_fname','MRI_breast_side','PET_breast_side','texture_autocorrelation']].tolist()

for ii in range(len(oh1)):
    print oh1[ii]
    print oh2[ii]


# print len(df_tp_glcmNbin)
# print len(df_tp_glcmNbin['MRN_side'].unique())

# # print df_tp_glcmNbin.ix[:,['pt_mrn','MRI_breast_side','PET_breast_side','Accession_num_PET','pt_accession_num','pt_id_x','PRIMARY_ID','MRI_pt_accession_num','MRN_MRI','MRN_PET','Lesion locations','texture_idmn']]
#
# # #TESTING
# # test_mrn = ['47712882','11854133','47454751','32140564','46540182','08075719','42704318','49520689']
# # for tt in test_mrn:
# #     print df_mri_pet.loc[df_mri_pet['pt_mrn'] == tt,['pt_mrn','PET_breast_side','Accession_num_PET','Side','Lesion locations']]
# #
# #     print df_img_features.loc[df_img_features['pt_mrn'] == tt,['pt_id','pt_mrn','MRI_breast_side','texture_idmn','glcm_Nbin','dce_MRI_time_point']]
# #
# #     print df_img_clinical_data.loc[df_img_clinical_data['pt_mrn'] == tt,['pt_mrn','MRI_breast_side','PET_breast_side','Accession_num_PET','Side','Lesion locations','texture_idmn','glcm_Nbin','dce_MRI_time_point']]
# #
# #     print df_tp_glcmNbin.loc[df_tp_glcmNbin['pt_mrn'] == tt,['pt_mrn','MRI_breast_side','PET_breast_side','Accession_num_PET','pt_accession_num','pt_id_x','PRIMARY_ID','MRI_pt_accession_num','MRN_MRI','MRN_PET','Lesion locations']]
#
# # print df_tp_glcmNbin.loc[:,['pt_mrn','MRI_breast_side','PET_breast_side','Accession_num_PET','pt_accession_num','pt_id_x','PRIMARY_ID','MRI_pt_accession_num','MRN_MRI','MRN_PET','Lesion locations']]
#
# drop un-needed columns
df_tp_glcmNbin = df_tp_glcmNbin.drop(['organ_mask','process_name','process_version','Last name','BP SUV ave','Corrected SUV','Side','MRI mass size','MRI_CR_AccessionSeq','Lesion locations','MRI_dce_fname','MRN_MRI','MRN_PET','texture_glcm_offset','anon_acc_num','exam_dir','exam_dir_num','dce_MRI_time_point','dce_series_dmi_fname','glcm_Nbin','img_norm_Nbin','pt_id_y'],axis=1)

# data munging: get the average of all the features values (over all offsets) and split minmax into separate variables
texture_cols = ['texture_autocorrelation','texture_cluster_prominence','texture_cluster_shade','texture_cluster_tendency','texture_contrast','texture_correlation','texture_diff_entropy','texture_dissimilarity',
                'texture_energy','texture_entropy','texture_homogeneity1','texture_homogeneity2','texture_idmn','texture_idn','texture_inv_var','texture_maxprob','texture_sum_avg','texture_sum_entropy','texture_sum_var']
for tc in texture_cols:
    df_tp_glcmNbin[tc +'_avg'] = df_tp_glcmNbin[tc].apply(np.mean)
    df_tp_glcmNbin = df_tp_glcmNbin.drop(tc,axis=1)

df_tp_glcmNbin['FOstats_min'] = df_tp_glcmNbin['FOstats_minmax'].apply(lambda x: x[0])
df_tp_glcmNbin['FOstats_max'] = df_tp_glcmNbin['FOstats_minmax'].apply(lambda x: x[1])
df_tp_glcmNbin = df_tp_glcmNbin.drop('FOstats_minmax',axis=1)

# re-format the radiomic features into columnes of feature name and feature values
pt_info_names = ['Accession_num_PET','Anon_Accession_num_PET','PRIMARY_ID','pt_id_x','pt_mrn','pt_accession_num','MRI_breast_side']
clinical_feature_names = ['Marjan_Histology','Marjan_Size (mm)','Sjoerd_ER','Sjoerd_HER2','Sjoerd_PR','Sjoerd_Grade']
img_feature_names = ['MRI_PE2ROI_PE2_PEAK','MRI_PE2ROI_SER_MEAN','MRI_SERROI_PE2_PEAK','MRI_SERROI_SER_MEAN','MRI_VOLUME_BLUE','MRI_VOLUME_GREEN','MRI_VOLUME_PURPLE','MRI_VOLUME_RED','MRI_VOLUME_TUMOR','MRI_VOLUME_TUM_BLU','MRI_VOLUME_WHITE','FOstats_energy',
                     'FOstats_entropy','FOstats_kurtosis','FOstats_mean','FOstats_min','FOstats_max','FOstats_skewness','FOstats_uniformity','FOstats_variance','ShapeSize_compactness1','ShapeSize_compactness2','ShapeSize_max_euc_dis',
                     'ShapeSize_spherical_disproportion','ShapeSize_sphericity','ShapeSize_surf_area_cm2','ShapeSize_surface2volratio','ShapeSize_vol_cm3'] + [ss + '_avg' for ss in texture_cols]

df_flatten_data = pd.melt(df_tp_glcmNbin,id_vars=pt_info_names + clinical_feature_names,value_vars=img_feature_names)
df_flatten_data['pt_id_side'] = df_flatten_data.apply(lambda row: '{}_{}'.format(row['pt_id_x'],row['MRI_breast_side']),axis=1)
print len(df_flatten_data['pt_id_side'].unique())
print len(df_flatten_data['pt_id_side'])


print df_flatten_data.columns
print df_flatten_data.index
df_pivot = df_flatten_data.pivot(index='variable',columns='pt_id_side',values='value')
print df_pivot


df_flatten_data = pd.melt(df_tp_glcmNbin,id_vars=pt_info_names + clinical_feature_names,value_vars=img_feature_names)
print df_flatten_data
print df_flatten_data.columns
print df_flatten_data.index

# pivot_table to prepare for clustermap, make sure dropna() before clustermap (or else there will be an error of 'Z' contains negative indices...==> pivot_table aggregate the values with np.mean by default, not what i intended to do!!
df_pivot_table = pd.pivot_table(df_flatten_data,values='value',columns=['Sjoerd_ER','Sjoerd_HER2','Sjoerd_PR','Sjoerd_Grade','Marjan_Histology'],index='variable')


img_feature_vars = df_pivot_table.index.tolist()

# categorical label for the image feature vars
img_feature_labels = []
for ss in img_feature_vars:
    if re.search(r'FOstats_(.+)',ss):
        img_feature_labels.append('First-order Stats')
    elif re.search(r'MRI_(.+)',ss):
        img_feature_labels.append('Marjan MRI features')
    elif re.search(r'ShapeSize_(.+)',ss):
        img_feature_labels.append('Shape and Size')
    elif re.search(r'texture_(.+)',ss):
        img_feature_labels.append('GLCM texture features')

# img_feature_pal = sns.cubehelix_palette(len(set(img_feature_labels)),
#                                             light=.9, dark=.1, reverse=True,
#                                             start=1, rot=-2)
# img_feature_pal = sns.cubehelix_palette(len(set(img_feature_labels)))

# img_feature_pal = sns.color_palette('PuBuGn_d', len(set(img_feature_labels)))
img_feature_pal = sns.light_palette('navy',n_colors=len(set(img_feature_labels)),reverse=True)

img_feature_lut = dict(zip(map(str, list(set(img_feature_labels))), img_feature_pal))
img_feature_colors = pd.Series(img_feature_labels,index=df_pivot_table.index).map(img_feature_lut)
df_img_feature_colors = pd.DataFrame(dict(ImageFeature=img_feature_colors))
df_img_feature_colors.columns = ['']
# print df_img_feature_colors.columns


# categorize the clinical data into sub-group for visualization (set col_colors for the clustermap)
the_pal = sns.color_palette('Set2',10)
ER_colors = pd.Series(df_pivot_table.columns.get_level_values('Sjoerd_ER'),index=df_pivot_table.columns).map({0:(1.,1.,1.),1:the_pal[0]})
HER2_colors = pd.Series(df_pivot_table.columns.get_level_values('Sjoerd_HER2'),index=df_pivot_table.columns).map({0:(1.,1.,1.),1:the_pal[1]})
PR_colors = pd.Series(df_pivot_table.columns.get_level_values('Sjoerd_PR'),index=df_pivot_table.columns).map({0:(1.,1.,1.),1:the_pal[3]})

light_pal = sns.light_palette((210, 90, 60), input="husl")
grade_pal = light_pal[1:len(light_pal):2]
grade_labels = df_pivot_table.columns.get_level_values('Sjoerd_Grade')
grade_lut = dict(zip(np.sort(grade_labels.unique()), grade_pal))
Grade_colors = pd.Series(grade_labels,index=df_pivot_table.columns).map(grade_lut)


colors = ['windows blue', 'amber', 'greyish','brick']
histology_pal = sns.xkcd_palette(colors)
histology_labels =  df_pivot_table.columns.get_level_values('Marjan_Histology')
histology_lut = dict(zip(map(str, histology_labels.unique()), histology_pal))
histology_colors = pd.Series(histology_labels,index=df_pivot_table.columns).map(histology_lut)

clinical_data_colors = pd.DataFrame(dict(ER=ER_colors,HER2=HER2_colors,PR=PR_colors,Grade=Grade_colors,Histology=histology_colors))

# # Create a custom colormap for the heatmap values
cmap = sns.diverging_palette(h_neg=10, h_pos=220, s=80, n=7, as_cmap=True)

# plot clustermap via seaborn (hierarchical cluster)
sns.set()
sns.set(font="monospace")
cluster_method = 'weighted'  # try 'average, 'complete' (max of distance), 'weighted'; method='centroid','ward', 'median' didn't work
cluster_metric = ['euclidean','minkowski','cityblock','seuclidean','sqeuclidean','cosine','correlation','chebyshev','canberra','braycurtis']

for mm in cluster_metric:
    try:
        g = sns.clustermap(df_pivot_table,method=cluster_method,metric=mm,row_colors=df_img_feature_colors,col_colors=clinical_data_colors,
                           z_score=0,cmap=cmap,linewidths=0, xticklabels=False, yticklabels=False,
                           figsize=(15,15))

        # display the legend for the image feature colormap by adding a empty bar plot but display the legend
        for label in list(set(img_feature_labels)):
            g.ax_row_dendrogram.bar(0, 0, color=img_feature_lut[label],label=label, linewidth=0)
        g.ax_row_dendrogram.legend(bbox_to_anchor=(1.4,1.2),loc='best',title='Image Features')

        for label in list(set(histology_labels)):
            g.ax_col_dendrogram.bar(0, 0, color=histology_lut[label],label=label, linewidth=0)
        g.ax_col_dendrogram.legend(bbox_to_anchor=(0.95,1.62),loc='best',title='Tumor histology')

        for label in list(set(grade_labels)):
            g.ax_heatmap.bar(0, 0, color=grade_lut[label],label=label, linewidth=0)
        g.ax_heatmap.legend(bbox_to_anchor=(0.5,1.34),loc='best',title='Tumor grade')


        # position the heatmap colorbar appropriately
        g.cax.set_position([0.92, .2, .03, .45])
        g.ax_heatmap.set(xlabel='',ylabel='')

        g.fig.suptitle('cluster method = {}, metric = {}'.format(cluster_method,mm),fontweight='bold')

        clusterplotfname = '{}/her2_Analysis/clustermap_{}_{}.pdf'.format(rootdir,cluster_method,mm)
        g.savefig(clusterplotfname)

        plt.show()
    except:
        print 'error in cluster method: {}, metric: {}'.format(cluster_method,mm)
        continue








