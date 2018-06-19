# -*- coding: utf-8 -*-
"""
Created on 3/7/17 1:05 PM

@author: shuang Shih-ying Huang
@goal: organize the image, mask, and other data for her2 dataset (MRI, PET, and patient outcome if there's any)

"""

import pandas as pd
import re
import numpy as np
import os
import shutil,errno
import glob

def makedirs(file_path):
    if os.path.exists(file_path):
        if not os.path.isdir(file_path):
            raise IOError('O_O: {} is not a directory!'.format(file_path))
    else:
        os.makedirs(file_path)

def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc: # python >2.5
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise



rootdir = '/data/francgrp1/clare_work/Data'
feature_data_dir = '{}/her2_ImageFeatures/IsoVoxelSize'.format(rootdir)
pt_data_dir = '{}/her2_ClinicalData'.format(rootdir)

# MRI_PET_ImgFeature json
fname = '{}/MRI_PET_imgfeat_data.json'.format(pt_data_dir)
df_mri_pet_imgfeat = pd.read_json(fname,dtype={'pt_mrn':str,'PRIMARY_ID':str,'MRN':str,'Anon_Accession_Num':str})
df_mri_pet_imgfeat['MRN'] = df_mri_pet_imgfeat['MRN'].apply(lambda x: '{0:0>8}'.format(x))

# fill in the clinical data like ER, HER2 HR status with nan for '8888' or '9999'
df_mri_pet_imgfeat['Sjoerd_ER'] = df_mri_pet_imgfeat['Sjoerd_ER'].replace([8888,9999],np.nan)
df_mri_pet_imgfeat['Sjoerd_HER2'] = df_mri_pet_imgfeat['Sjoerd_HER2'].replace([8888,9999],np.nan)
df_mri_pet_imgfeat['Sjoerd_PR'] = df_mri_pet_imgfeat['Sjoerd_PR'].replace([8888,9999],np.nan)
df_mri_pet_imgfeat['Sjoerd_Grade'] = df_mri_pet_imgfeat['Sjoerd_Grade'].replace([8888,9999],np.nan)
col_name_lut = {'Sjoerd_Grade':'Tumor_Grade', 'Marjan_Histology': 'Tumor_Histology',
                                   'Sjoerd_ER': 'ER_status', 'Sjoerd_HER2': 'HER2_status','Sjoerd_PR':'PR_status'}
df_mri_pet_imgfeat.rename(columns=col_name_lut,inplace=True)

# pick the needed columns for cluster analysis
clinical_colname = ['Anon_Accession_Num','Laterality','PRIMARY_ID','MRN','CR_AccessionSeq','MRI_VOLUME_TUM_BLU','MRI_VOLUME_TUMOR',
                    'MRI_VOLUME_WHITE','MRI_VOLUME_RED','MRI_VOLUME_GREEN','MRI_VOLUME_PURPLE','MRI_VOLUME_BLUE','MRI_SERROI_SER_MEAN',
                    'MRI_PE2ROI_PE2_PEAK','MRI_SERROI_PE2_PEAK','MRI_PE2ROI_SER_MEAN','SUV max'] + col_name_lut.values()

df_clinicaldata = df_mri_pet_imgfeat.ix[:,clinical_colname]

img_feature_names = ['MRI_PE2ROI_PE2_PEAK','MRI_PE2ROI_SER_MEAN','MRI_SERROI_PE2_PEAK','MRI_SERROI_SER_MEAN','MRI_VOLUME_BLUE',
                     'MRI_VOLUME_GREEN','MRI_VOLUME_PURPLE','MRI_VOLUME_RED','MRI_VOLUME_TUMOR','MRI_VOLUME_TUM_BLU','MRI_VOLUME_WHITE',
                     'FOstats_energy','FOstats_entropy','FOstats_kurtosis','FOstats_mean','FOstats_min','FOstats_max','FOstats_skewness',
                     'FOstats_uniformity','FOstats_variance','ShapeSize_compactness1','ShapeSize_compactness2','ShapeSize_max_euc_dis',
                     'ShapeSize_spherical_disproportion','ShapeSize_sphericity','ShapeSize_surf_area_cm2','ShapeSize_surface2volratio',
                     'ShapeSize_vol_cm3','texture_autocorrelation','texture_cluster_prominence','texture_cluster_shade','texture_cluster_tendency',
                     'texture_contrast','texture_correlation','texture_diff_entropy','texture_dissimilarity','texture_energy','texture_entropy',
                     'texture_homogeneity1','texture_homogeneity2','texture_idmn','texture_idn','texture_inv_var','texture_maxprob',
                     'texture_sum_avg','texture_sum_entropy','texture_sum_var']

# read all the PET Image Feature data
all_jsons_pet = glob.glob('{}/PET*.json'.format(feature_data_dir))
df_img_features_pet = pd.DataFrame()
for fj in all_jsons_pet:
    df = pd.read_json(fj,dtype={'pt_id':str,'pt_mrn':str,'pt_accession_num':str}) #make sure pt_id and MRN are read in as string
    df_img_features_pet = df_img_features_pet.append(df,ignore_index=True)

img_feats = ['FOstats_minmax','FOstats_energy','FOstats_entropy','FOstats_kurtosis','FOstats_mean','FOstats_skewness',
                 'FOstats_uniformity','FOstats_variance','ShapeSize_compactness1','ShapeSize_compactness2','ShapeSize_max_euc_dis',
                 'ShapeSize_spherical_disproportion','ShapeSize_sphericity','ShapeSize_surf_area_cm2','ShapeSize_surface2volratio',
                 'ShapeSize_vol_cm3','texture_autocorrelation','texture_cluster_prominence','texture_cluster_shade','texture_cluster_tendency',
                 'texture_contrast','texture_correlation','texture_diff_entropy','texture_dissimilarity','texture_energy','texture_entropy',
                 'texture_homogeneity1','texture_homogeneity2','texture_idmn','texture_idn','texture_inv_var','texture_maxprob',
                 'texture_sum_avg','texture_sum_entropy','texture_sum_var']

df_img_features_pet = df_img_features_pet.drop(img_feats,axis=1)

all_jsons_mri = glob.glob('{}/MRI*.json'.format(feature_data_dir))
df_img_features_mri = pd.DataFrame()
for fj in all_jsons_mri:
    df = pd.read_json(fj,dtype={'pt_id':str,'pt_mrn':str}) #make sure pt_id and MRN are read in as string
    df['dce_MRI_time_point'] = df['dce_series_dmi_fname'].apply(lambda x: re.search(r'(.+)tp(\d+).dmi',x).group(2))
    df['Laterality'] = df['dce_series_dmi_fname'].apply(lambda x: re.search(r'(.+)sag(\w+)_tp(\d+).dmi',x).group(2))
    df_img_features_mri = df_img_features_mri.append(df,ignore_index=True)

df_img_features_mri = df_img_features_mri.drop(img_feats,axis=1)

the_tp = '2'
the_glcm_Nbin = 128
df_mri = df_img_features_mri[(df_img_features_mri['dce_MRI_time_point'] == the_tp) & (df_img_features_mri['glcm_Nbin'] == the_glcm_Nbin)]

the_img_norm_Nbin = 128
the_glcm_Nbin_pet = 64
df_pet = df_img_features_pet[(df_img_features_pet['glcm_Nbin'] == the_glcm_Nbin_pet) & (df_img_features_pet['img_norm_Nbin']==the_img_norm_Nbin)]

tmp1 = pd.merge(df_mri,df_clinicaldata,how='left',left_on=['pt_mrn','Laterality'],right_on=['MRN','Laterality'])
tmp2 = pd.merge(tmp1,df_pet,left_on=['Anon_Accession_Num','Laterality'],right_on=['pt_accession_num','breast_side'],suffixes=['_mri','_pet'])
tmp2['ptid_side'] = tmp2.apply(lambda x: '{}_{}'.format(x['pt_id_pet'],x['Laterality']),axis=1)
        # print df_tp_glcmNbin.loc[df_tp_glcmNbin['pt_id'] == '96','ptid_side']
# print tmp2.ix[tmp2['pt_accession_num_pet'].isnull(),['Anon_Accession_Num','MRN','Laterality','pt_id_mri','pt_id_pet']]

# generate outcome variables the dataset
# outcome_fname = '{}/her2_ClinicalData/her2_outcome.csv'.format(rootdir)
outcome_fname = '{}/her2_ClinicalData/her2_outcome_NJupdated.csv'.format(rootdir)
her2_outcome_df = pd.read_csv(outcome_fname,dtype={'MRN': str,'Recurrence Type Summary': str})
her2_outcome_df.rename(columns={'Brith Date':'Birth Date'},inplace=True)
tmp3 = pd.merge(tmp2,her2_outcome_df,how='left',on='MRN')
tmp3['tumor_id'] = np.array(tmp3.index.tolist()) + 1

sens_feat_name = ['Laterality','dce_MRI_time_point','dce_series_dmi_fname','glcm_Nbin_mri','img_norm_Nbin_mri','organ_mask_mri',
                  'process_name_mri','process_version_mri','pt_accession_num_mri','pt_mrn_mri','texture_glcm_offset_mri',
                  'voxel_size_mm3_mri','Anon_Accession_Num','PRIMARY_ID','MRN','CR_AccessionSeq','glcm_Nbin_pet','img_norm_Nbin_pet',
                  'organ_mask_pet','pet_series_fdir','process_name_pet','process_version_pet','pt_accession_num_pet',
                  'pt_mrn_pet','texture_glcm_offset_pet','voxel_size_mm3_pet','Last name','Accession/Seq']

# make sure 'tumor_id' column is the first column instead of the last
tmp4 = tmp3.drop(sens_feat_name,axis=1)
cols = tmp4.columns.tolist()
the_idx_col = 'tumor_id'
cols.pop(cols.index(the_idx_col))
cols.insert(0,the_idx_col)
tmp4 = tmp4[cols]

# pet/mri/mask data re-org
data_dir = '/data/francgrp1/breast_radiomics/her2/Data'
# makedirs(data_dir)

# save the de-identified outcome data
outcome_deid_fname = '{}/her2_outcomes_anon.csv'.format(data_dir)
tmp4.to_csv(outcome_deid_fname, index=False)

# pet_data_dir = '/data/francgrp1/breast_radiomics/her2/ALL_TUMORS'
# pet_mask_dir = '{}/mevis_manual_segmentation'.format(pet_data_dir)
# fragmentedpet_dcm_ids = [16,17,18,22,25,27,30]
#
# mri_data_dir = '/data/francgrp1/breast_radiomics/her2/MRI'
# mri_info_dir = '{}/MORE_INFO'.format(mri_data_dir)
#
# tumor_id = tmp2.index.tolist()
# for ii in tumor_id:
#     # create tumor_id dir
#     tumor_dir = '{}/{:0>3d}'.format(data_dir,ii+1)
#     makedirs(tumor_dir)
#
#     # pull in PET image related data
#     pet_dir = '{}/PET'.format(tumor_dir)
#     pet_series_fdir = str(tmp2.ix[ii,'pet_series_fdir'])
#     pt_id_pet = int(tmp2.ix[ii, 'pt_id_pet'])
#     breast_side = str(tmp2.ix[ii, 'breast_side'])
#     shutil.copytree(pet_series_fdir,pet_dir)
#
#     # copy over corrected PET volume NRRD files for the fragmented PET volumes
#     if pt_id_pet in fragmentedpet_dcm_ids:
#         find_pet_mevis_nrrd = glob.glob('{}/{:0>3d}/*/PET_voxshiftON_rescale.nrrd'.format(pet_data_dir,pt_id_pet))
#         if find_pet_mevis_nrrd:
#             pet_mevis_nrrd = find_pet_mevis_nrrd[0]
#             shutil.copy(pet_mevis_nrrd,pet_dir)
#         else:
#             print 'pt_id,side:{},{}, cannot find PET_voxshiftON_rescale.nrrd'.format(pt_id_pet,breast_side)
#
#     # copy the PET mask over
#     mask_dir = '{}/MASK'.format(tumor_dir)
#     makedirs(mask_dir)
#
#     check_mask_fname = glob.glob('{}/{:0>3d}*sag{}*_uc_origin.nrrd'.format(pet_mask_dir, pt_id_pet, breast_side))
#     if check_mask_fname:
#         mask_fname = check_mask_fname[0]
#     else:
#         print 'pt_id_pet, side = {}, {} Cannot find the mask .nrrd file with the breast side!'.format(pt_id_pet,breast_side)
#     shutil.copy(mask_fname,mask_dir)
#
#     # pull over MRI image related data
#     mri_dir = '{}/MRI'.format(tumor_dir)
#     pt_id_mri = int(tmp2.ix[ii,'pt_id_mri'])
#     mri_dce_dir = glob.glob('{}/{:0>3d}/*/*/dce_anonymized'.format(mri_data_dir, pt_id_mri))[0]
#     shutil.copytree(mri_dce_dir,mri_dir)
#
#     find_info_fname = glob.glob('{}/{}_{}*.info.json'.format(mri_info_dir,pt_id_mri,breast_side))
#     if find_info_fname:
#         info_fname = find_info_fname[0]
#         shutil.copy(info_fname,mask_dir)
#     else:
#         print 'pt_id,side: {},{}, cannot find .info.json file!'.format(pt_id_mri,breast_side)
#
#     # Get the SER mask (.dmi) from the brtool output
#     findserdmi = glob.glob('{}/{}_{}*.SER_MASK.dmi'.format(mri_info_dir, pt_id_mri, breast_side))
#     if findserdmi:
#         ser_mask_fname = findserdmi[0]
#     else:
#         print 'pt_id,side={},{}: SER mask DMI file was NOT found!'.format(pt_id_mri,breast_side)
#     shutil.copy(ser_mask_fname,mask_dir)








