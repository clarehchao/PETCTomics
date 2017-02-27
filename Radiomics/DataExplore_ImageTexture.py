# -*- coding: utf-8 -*-
"""
Created on 9/22/16 9:11 AM

@author: shuang Shih-ying Huang
@goal: explore the image texture from the image texture extracted from the MRI images

"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import re
import os
import ITKImageHelper
import ImageProcess as ImSeg
import dmi
import image_geometry
import numpy as np


# take a look at the distribution of the texture features for different settings
rootdir = '/data/francgrp1/clare_work/Data'
feature_data_dir = '{}/her2_ImageFeatures'.format(rootdir)

# read all the Image Feature data
all_jsons = glob.glob('{}/*.json'.format(feature_data_dir))
df_img_features = pd.DataFrame()
for fj in all_jsons:
    df = pd.read_json(fj,dtype={'pt_id':str,'pt_mrn':str}) #make sure pt_id and MRN are read in as string
    df['dce_MRI_time_point'] = df['dce_series_dmi_fname'].apply(lambda x: re.search(r'(.+)tp(\d+).dmi',x).group(2))
    df['Laterality'] = df['dce_series_dmi_fname'].apply(lambda x: re.search(r'(.+)sag(\w+)_tp(\d+).dmi',x).group(2))
    df_img_features = df_img_features.append(df,ignore_index=True)

df_img_features['homogeneity1_offset1'] = df_img_features['texture_homogeneity1'].apply(lambda x: x[0])
df_img_features['homogeneity1_offset5'] = df_img_features['texture_homogeneity1'].apply(lambda x: x[5])
df_img_features['correlation_offset1'] = df_img_features['texture_correlation'].apply(lambda x:x[3])
df_img_features['dissimilarity_offset1'] = df_img_features['texture_dissimilarity'].apply(lambda x:x[0])

g1 = sns.FacetGrid(df_img_features,row='dce_MRI_time_point',col='glcm_Nbin',margin_titles=True)
g1.map(plt.hist,'dissimilarity_offset1',color='steelblue')
plt.show()

g2 = sns.FacetGrid(df_img_features,row='dce_MRI_time_point',col='glcm_Nbin',margin_titles=True)
g2.map(plt.hist,'homogeneity1_offset1',color='steelblue')
plt.show()

# take a look at the voxel size of the DCE-MRI imagegs
rootdir = '/data/francgrp1'
imdir = '{}/breast_radiomics/her2/MRI'.format(rootdir)
im_info_dir = '{}/MORE_INFO'.format(imdir)
all_dir = [d for d in glob.glob('{}/*'.format(imdir)) if os.path.isdir(d)]
all_pt_id = [int(re.search(r'{}/(\d+)'.format(imdir), ss).group(1)) for ss in all_dir if
             re.search(r'{}/(\d+)'.format(imdir), ss)]

# Go through all the patient data via .dmi
ids = [ii for ii in all_pt_id if ii != 88]  # don't have pt 88 json file for now...

voxsize_df = pd.DataFrame()
for pt_id in ids:
    print 'patient id: {}'.format(pt_id)

    # VOI coordinate info and the time-sorted .DMI files (brtool output)
    findjson = glob.glob('{}/{}_*.json'.format(im_info_dir,pt_id))

    if findjson:
        voi_lps_paths = findjson
    else:
        print 'pt_id {}: VOI JSON file was NOT found!'.format(pt_id)
        continue

    for jj in voi_lps_paths:
        # determine breast side
        info_fname = os.path.basename(jj)
        check = re.search(r'{}_(\w+)_(\d+).info.json'.format(pt_id), info_fname)
        if check:
            breast_side = check.group(1)
        else:
            print 'pt_id {}: cannot determine breast side!'.format(pt_id)
        print 'breast_side: {}'.format(breast_side)

        # initialize a Dataframe for each patient data chunk
        CENTER, HALFLENGTHS,time_sorted_IMNAMES = ImSeg.GetImInfo(jj)

        # DMI MRI volumes (from brtool output)
        pt_series_dir = glob.glob('{}/{:0>3d}/*/*'.format(imdir, pt_id))[0]
        pt_dmi_list = ['{}/dce/{}'.format(pt_series_dir, ff) for ff in time_sorted_IMNAMES]

        print 'image file dmi name: {}'.format(pt_dmi_list[0])
        dceSeries = dmi.DMI(pt_dmi_list[0])
        mri_itk_img = ITKImageHelper.generate_oriented_itkImage(dceSeries)
        IG_mri_img = image_geometry.ImageGeometry(mri_itk_img)
        voxsize = IG_mri_img.samplingRCS

        for ii in range(len(voxsize)):
            tmp_dict = {}
            tmp_dict['pt_id'] = pt_id
            tmp_dict['breast_side'] = breast_side
            tmp_dict['voxel_size'] = voxsize[ii]
            tmp_dict['direction_ii'] = ii
            voxsize_df = voxsize_df.append(tmp_dict,ignore_index=True)

dir_dict = {0:'R',1:'C',2:'S'}
voxsize_df['direction'] = voxsize_df['direction_ii'].apply(lambda x: dir_dict[x])
json_fname = '/data/francgrp1/clare_work/Data/her2_Analysis/DCEMRI_VoxelSize_all.json'
voxsize_df.to_json(json_fname)


# read in the voxsize json file
json_fname = '/data/francgrp1/clare_work/Data/her2_Analysis/DCEMRI_VoxelSize_all.json'
voxsize_df = pd.read_json(json_fname)
print voxsize_df

# perform some aggregate
table = pd.pivot_table(voxsize_df,values='voxel_size',columns='direction',aggfunc=np.amin)
print table

g = sns.FacetGrid(voxsize_df,row='direction',margin_titles=True)
g.map(plt.hist,'voxel_size',color='steelblue')
plt.show()