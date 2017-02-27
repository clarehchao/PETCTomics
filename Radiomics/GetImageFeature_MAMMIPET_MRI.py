# -*- coding: utf-8 -*-
"""
Created on 11/28/16 10:07 AM

@author: shuang Shih-ying Huang
@goal: extract features for MAMMIPET and MRI

"""

import dicom_series
import image_geometry
import ITKImageHelper
import ImageFeature as ImF
import ImageProcess as ImP
import pandas as pd
import nrrd

rootdir = '/data/francgrp1/breast_radiomics/MAMMIPET/PT001'

# # load in MRI images
# Rneighbor = 1
# glcm_Nbin = [128, 256]
# img_norm_Nbin = [64, 128, 256]
# tumor_fname_tag = 'DCE601'
# img_type = 'MRI'
# img_fname = '{}/MRI/MRI_voxshiftON_rescale.nrrd'.format(rootdir)
# itk_img = ITKImageHelper.itkImage_read(img_fname)
# itk_img_orient = ITKImageHelper.itkImage_orient_to_axial(itk_img)
# ig_img = image_geometry.ImageGeometry(itk_img_orient)


# load MAMMIPET images
Rneighbor = 1
glcm_Nbin = [128, 256]
img_norm_Nbin = [64, 128, 256]
tumor_fname_tag = 'PET_IM-0001-0359'
img_type = 'PET'
img_fname = '{}/PET/PET_voxshiftON_rescale.nrrd'.format(rootdir)
itk_img = ITKImageHelper.itkImage_read(img_fname)
itk_img_orient = ITKImageHelper.itkImage_orient_to_axial(itk_img)
ig_img = image_geometry.ImageGeometry(itk_img_orient)


# load the tumor mask for MRI
iTumors = [1,2]
for ii in iTumors:
    tumor_fname = '{}/TumorMask/{}_Tumor{}.nrrd'.format(rootdir,tumor_fname_tag,ii)
    itk_tumor_mask = ITKImageHelper.itkImage_read(tumor_fname)
    itk_tumor_orient = ITKImageHelper.itkImage_orient_to_axial(itk_tumor_mask)
    ig_tumor = image_geometry.ImageGeometry(itk_tumor_orient)
    # print ig_tumor.samplingRCS
    # print ig_mri.direction_cosine_matrix

    # # visualize the volumn
    # ImP.ITK_Image_OverlayPlot(itk_img_orient,itk_tumor_orient,'MAMMIPET patient {} images, Tumor #{}'.format(img_type,ii))

    # get image feature
    feature_list = ['autocorrelation', 'cluster_prominence', 'cluster_shade', 'cluster_tendency', 'contrast', 'correlation',
                        'diff_entropy', 'dissimilarity', 'energy', 'entropy', 'homogeneity1', 'homogeneity2', 'idmn',
                        'idn', 'inv_var', 'maxprob', 'sum_avg', 'sum_entropy', 'sum_var']

    theimf = ImF.ImageFeature(itk_img_orient,feature_list,itk_tumor_orient)


    pt_features_data = pd.DataFrame()
    for aa in img_norm_Nbin:
        for bb in glcm_Nbin:
            print 'img_Nbin: {}, glcm Nbin: {}'.format(aa,bb)
            theimf._compute_texture_features(Rneighbor, GSnorm_Nbin=aa, glcm_Nbin=bb)
            theimf._compute_first_order_stats()
            theimf._compute_shape_size_features()

            tmp_df = theimf.feature_output
            tmp_df['glcm_Nbin'] = bb
            tmp_df['img_norm_Nbin'] = aa
            tmp_df['organ_mask'] = 'breast tumor'
            tmp_df['process_name'] = 'GetImageFeature_VOI'
            tmp_df['process_version'] = '1.0'
            tmp_df['voxel_size_mm3'] = [tuple(ig_img.samplingRCS)] * len(tmp_df)

            # TODO: instead of writing to a dataframe, consider write to the mongodb db
            pt_features_data = pt_features_data.append(tmp_df, ignore_index=True)

    dataoutfname = '{}/{}_features_data_Tumor{}.csv'.format(rootdir,img_type,ii)
    pt_features_data.to_csv(dataoutfname)



