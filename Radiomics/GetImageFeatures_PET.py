#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 9/29/16 3:12 PM

@author: shuang Shih-ying Huang
@goal: Compute image features from the PET images with the manual segmenated volume (via MeVisLab by R. Harnish)

"""


import ImageFeature as ImF
import glob
import ImageProcess as ImP
import dicom_series
import ITKImageHelper
import re
import os
import image_geometry
import sys
import ITKImageFilters as itkif
import resource
import pet_util
import numpy as np
from dicom.tag import Tag
import pandas as pd


if __name__ == '__main__':

    pid_start = int(sys.argv[1])
    pid_end = int(sys.argv[2])
    is_vis_check = int(sys.argv[3])

    # Lassen
    rootdir = '/data/francgrp1'

    # get all patient ID's in the PETCT directory
    imdir = '{}/breast_radiomics/her2/ALL_TUMORS'.format(rootdir)
    all_dir = [d for d in glob.glob('{}/*'.format(imdir)) if os.path.isdir(d)]
    all_pt_id = [int(re.search(r'{}/(\d+)'.format(imdir), ss).group(1)) for ss in all_dir if
                 re.search(r'{}/(\d+)'.format(imdir), ss)]
    mask_dir = '{}/mevis_manual_segmentation'.format(imdir)

    # image standardization => rescale the raw image into 0-Nbin range
    # img_norm_Nbin = 256
    # img_norm_Nbin = [32,64,128]
    # define the texture features list to compute
    # feature_list = ['autocorrelation', 'cluster_prominence', 'cluster_shade', 'cluster_tendency', 'contrast', 'correlation',
    #                 'diff_entropy', 'dissimilarity', 'energy', 'entropy', 'homogeneity1', 'homogeneity2', 'idmn',
    #                 'idn', 'inv_var', 'maxprob', 'sum_avg', 'sum_entropy', 'sum_var']

    feature_list = ['autocorrelation', 'cluster_prominence', 'cluster_shade', 'cluster_tendency', 'contrast',
                    'correlation','diff_entropy', 'dissimilarity', 'energy', 'entropy', 'homogeneity1', 'homogeneity2',
                    'idmn','idn', 'inv_var', 'maxprob', 'sum_avg', 'sum_entropy', 'sum_var', 'imc1', 'imc2', 'diff_avg',
                    'diff_var', 'avg_intensity', 'sum_squares']

    # parameters for the glcm texture analysis
    Rneighbor = 1
    # glcm_Nbin = [32,64,128]
    # glcm_Nbin = [128]
    glcm_bin_width = [0.1]
    proc_version = 2

    the_img_spacingRCS = np.array([2.,2.,2.]) #in millimeter
    fragmented_dcm_ids = [16,17,18,22,25,27,30]

    mem_use = []
    the_select_ids = [dd for dd in all_pt_id if dd >= pid_start and dd < pid_end]
    save_dir = '{}/clare_work/Data/her2_ImageFeatures/IsoVoxelSize'.format(rootdir)

    lst_img_vox_size = []
    col_names = ['pt_id', 'MRN', 'voxsize_mm_x', 'voxsize_mm_y', 'voxsize_mm_z']

    for pt_id in the_select_ids:
        print 'patient id: {}'.format(pt_id)

        # PETCT image dir
        pet_fdirs = [d for d in glob.glob('{}/{:0>3d}/*/*'.format(imdir, pt_id)) if os.path.isdir(d) and len(d.split('/')[-1]) > 40]

        for pet_fdir in pet_fdirs:
            check = re.search(r'(.+)/images_(\w+)/(.+)', pet_fdir)
            if check:
                breast_side = check.group(2)
                print 'breast side: {}'.format(breast_side)
            else:
                print 'Cannot determine the breast side of the PET images'
                continue

            # Tumor mask .nrrd
            check_mask_fname = glob.glob('{}/{:0>3d}*sag{}*_uc_origin.nrrd'.format(mask_dir,pt_id, breast_side))
            if check_mask_fname:
                mask_fname = check_mask_fname[0]
            else:
                print 'Cannot find the mask .nrrd file with the breast side!'
                continue

            # determine the correct PET series
            pet_img_series = dicom_series.read_files(pet_fdir)
            pet_ds_list = [ss for ss in pet_img_series if len(ss.shape) > 2]
            print(pet_ds_list)
            if len(pet_ds_list) > 1:
                print 'Warning! more than ONE dicom_series is found...{}'.format(pet_ds_list)
                # find the appropriate dicom_series to look at
                dc_lens = [ss.shape[0] for ss in pet_ds_list]
                pet_ds = pet_ds_list[dc_lens.index(max(dc_lens))]
            else:
                pet_ds = pet_ds_list[0]
            IG_pet = image_geometry.ImageGeometry(pet_ds)

            if pt_id not in fragmented_dcm_ids:
                # use pet_ds for image analysis since it's not a fragmented cases
                PIX_pet = pet_ds.get_pixel_array()
                pet_itk_img = ITKImageHelper.generate_oriented_itkImage(ig=IG_pet, pixarray=PIX_pet)
            else: # use the mevistlab exported NRRD file for the fragmented dicom series (dicom_series and mevislab did not parse the dicom in the same way for some reason...)
                print '::O_O:: Fragmented dicom series! read the MeVisLab exported NRRD instead of dicom series!'
                pet_nrrd_fname = '{}/PET_voxshiftON_rescale.nrrd'.format(os.path.dirname(pet_fdir))
                pet_itk_img = ITKImageHelper.itkImage_read(pet_nrrd_fname)

            # need to orient the PETCT to 'AXIAL' for correct mask and image registration
            petct_itk_img_orient = ITKImageHelper.itkImage_orient_to_axial(pet_itk_img)

            # Convert raw PET activity concentration to SUV
            suv_factor = pet_util.get_SUV_multiplier(pet_ds)
            if suv_factor == 0:
                print '::O_O:: cannot calculate the SUV factor!'
                continue
            else:
                print 'suv_factor = {}'.format(suv_factor)
                petct_suv_itk = itkif.ITKImageMultiplyConstant(petct_itk_img_orient,suv_factor)

            # get the mask ITK image
            mask_itk = ITKImageHelper.itkImage_read(mask_fname)
            # not sure if one needs to orient it to 'axial' for the mask_itk to work well...
            mask_itk_orient = ITKImageHelper.itkImage_orient_to_axial(mask_itk)

            # get patient info
            tag_mrn = Tag(0x0010, 0x0020)
            pt_mrn = pet_ds.info[tag_mrn].value

            tag_acc_num = Tag(0x0008, 0x0050)
            pt_acc_num = pet_ds.info[tag_acc_num].value

            tmp = dict(zip(col_names, [pt_id, pt_mrn] + IG_pet.samplingRCS))
            lst_img_vox_size.append(tmp)

            # if is_vis_check == 1:
            #     # visual checking of SUV-PET images and the tumor mask
            #     ImP.ITK_Image_OverlayPlot(petct_suv_itk,mask_itk,'SUV-PET image + Tumor mask')

    #         # resample the mask VOI to the desired voxel spacing & determine if the images need to be resampled..
    #         check_resample = (np.array(IG_pet.samplingRCS) > the_img_spacingRCS).all()
    #         if check_resample:
    #             print 'PET image voxsel spacing is coarser than the desired spacing, {}'.format(IG_pet.samplingRCS)
    #             final_itk_mask = itkif.ITKImageResample(mask_itk, the_img_spacingRCS,is_mask=True)
    #             final_itk_img = itkif.ITKImageResample(petct_suv_itk,the_img_spacingRCS,is_mask=False)
    #         else:
    #             print 'the image voxel spacing is the same as the desired voxel spacing! NO need to resample the image!'
    #             final_itk_img = petct_suv_itk
    #             final_itk_mask = mask_itk
    #
    #         if is_vis_check == 1:
    #             # visual checking of the final processed PET images and the tumor mask
    #             # ImP.ITK_Image_OverlayPlot(final_itk_img,final_itk_mask,'Post Resampling + SUV PET image + Tumor mask')
    #             output_dir = '{}/{:0>3d}/images_{}/img_overlay'.format(imdir,pt_id,breast_side)
    #             print(output_dir)
    #             ImP.ITK_Image_OverlayPlot(final_itk_img, final_itk_mask, 'Post Resampling + SUV PET image',
    #                                       vmax_s=0.3, output_dir=output_dir, filetag='PET_{}_{}'.format(pt_id, breast_side))
    #         # define texture feature list
    #         theimf = ImF.ImageFeature(final_itk_img,feature_list,final_itk_mask,'dict')
    #
    #         pt_features_data = pd.DataFrame()
    #         for bb in glcm_bin_width:
    #             print 'bin width: {}'.format(bb)
    #             theimf._compute_texture_features(Rneighbor, binWidth=bb)
    #             theimf._compute_first_order_stats(binWidth=bb)
    #             theimf._compute_shape_size_features()
    #
    #             tmp_dict = theimf.feature_output
    #             tmp_dict['pt_id'] = pt_id
    #             tmp_dict['pt_mrn'] = pt_mrn
    #             tmp_dict['pt_accession_num'] = pt_acc_num
    #             tmp_dict['pet_series_fdir'] = pet_fdir
    #             tmp_dict['bin_width'] = bb
    #             tmp_dict['organ_mask'] = 'breast tumor'
    #             tmp_dict['process_name'] = 'GetImageFeature_VOI'
    #             tmp_dict['process_version'] = proc_version
    #             tmp_dict['voxel_size_mm3'] = [tuple(the_img_spacingRCS)] * len(tmp_dict)
    #             tmp_dict['breast_side'] = breast_side
    #             # print theimf.feature_output
    #
    #             #TODO: instead of writing to a dataframe, consider write to the mongodb db
    #             pt_features_data = pt_features_data.append(tmp_dict, ignore_index=True)
    #             # print pt_features_data
    #
    #         # monitor the memory usage
    #         mem_use.append(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    #
    #         pt_features_data['pt_id'] = pt_features_data['pt_id'].astype('int')
    #         dataoutfname = '{}/PET_features_data_v{}_{}_{}.json'.format(save_dir,proc_version,pt_id,breast_side)
    #         pt_features_data.to_json(dataoutfname)
    #
    # # save the memory log
    # outfile = open('{}/PET_Features_v{}_memory_usage_{}_{}.txt'.format(save_dir,proc_version, pid_start,pid_end), 'w')
    # outfile.write('\n'.join([str(s) for s in mem_use]))
    # outfile.close()

    save_dir = '{}/clare_work/Data/her2_ImageFeatures/IsoVoxelSize'.format(rootdir)
    df_data_all = pd.DataFrame(lst_img_vox_size)
    df_data_all = df_data_all.ix[:, col_names]
    fname = '{}/PET_img_vox_size_mm.csv'.format(save_dir)
    df_data_all.to_csv(fname)

#TODO:
# check on 79_R and make sure the mask exists!!!! ==> MRI has a 79_R mask but it doesn't seem to cover the breast section
# one can run ImageDisplay_MRI.py with 79 80 1 as the input arg to test to see if the image and mask make sense

# TESTING the dicom series vs MeVisLab NRRD files
# pet_itk_fname1 = '/data/francgrp1/TMP/MeVisOriginObfuscation/084_correct_sub_voxel_shift_ON.nrrd'
# petct_itk_img1 = ITKImageHelper.itkImage_read(pet_itk_fname1)
#
# pet_itk_fname2 = '/data/francgrp1/TMP/MeVisOriginObfuscation/084_correct_sub_voxel_shift_OFF.nrrd'
# petct_itk_img2 = ITKImageHelper.itkImage_read(pet_itk_fname2)
#
# pet_itk_fname3 = '/data/francgrp1/TMP/MeVisOriginObfuscation/084_correct_sub_voxel_shift_ON_rescaled.nrrd'
# petct_itk_img3 = ITKImageHelper.itkImage_read(pet_itk_fname3)
#
# # save the numpy arrays to testing
# img_np_fname = '/data/francgrp1/clare_work/tmp/img_array_nrrd_voxshiftON.npy'
# img_array1 = ITKImageHelper.itkImage_to_ndarray(petct_itk_img1)
# im1 = np.zeros(img_array1.shape)
# im1[img_array1 != 0] = img_array1[img_array1 != 0]
# np.save(img_np_fname, im1)
#
# img_np_fname = '/data/francgrp1/clare_work/tmp/img_array_nrrd_voxshiftOFF.npy'
# img_array2 = ITKImageHelper.itkImage_to_ndarray(petct_itk_img2)
# im2 = np.zeros(img_array2.shape)
# im2[img_array2 != 0] = img_array2[img_array2 != 0]
# np.save(img_np_fname, im2)
#
# img_np_fname = '/data/francgrp1/clare_work/tmp/img_array_dicomseries.npy'
# img_array3 = ITKImageHelper.itkImage_to_ndarray(petct_itk_img)
# im3 = np.zeros(img_array3.shape)
# im3[img_array3 != 0] = img_array3[img_array3 != 0]
# np.save(img_np_fname, im3)
#
# mask_np_fname = '/data/francgrp1/clare_work/tmp/mask_array.npy'
# im4 = np.zeros(mask_array.shape)
# im4[mask_array != 0] = mask_array[mask_array != 0]
# np.save(mask_np_fname, im4)
#
# img_np_fname = '/data/francgrp1/clare_work/tmp/img_array_nrrd_voxshiftON_rescale.npy'
# img_array4 = ITKImageHelper.itkImage_to_ndarray(petct_itk_img3)
# im4 = np.zeros(img_array4.shape)
# im4[img_array4 != 0] = img_array4[img_array4 != 0]
# np.save(img_np_fname, im4)
