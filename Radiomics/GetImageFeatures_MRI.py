#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 9/27/16 10:29 AM

@author: shuang Shih-ying Huang
@goal: Compute image features for a given image input (MRI, PET, etc.)

"""

import ImageFeature as ImF
import glob
import ImageProcess as ImP
import dicom
import dmi
import ITKImageHelper
import pandas as pd
import re
import os
import image_geometry
import sys
import ITKImageFilters as itkif
import resource
import numpy as np


if __name__ == '__main__':

    pid_start = int(sys.argv[1])
    pid_end = int(sys.argv[2])
    is_check_vis = int(sys.argv[3])

    # Lassen
    rootdir = '/data/francgrp1'

    # find patient_id in a list of INT (not string)
    # need to make sure all sub-dir is a directory (not files, as it will pick up duplicate patient id..)
    imdir = '{}/breast_radiomics/her2/MRI'.format(rootdir)
    im_info_dir = '{}/MORE_INFO'.format(imdir)
    all_dir = [d for d in glob.glob('{}/*'.format(imdir)) if os.path.isdir(d)]
    all_pt_id = [int(re.search(r'{}/(\d+)'.format(imdir), ss).group(1)) for ss in all_dir if
                 re.search(r'{}/(\d+)'.format(imdir), ss)]

    # Go through all the patient data via .dmi
    ids = [ii for ii in all_pt_id if ii != 88]  # don't have pt 88 json file for now...

    # image standardization => rescale the raw image into 0-Nbin range
    # img_norm_Nbin = [64, 128, 256]
    # img_norm_Nbin = [256]

    # define the texture features list to compute
    feature_list = ['autocorrelation', 'cluster_prominence', 'cluster_shade', 'cluster_tendency', 'contrast',
                    'correlation', 'diff_entropy', 'dissimilarity', 'energy', 'entropy', 'homogeneity1', 'homogeneity2',
                    'idmn', 'idn', 'inv_var', 'maxprob', 'sum_avg', 'sum_entropy', 'sum_var', 'imc1', 'imc2',
                    'diff_avg','diff_var', 'avg_intensity', 'sum_squares']

    # parameters for the glcm texture analysis
    Rneighbor = 1
    # glcm_Nbin = [64, 128, 256]
    glcm_bin_width = [5]
    the_img_spacing = np.array([0.5,0.5,0.5])  #in millimeter
    proc_version = 2

    mem_use = []
    the_select_ids = [dd for dd in ids if dd >= pid_start and dd < pid_end]

    lst_img_vox_size = []
    col_names = ['pt_id','MRN','img_time_point','voxsize_mm_x','voxsize_mm_y', 'voxsize_mm_z']

    for pt_id in the_select_ids:
        print 'patient id: {}'.format(pt_id)

        # VOI coordinate info and the time-sorted .DMI files (brtool output)
        findjson = glob.glob('{}/{}_*.json'.format(im_info_dir,pt_id))

        if findjson:
            voi_lps_paths = findjson
        else:
            print 'pt_id {}: VOI JSON file was NOT found!'.format(pt_id)
            continue

        for jj in voi_lps_paths:
            # initialize a Dataframe for each patient data chunk
            pt_features_data = pd.DataFrame()

            # Get the images and files needed for MRI
            # determine breast side
            info_fname = os.path.basename(jj)
            check = re.search(r'{}_(\w+)_(\d+).info.json'.format(pt_id),info_fname)
            if check:
                breast_side = check.group(1)
            else:
                print 'pt_id {}: cannot determine breast side!'.format(pt_id)

            if breast_side == 'R' and pt_id == 89:
                print 'pt_id {}, breast_side {}, this case has bad image input..'.format(pt_id,breast_side)
                continue
            else:
                CENTER, HALFLENGTHS,time_sorted_IMNAMES = ImP.GetImInfo(jj)

                # Get the SER mask (.dmi) from the brtool output
                findserdmi = glob.glob('{}/{}_{}*.SER_MASK.dmi'.format(im_info_dir,pt_id,breast_side))
                if findserdmi:
                    ser_mask_fname = findserdmi[0]
                    print 'ser_mask_fname: {}'.format(ser_mask_fname)
                else:
                    print 'pt_id {}: SER mask DMI file was NOT found!'.format(pt_id)
                    continue

                ser_mask_dmi = dmi.DMI(ser_mask_fname)
                ser_mask_itk = ITKImageHelper.generate_oriented_itkImage(ser_mask_dmi)
                ser_mask_sagittal = ITKImageHelper.itkImage_orient_to_sagittal(ser_mask_itk)

                # DMI MRI volumes (from brtool output)
                pt_series_dir = glob.glob('{}/{:0>3d}/*/*'.format(imdir, pt_id))[0]
                pt_dmi_list = ['{}/dce/{}'.format(pt_series_dir, ff) for ff in time_sorted_IMNAMES]

                for ii in range(len(pt_dmi_list)):
                    print 'image file dmi name: {}'.format(pt_dmi_list[ii])
                    dceSeries = dmi.DMI(pt_dmi_list[ii])
                    mri_itk_img = ITKImageHelper.generate_oriented_itkImage(dceSeries)
                    mri_itk_sagittal = ITKImageHelper.itkImage_orient_to_sagittal(mri_itk_img)
                    IG_mri_img = image_geometry.ImageGeometry(mri_itk_sagittal)
                    print 'MRI voxel size (mm^3): {}'.format(IG_mri_img.samplingRCS)
                    # # visual checking of images and mask
                    # if is_check_vis == 1:
                    #     ar = IG_mri_img.samplingSRC[2] / IG_mri_img.samplingSRC[1]
                    #     mri_array = ITKImageHelper.itkImage_to_ndarray(mri_itk_sagittal)
                    #     ser_mask_array = ITKImageHelper.itkImage_to_ndarray(ser_mask_sagittal)
                    #     ImP.display_overlay_volume(mri_array, ser_mask_array, 'DCE-MRI image + SER mask', aspect=ar)

                    # figure out the patient MRN
                    tag_mrn = dicom.tag.Tag(0x0010, 0x0020)
                    pt_mrn = dceSeries._DS[tag_mrn].value

                    tag_acc_num = dicom.tag.Tag(0x0008,0x0050)
                    pt_acc_num = dceSeries._DS[tag_acc_num].value
                    print 'pt_mrn: {}, pt_acc_num: {}'.format(pt_mrn,pt_acc_num)

                    # figure out the MR time point
                    fname = pt_dmi_list[ii]
                    the_tp = os.path.basename(os.path.splitext(fname)[0]).split('_')[-1]
                    tmp = dict(zip(col_names, [pt_id, pt_mrn, the_tp] + IG_mri_img.samplingRCS))
                    lst_img_vox_size.append(tmp)

        #             if ii == 0:
        #                 # get VOI info
        #                 voi_ixyz,voi_size = ImP.GetVOIinfo(CENTER,HALFLENGTHS,IG_mri_img)
        #
        #                 # Get mask VOI (only need to do it once)
        #                 ser_mask_sagittal_roi = ImP.GetITKVOI(ser_mask_sagittal, voi_size, voi_ixyz)
        #
        #                 # resample the mask VOI to the desired voxel spacing
        #                 resampled_itk_mask = itkif.ITKImageResample(ser_mask_sagittal_roi,the_img_spacing,is_mask=True)
        #
        #
        #             # get the image VOI
        #             mri_itk_sagittal_roi = ImP.GetITKVOI(mri_itk_sagittal, voi_size, voi_ixyz)
        #
        #             # ## visual checking of images and mask
        #             # if is_check_vis == 1:
        #             #     IG_roi_itk_img = image_geometry.ImageGeometry(mri_itk_sagittal_roi)
        #             #     ar = IG_roi_itk_img.samplingSRC[2] / IG_roi_itk_img.samplingSRC[1]
        #             #     img_array = ITKImageHelper.itkImage_to_ndarray(mri_itk_sagittal_roi)
        #             #     mask_array = ITKImageHelper.itkImage_to_ndarray(ser_mask_sagittal_roi)
        #             #     ImP.display_overlay_volume(img_array, mask_array,'DCE-MRI image VOI + SER mask VOI', aspect=ar)
        #
        #             # print out input image geometry info
        #             ITKImageHelper.itkImage_print_image_geometry(mri_itk_sagittal_roi)
        #
        #             # resample the image VOI to the desired voxel spacing
        #             resampled_itk_img = itkif.ITKImageResample(mri_itk_sagittal_roi,the_img_spacing,is_mask=False)
        #             IG_resampled_itk_img = image_geometry.ImageGeometry(resampled_itk_img)
        #
        #             # visual checking of images and mask
        #             if ii == 0 and is_check_vis == 1:
        #                 output_dir = '{}/{:0>3d}/img_overlay'.format(imdir, pt_id)
        #                 print(output_dir)
        #                 ImP.ITK_Image_OverlayPlot(resampled_itk_img, resampled_itk_mask, ' Resampled DCE-MRI image + SER mask',
        #                                           output_dir=output_dir, filetag='MR_DCE_{}_{}'.format(pt_id, breast_side))
        #
        #             print 'mask resampled voxel size (mm^3): {}'.format(IG_resampled_itk_img.samplingRCS)
        #             # define texture feature list
        #             theimf = ImF.ImageFeature(resampled_itk_img,feature_list,resampled_itk_mask,'dict')
        #
        #             for bb in glcm_bin_width:
        #                 print 'glcm bin width: {}'.format(bb)
        #                 theimf._compute_texture_features(Rneighbor, binWidth=bb)
        #                 theimf._compute_first_order_stats(binWidth=bb)
        #                 theimf._compute_shape_size_features()
        #
        #                 tmp_dict = theimf.feature_output
        #                 tmp_dict['pt_id'] = pt_id
        #                 tmp_dict['pt_mrn'] = pt_mrn
        #                 tmp_dict['pt_accession_num'] = pt_acc_num
        #                 tmp_dict['dce_series_dmi_fname'] = pt_dmi_list[ii]
        #                 tmp_dict['glcm_bin_width'] = bb
        #                 tmp_dict['organ_mask'] = 'breast tumor'
        #                 tmp_dict['process_name'] = 'GetImageFeature_VOI'
        #                 tmp_dict['process_version'] = proc_version
        #                 tmp_dict['voxel_size_mm3'] = the_img_spacing
        #
        #                 #TODO: instead of writing to a dataframe, consider write to the mongodb db
        #                 pt_features_data = pt_features_data.append(tmp_dict, ignore_index=True)
        #                 # print pt_features_data
        #
        #     # monitor the memory usage
        #     mem_use.append(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        #
        #     pt_features_data['pt_id'] = pt_features_data['pt_id'].astype('int')
        #     save_dir = '{}/clare_work/Data/her2_ImageFeatures/IsoVoxelSize'.format(rootdir)
        #     dataoutfname = '{}/MRI_features_data_v{}_{}_{}.json'.format(save_dir,proc_version, pt_id,breast_side)
        #     pt_features_data.to_json(dataoutfname)
        #
        # # save the memory log
        # outfile = open('{}/memory_usage_v{}_{}_{}.txt'.format(save_dir,proc_version, pid_start,pid_end), 'w')
        # outfile.write('\n'.join([str(s) for s in mem_use]))
        # outfile.close()

        save_dir = '{}/clare_work/Data/her2_ImageFeatures/IsoVoxelSize'.format(rootdir)
        df_data_all = pd.DataFrame(lst_img_vox_size)
        df_data_all = df_data_all.ix[:, col_names]
        fname = '{}/MR_img_vox_size_mm.csv'.format(save_dir)
        df_data_all.to_csv(fname)