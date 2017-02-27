#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 7/12/16 11:17 AM

@author: shuang Shih-ying Huang

"""

import ImageFeature as ImF
import glob
import ImageProcess as ImSeg
import dicom
import dmi
import ITKImageHelper
import pandas as pd
import re
import os
import image_geometry
import sys


if __name__ == '__main__':

    pid_start = int(sys.argv[1])

    # # mac computer
    # rootdir = '/Users/shuang/Documents/Proj_Radiomics/IM'

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

    # pt_id's with bilateral tumors
    # pt_id = 79 doesn't have R side VOI info and etc. and its L side VOI info is not readable for some reason...
    # pt_id = 89, the R side mask is not good or brought up 'a masked constant error message during feature computation
    # ids = [84]

    # image standardization => rescale the raw image into 0-Nbin range
    img_norm_Nbin = 256

    # define the texture features list to compute
    feature_list = ['autocorrelation', 'cluster_prominence', 'cluster_shade', 'cluster_tendency', 'contrast', 'correlation',
                    'diff_entropy', 'dissimilarity', 'energy', 'entropy', 'homogeneity1', 'homogeneity2', 'idmn',
                    'idn', 'inv_var', 'maxprob', 'sum_avg', 'sum_entropy', 'sum_var']

    # parameters for the glcm texture analysis
    Rneighbor = 1
    glcm_Nbin = [64, 128, 256]

    the_select_ids = [dd for dd in ids if dd > pid_start]
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

            # determine breast side
            info_fname = os.path.basename(jj)
            check = re.search(r'{}_(\w+)_(\d+).info.json'.format(pt_id),info_fname)
            if check:
                breast_side = check.group(1)
            else:
                print 'pt_id {}: cannot determine breast side!'.format(pt_id)

            CENTER, HALFLENGTHS,time_sorted_IMNAMES = ImSeg.GetImInfo(jj)

            # Get the SER mask (.dmi) from the brtool output
            findserdmi = glob.glob('{}/{}_{}*.SER_MASK.dmi'.format(im_info_dir,pt_id,breast_side))
            if findserdmi:
                ser_mask_fname = findserdmi[0]
            else:
                print 'pt_id {}: SER mask DMI file was NOT found!'.format(pt_id)
                continue

            ser_mask_itk = ITKImageHelper.generate_oriented_itkImage(dmi.DMI(ser_mask_fname))
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
                # print IG_mri_img.samplingRCS

                # display image with mask for checking
                ar = IG_mri_img.samplingSRC[2] / IG_mri_img.samplingSRC[1]

                mri_array = ITKImageHelper.itkImage_to_ndarray(mri_itk_sagittal)
                ser_mask_array = ITKImageHelper.itkImage_to_ndarray(ser_mask_sagittal)

                # visual checking of images and mask
                ImSeg.display_overlay_volume(mri_array, ser_mask_array, 'DCE-MRI image + SER mask', aspect=ar)


                if ii == 0:
                    # figure out the patient MRN
                    tag_mrn = dicom.tag.Tag(0x0010, 0x0020)
                    pt_mrn = dceSeries._DS[tag_mrn].value

                    tag_acc_num = dicom.tag.Tag(0x0008,0x0050)
                    pt_acc_num = dceSeries._DS[tag_acc_num].value

                    # get VOI info
                    voi_ixyz,voi_size = ImSeg.GetVOIinfo(CENTER,HALFLENGTHS,IG_mri_img)

                    # Get mask VOI (only need to do it once)
                    ser_mask_sagittal_roi = ImSeg.GetITKVOI(ser_mask_sagittal, voi_size, voi_ixyz)
                    ser_mask_array = ITKImageHelper.itkImage_to_ndarray(ser_mask_sagittal_roi)

                # get the image VOI
                mri_itk_sagittal_roi = ImSeg.GetITKVOI(mri_itk_sagittal, voi_size, voi_ixyz)
                mri_array = ITKImageHelper.itkImage_to_ndarray(mri_itk_sagittal_roi)

                # rescale the MRI VOI image
                mri_roi_rescale_itk = ImSeg.ITKStandardizeImageIntensity(mri_itk_sagittal_roi,img_norm_Nbin)

                # # visual checking of images and mask
                # ImSeg.display_overlay_volume(mri_array, ser_mask_array, 'DCE-MRI image VOI + SER mask VOI', aspect=ar)

                # define texture feature list
                theimf = ImF.ImageFeature(mri_roi_rescale_itk, feature_list, ser_mask_sagittal_roi)

                for bb in glcm_Nbin:
                    print 'glcm Nbin: {}'.format(bb)
                    theimf._compute_texture_features(Rneighbor,bb)
                    theimf._compute_first_order_stats()
                    theimf._compute_shape_size_features()

                    tmp_dict = theimf.feature_output
                    tmp_dict['pt_id'] = pt_id
                    tmp_dict['pt_mrn'] = pt_mrn
                    tmp_dict['pt_accession_num'] = pt_acc_num
                    tmp_dict['dce_series_dmi_fname'] = pt_dmi_list[ii]
                    tmp_dict['glcm_Nbin'] = bb
                    tmp_dict['organ_mask'] = 'breast tumor'
                    tmp_dict['process_name'] = 'GetImageFeature_VOI'
                    tmp_dict['process_version'] = '1.0'
                    tmp_dict['img_norm_Nbin'] = img_norm_Nbin

                    #TODO: instead of writing to a dataframe, consider write to the mongodb db
                    pt_features_data = pt_features_data.append(tmp_dict, ignore_index=True)

            pt_features_data['pt_id'] = pt_features_data['pt_id'].astype('int')
            dataoutfname = '{}/clare_work/Data/MRI_features_data_{}_{}.json'.format(rootdir,pt_id,breast_side)
            pt_features_data.to_json(dataoutfname)

