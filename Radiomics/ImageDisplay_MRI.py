# -*- coding: utf-8 -*-
"""
Created on 9/11/17 2:52 PM

@author: shuang Shih-ying Huang
@goal:

"""


import glob
import ImageProcess as ImP
import dmi
import ITKImageHelper
import re
import os
import image_geometry
import sys



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
    the_select_ids = [dd for dd in ids if dd >= pid_start and dd < pid_end]
    for pt_id in the_select_ids:
        print 'patient id: {}'.format(pt_id)

        # Get the SER mask (.dmi) from the brtool output
        breast_side = 'R'
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
        pt_dmi_list = ['{}/dce/079_E03590S004_sagR_tp1.dmi'.format(pt_series_dir)]

        for ii in range(len(pt_dmi_list)):
            print 'image file dmi name: {}'.format(pt_dmi_list[ii])
            dceSeries = dmi.DMI(pt_dmi_list[ii])
            mri_itk_img = ITKImageHelper.generate_oriented_itkImage(dceSeries)
            mri_itk_sagittal = ITKImageHelper.itkImage_orient_to_sagittal(mri_itk_img)
            IG_mri_img = image_geometry.ImageGeometry(mri_itk_sagittal)
            print 'MRI voxel size (mm^3): {}'.format(IG_mri_img.samplingRCS)

            # visual checking of images and mask
            if is_check_vis == 1:
                ar = IG_mri_img.samplingSRC[2] / IG_mri_img.samplingSRC[1]
                mri_array = ITKImageHelper.itkImage_to_ndarray(mri_itk_sagittal)
                ser_mask_array = ITKImageHelper.itkImage_to_ndarray(ser_mask_sagittal)
                ImP.display_overlay_volume(mri_array, ser_mask_array, 'DCE-MRI image + SER mask', aspect=ar)