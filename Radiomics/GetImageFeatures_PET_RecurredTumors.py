#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 5/30/18

@author: shuang
@goal: get the radiomics of PET images of the recurred tumors

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
from pydicom.tag import Tag
import resource
import pet_util
import numpy as np
import pandas as pd


if __name__ == '__main__':
    rootdir = sys.argv[1]
    is_vis_check = int(sys.argv[2])

    feature_list = ['autocorrelation', 'cluster_prominence', 'cluster_shade', 'cluster_tendency', 'contrast',
                    'correlation','diff_entropy', 'dissimilarity', 'energy', 'entropy', 'homogeneity1', 'homogeneity2',
                    'idmn','idn', 'inv_var', 'maxprob', 'sum_avg', 'sum_entropy', 'sum_var', 'imc1', 'imc2', 'diff_avg',
                    'diff_var', 'avg_intensity', 'sum_squares']

    # parameters for the glcm texture analysis
    Rneighbor = 1
    glcm_bin_width = [0.1]
    proc_version = 2
    the_img_spacingRCS = np.array([2.,2.,2.]) #in millimeter

    mem_use = []
    lst_img_vox_size = []
    col_names = ['pt_id', 'MRN', 'voxsize_mm_x', 'voxsize_mm_y', 'voxsize_mm_z']

    # get all the PET image dirs
    imdir = '{}/Data/her2/her2_IM/breastCA_pet_data'.format(rootdir)

    #data saving directory
    save_dir = '{}/Data/her2/her2_ImageFeatures/IsoVoxelSize'.format(rootdir)

    all_img_dirs = os.walk(imdir).next()[1]
    for ptid_side in all_img_dirs:
        print('patient id + side: {}'.format(ptid_side))
        pt_id, breast_side = ptid_side.split('_')

        # PETCT image dir
        pet_fdir = '{}/{}/PET'.format(imdir, ptid_side)

        pet_img_series = dicom_series.read_files(pet_fdir)
        pet_ds_list = [ss for ss in pet_img_series if len(ss.shape) > 2]
        if len(pet_ds_list) > 1:
            print('Warning! more than ONE dicom_series is found...{}'.format(pet_ds_list))
            # find the appropriate dicom_series to look at
            dc_lens = [ss.shape[0] for ss in pet_ds_list]
            pet_ds = pet_ds_list[dc_lens.index(max(dc_lens))]
        else:
            pet_ds = pet_ds_list[0]

        print('pet dicom series: {}'.format(pet_ds))
        series_des = pet_ds.info.SeriesDescription

        # # PET image geometry
        # IG_pet = image_geometry.ImageGeometry(pet_ds)
        #
        # # PET image ITK instance
        # PIX_pet = pet_ds.get_pixel_array()
        # pet_itk_img = ITKImageHelper.generate_oriented_itkImage(ig=IG_pet, pixarray=PIX_pet)
        # # print(ITKImageHelper.itkImage_print_image_geometry(pet_itk_img))

        # read the PET nrrd file from 3dslicer
        # NOTE: use the nrrd file for PET images ouputted from 3dslicer so the voxel size is consistent
        # (it caused inconsistent re-sampling when using the PET dicom images via dicom_series directly
        slicer_dir = '{}/{}/3dslicer'.format(imdir, ptid_side)
        all_files = glob.glob('{}/*.nrrd'.format(slicer_dir))
        tmp_lst = [ff for ff in all_files if re.search('{}/([\w\s.]+){}_f32.nrrd'.format(slicer_dir,re.escape(series_des)), ff)]
        if not tmp_lst:
            print('cannot find any PET nrrd file!')
            continue
        else:
            pet_nrrd_fname = tmp_lst[0]
	    print(pet_nrrd_fname)
            pet_itk_img = ITKImageHelper.itkImage_read(pet_nrrd_fname)
            pet_itk_img_orient = ITKImageHelper.itkImage_orient_to_axial(pet_itk_img)
            print(ITKImageHelper.itkImage_print_image_geometry(pet_itk_img_orient))

            IG_pet = image_geometry.ImageGeometry(pet_itk_img_orient)

            # get patient info from PET
            tag_mrn = Tag(0x0010, 0x0020)
            pt_mrn = pet_ds.info[tag_mrn].value

            tag_acc_num = Tag(0x0008, 0x0050)
            pt_acc_num = pet_ds.info[tag_acc_num].value

            tmp = dict(zip(col_names, [pt_id, pt_mrn] + IG_pet.samplingRCS))
            lst_img_vox_size.append(tmp)

            # Convert raw PET activity concentration to SUV
            suv_factor = pet_util.get_SUV_multiplier(pet_ds)
            if suv_factor == 0:
                print('::O_O:: cannot calculate the SUV factor!')
                continue
            else:
                print('suv_factor = {}'.format(suv_factor))
                pet_suv_itk = itkif.ITKImageMultiplyConstant(pet_itk_img_orient,suv_factor)

            check_resample = (np.array(IG_pet.samplingRCS) > the_img_spacingRCS).all()
            if check_resample:
                print('PET image voxsel spacing is coarser than the desired spacing, {}'.format(IG_pet.samplingRCS))
                final_itk_img = itkif.ITKImageResample(pet_suv_itk, the_img_spacingRCS, is_mask=False)
            else:
                print('the image voxel spacing is the same as the desired voxel spacing! NO need to resample the image!')
                final_itk_img = pet_suv_itk

            # get all the segmented tumors
            tumor_dir = '{}/{}/3dslicer'.format(imdir, ptid_side)
            all_tumors = glob.glob('{}/*-label.nrrd'.format(tumor_dir))

            for tumor_nrrd_fname in all_tumors:
                match = re.search(r'{}/Tumor-([.\w]+)-label.nrrd'.format(tumor_dir),tumor_nrrd_fname)
                if match:
                    tumor_tag = match.group(1)
                else:
                    print('Cannot determine tumor location from the .nrrd file!')

                print('tumor_tag: {}'.format(tumor_tag))

                mask_itk = ITKImageHelper.itkImage_read(tumor_nrrd_fname)
                # not sure if one needs to orient it to 'axial' for the mask_itk to work well...
                mask_itk_orient = ITKImageHelper.itkImage_orient_to_axial(mask_itk)
                print(ITKImageHelper.itkImage_print_image_geometry(mask_itk_orient))

                IG_mask = image_geometry.ImageGeometry(mask_itk_orient)

                # resample the mask VOI to the desired voxel spacing & determine if the images need to be resampled..
                check_resample = (np.array(IG_mask.samplingRCS) > the_img_spacingRCS).all()
                if check_resample:
                    print('PET image voxsel spacing is coarser than the desired spacing, {}'.format(IG_mask.samplingRCS))
                    final_itk_mask = itkif.ITKImageResample(mask_itk_orient, the_img_spacingRCS,is_mask=True)
                else:
                    print('the image voxel spacing is the same as the desired voxel spacing! NO need to resample the image!')
                    final_itk_mask = mask_itk_orient

                if is_vis_check == 1:
                    # visual checking of the final processed PET images and the tumor mask
                    # ImP.ITK_Image_OverlayPlot(final_itk_img,final_itk_mask,'Post Resampling + SUV PET image + Tumor mask')
                    print('in vis check statement!')
                    ImP.ITK_Image_OverlayPlot(final_itk_img, final_itk_mask, 'SUV-PET image + Tumor mask')

                print(ITKImageHelper.itkImage_print_image_geometry(final_itk_img))
                print(ITKImageHelper.itkImage_print_image_geometry(final_itk_mask))

                # define texture feature list
                theimf = ImF.ImageFeature(final_itk_img,feature_list,final_itk_mask,'dict')

                pt_features_data = pd.DataFrame()
                for bb in glcm_bin_width:
                    print('bin width: {}'.format(bb))
                    theimf._compute_texture_features(Rneighbor, binWidth=bb)
                    theimf._compute_first_order_stats(binWidth=bb)
                    theimf._compute_shape_size_features()

                    tmp_dict = theimf.feature_output
                    tmp_dict['pt_id'] = pt_id
                    tmp_dict['pt_mrn'] = pt_mrn
                    tmp_dict['pt_accession_num'] = pt_acc_num
                    tmp_dict['pet_series_fdir'] = pet_fdir
                    tmp_dict['bin_width'] = bb
                    tmp_dict['organ_mask'] = 'breast tumor'
                    tmp_dict['process_name'] = 'GetImageFeature_PET_RecurredTumors'
                    tmp_dict['process_version'] = proc_version
                    tmp_dict['voxel_size_mm3'] = [tuple(the_img_spacingRCS)] * len(tmp_dict)
                    tmp_dict['breast_side'] = breast_side
                    tmp_dict['tumor_tag'] = tumor_tag
                    # print theimf.feature_output

                    pt_features_data = pt_features_data.append(tmp_dict, ignore_index=True)
                    # print pt_features_data

                # monitor the memory usage
                mem_use.append(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)

                pt_features_data['pt_id'] = pt_features_data['pt_id'].astype('int')
                dataoutfname = '{}/PET_features_data_v{}_{}_{}_RecurredTumor_{}.json'.format(save_dir,proc_version,pt_id,breast_side, tumor_tag)
                pt_features_data.to_json(dataoutfname)
                print('write json file: {}!'.format(dataoutfname))

    # save the memory log
    outfile = open('{}/PET_Features_v{}_memory_usage_recur_tumor.txt'.format(save_dir,proc_version), 'w')
    outfile.write('\n'.join([str(s) for s in mem_use]))
    outfile.close()
    print('write .txt file: {}!'.format(outfile))

    df_data_all = pd.DataFrame(lst_img_vox_size)
    df_data_all = df_data_all.ix[:, col_names]
    fname = '{}/PET_img_vox_size_mm_RecurredTumor.csv'.format(save_dir)
    df_data_all.to_csv(fname)
    print('write .csv file: {}!'.format(fname))
