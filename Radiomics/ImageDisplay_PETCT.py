# -*- coding: utf-8 -*-
"""
Created on 9/13/16 3:00 PM

@author: shuang Shih-ying Huang
@goal:

"""

import ImageProcess as ImP
import ITKImageHelper
import image_geometry
import dicom_series
import glob
import re
import os
import numpy as np
import sys

if __name__ == '__main__':
    pid_start = int(sys.argv[1])
    pid_end = int(sys.argv[2])

    # # Lassen
    rootdir = '/data/francgrp1'
    imdir = '{}/breast_radiomics/her2/ALL_TUMORS'.format(rootdir)
    all_dir = [d for d in glob.glob('{}/*'.format(imdir)) if os.path.isdir(d)]
    all_pt_id = [int(re.search(r'{}/(\d+)'.format(imdir), ss).group(1)) for ss in all_dir if
                 re.search(r'{}/(\d+)'.format(imdir), ss)]
    mask_dir = '{}/mevis_manual_segmentation'.format(imdir)

    problem_id = [15,16,17,18,22,25,27,30,50,62,79,83,87,89,105,109,122,136,139,140]
    inconsistent_series = [16,17,25,27,30,87]

    # the_ids = sorted([dd for dd in all_pt_id if dd >= pid_start and dd <= pid_end])
    # the_ids = [15,25,27,30]
    pet_img_samplingRCS = []
    dcimg_diff_size_id = []
    the_ids = [17]
    for id in the_ids:
        print 'patient ID: {}'.format(id)

        # PETCT image dir
        image_fdirs = [d for d in glob.glob('{}/{:0>3d}/*/*'.format(imdir,id)) if os.path.isdir(d)]

        for image_fdir in image_fdirs:
            print image_fdir
            check = re.search(r'(.+)/images_(\w+)/(.+)', image_fdir)
            if check:
                breast_side = check.group(2)
                print 'breast side: {}'.format(breast_side)
            else:
                print 'Cannot determine the breast side of the PET images'
                continue

            # Tumor mask .nrrd
            check_mask_fname = glob.glob('{}/{:0>3d}*sag{}*_uc_origin.nrrd'.format(mask_dir, id, breast_side))
            # check_mask_fname = glob.glob('{}/{:0>3d}*sag{}*_unsigned_char.nrrd'.format(mask_dir, id, breast_side))
            if check_mask_fname:
                mask_fname = check_mask_fname[0]
            else:
                print 'Cannot find the mask .nrrd file with the breast side!'
                continue

            petct_img_series = dicom_series.read_files(image_fdir)
            petct_ds_list = [ss for ss in petct_img_series if len(ss.shape) > 2]
            if len(petct_ds_list) > 1:
                print 'Warning! more than ONE dicom_series is found...{}'.format(petct_ds_list)
                # find the appropriate dicom_series to look at
                dc_lens = [ss.shape[0] for ss in petct_ds_list]
                petct_ds = petct_ds_list[dc_lens.index(max(dc_lens))]
            else:
                petct_ds = petct_ds_list[0]

            print petct_ds

            ig_pet = image_geometry.ImageGeometry(petct_ds)
            pix_pet = petct_ds.get_pixel_array()
            petct_itk_img = ITKImageHelper.generate_oriented_itkImage(ig=ig_pet, pixarray=pix_pet)
            pet_img_samplingRCS.append(ig_pet.samplingRCS)

            # need to orient the PETCT to 'AXIAL' for correct mask and image registration
            petct_itk_img_orient = ITKImageHelper.itkImage_orient_to_axial(petct_itk_img)
            petct_array = ITKImageHelper.itkImage_to_ndarray(petct_itk_img_orient)

            mask_itk = ITKImageHelper.itkImage_read(mask_fname)
            # mask_itk_orient = ITKImageHelper.itkImage_orient_to_axial(mask_itk)
            ig_mask = image_geometry.ImageGeometry(mask_itk)
            mask_array = ITKImageHelper.itkImage_to_ndarray(mask_itk)
            if petct_array.shape != mask_array.shape:
                dcimg_diff_size_id.append('{}_{}'.format(id,breast_side))

            if petct_array.shape == mask_array.shape:
                print 'image and mask arrays have the same size!'
                # find the slice where tumor mask exists
                check = np.nonzero(mask_array)
                s1_idx = np.min(check[0]) - 2
                s2_idx = np.max(check[0]) + 2
                s1_idx = s1_idx if s1_idx >= 0 else 0
                s2_idx = s2_idx if s2_idx <= petct_array.shape[0] else petct_array.shape[0]

                ImP.display_overlay_volume(petct_array[s1_idx:s2_idx], mask_array[s1_idx:s2_idx],
                                           'PID {} side {}: PETCT image + Tumor mask'.format(id, breast_side),
                                           aspect=ig_pet.samplingRCS[0]/ig_pet.samplingRCS[1])
            else:
                print 'image and mask arrays have DIFFERENT size! O_O, pet img size: {} vs mask size: {}'.format(petct_array.shape,mask_array.shape)
                # retro-fit image geometry of mask to image geometry of the dicom series PET data
                ss, rr, cc = np.nonzero(mask_array)
                src_lst = zip(ss, rr, cc)
                mask_array_retrofit2pet = np.zeros(petct_array.shape)
                for src in src_lst:
                    world_coord = ig_mask.idx_to_coords(src[2], src[1], src[0])
                    pet_img_idx = ig_pet.coords_to_idx(world_coord[0], world_coord[1], world_coord[2])
                    print 'world_coord: {}, pet_img_idx: {}'.format(world_coord, pet_img_idx)
                    mask_array_retrofit2pet[pet_img_idx[2], pet_img_idx[0], pet_img_idx[1]] = 1

                # plot the retro-fit-to-PET mask with the PET dicom series data
                ImP.display_overlay_volume(petct_array, mask_array_retrofit2pet, 'PID {} side {}, Retrofit PET image + Tumor mask'.format(id,breast_side),
                                           aspect=ig_pet.samplingRCS[0] / ig_pet.samplingRCS[1])

    # print dcimg_diff_size_id
    # pet_img_samplingRCS_array = np.array(pet_img_samplingRCS)
    # np.save('/data/francgrp1/clare_work/tmp/pet_img_samplingRCS.npy', pet_img_samplingRCS_array)
    # print pet_img_samplingRCS_array
    # print np.amin(pet_img_samplingRCS_array, axis=0)
    # print np.amax(pet_img_samplingRCS_array, axis=0)




                # # display the MRI and MAMMIPET
    # # img_fname = '{}/MAMMIPET/PT001/DCE601.nrrd'.format(rootdir)
    # img_fname = '{}/MAMMIPET/PT001/IM-0001-0359.dcm'.format(rootdir)
    # itk_img = ITKImageHelper.itkImage_read(img_fname)
    #
    # # for image display
    # IG_img = image_geometry.ImageGeometry(itk_img)
    # ar = IG_img.samplingSRC[0] / IG_img.samplingSRC[1]
    # img_array = ITKImageHelper.itkImage_to_ndarray(itk_img)
    # ImP.display_volume(img_array,fig_title='MAMMIPET MRI images',aspect=ar)

    # # TODO: MAMMIPET dcm requires some major modification in dicom_series and image_geoemtry since the dicom tags aren't at the same place
    # image_fdir = '/data/francgrp1/MAMMIPET/PT001/PET'
    # petct_img_series = dicom_series.read_files(image_fdir)
    # petct_ds = petct_img_series[0]
    # print petct_ds.info
    # print petct_ds.description
    # ig = image_geometry.ImageGeometry(petct_ds)
    # pix = petct_ds.get_pixel_array()
    # ar = ig.samplingRCS[0]/ig.samplingRCS[1]
    # petct_itk_img = ITKImageHelper.generate_oriented_itkImage(ig=ig, pixarray=pix)
    #
    # # need to orient the PETCT to 'AXIAL' for correct mask and image registration
    # petct_itk_img_orient = ITKImageHelper.itkImage_orient_to_axial(petct_itk_img)
    # petct_array = ITKImageHelper.itkImage_to_ndarray(petct_itk_img_orient)
    #
    # mask_fname = '/data/francgrp1/MAMMIPET/PT001/TumorMask/IM-0001-0359_mevis_manual_seg_1_20.nrrd'
    # mask_itk = ITKImageHelper.itkImage_read(mask_fname)
    # mask_array = ITKImageHelper.itkImage_to_ndarray(mask_itk)
    #
    # # find the slice where tumor mask exists
    # check = np.nonzero(mask_array)
    # s1_idx = np.min(check[0])
    # s2_idx = np.max(check[0])
    #
    # ImP.display_overlay_volume(petct_array[s1_idx:s2_idx],mask_array[s1_idx:s2_idx], 'PETCT image + Tumor mask', aspect=ar)

