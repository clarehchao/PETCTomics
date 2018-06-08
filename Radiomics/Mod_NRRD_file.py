#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 6/4/18

@author: shuang
@goal: modify image format of a nrrd file and save to another nrrd file

"""

import nrrd
import os
import dicom_series
import glob
import re
import sys

rootdir = '/Users/shuang/Documents/Proj_Radiomics'

imdir = '{}/Data/her2/her2_IM/breastCA_pet_data'.format(rootdir)
all_img_dirs = os.walk(imdir).next()[1]

for ptid_side in all_img_dirs:
    print(ptid_side)

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
    print(series_des)

    slicer_dir = '{}/{}/3dslicer'.format(imdir, ptid_side)
    all_files = glob.glob('{}/*.nrrd'.format(slicer_dir))
    tmp = [ff for ff in all_files if re.search('{}/([\w\s.]+){}.nrrd'.format(slicer_dir,re.escape(series_des)), ff)]
    if not tmp:
        print('cannot find any PET nrrd file!')
        continue
    else:
        pet_nrrd_fname = tmp[0]
        data, header = nrrd.read(pet_nrrd_fname)
        data1 = data.astype('float32') #cast double to float
        header['type'] = 'float'
        fdir = os.path.dirname(pet_nrrd_fname)
        fname = os.path.basename(pet_nrrd_fname)
        tmp_lst = list(os.path.splitext(fname))
        tmp_lst.insert(1,'_f32')
        fname_new = '{}/{}'.format(fdir, ''.join(tuple(tmp_lst)))
        nrrd.write(fname_new, data1, options=header)
        print('write new nrrd file: {}'.format(fname_new))

