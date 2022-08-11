# -*- coding: utf-8 -*-
"""
Created on 10/11/16 3:52 PM

@author: shuang Shih-ying Huang
@goal: cast the data type for the nrrd files for the PET tumor mask to a given datatype for ease of I/O for radiomics processing

"""

import ImageProcess as ImP
import glob



rootdir = '/data/francgrp1'
mask_dir = '{}/breast_radiomics/her2/ALL_TUMORS/mevis_manual_segmentation'.format(rootdir)
all_nrrd_fname = glob.glob('{}/*.nrrd'.format(mask_dir))

the_cast_dtype = 'unsigned char'
for ff in all_nrrd_fname:
    print 'cast nrrd image {}'.format(ff)
    ImP.NRRDImageCast(ff,the_cast_dtype)
