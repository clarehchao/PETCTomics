# -*- coding: utf-8 -*-
"""
Created on 6/21/16 5:42 PM

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
# import numpy as np
# import matplotlib.pyplot as plt


# mac computer
rootdir = '/Users/shuang/Documents/Proj_Radiomics/IM'

# # Lassen
# rootdir = '/data/francgrp1'

# find patient_id in a list of INT (not string)
# need to make sure all sub-dir is a directory (not files, as it will pick up duplicate patient id..)
imdir = '{}/breast_radiomics/her2/MRI'.format(rootdir)
all_dir = [d for d in glob.glob('{}/*'.format(imdir)) if os.path.isdir(d)]
all_pt_id = [int(re.search(r'{}/(\d+)'.format(imdir), ss).group(1)) for ss in all_dir if re.search(r'{}/(\d+)'.format(imdir), ss)]

# Go through all the patient data via .dmi
ids = [ii for ii in all_pt_id if ii != 88]  # don't have pt 88 json file for now...

# define the texture features list to compuate
feature_list = ['autocorrelation','cluster_prominence','cluster_shade','cluster_tendency','contrast', 'correlation','diff_entropy','dissimilarity', 'energy','entropy', 'homogeneity1','homogeneity2', 'idmn',
                'idn', 'inv_var','maxprob','sum_avg','sum_entropy','sum_var']
                
ids = [12]
pt_features_data = pd.DataFrame()
for pt_id in ids:
    print 'patient id: {}'.format(pt_id)

    # find the MRI VOI json file
    info_dir = '{}/MORE_INFO'.format(imdir)

    # VOI coordinate info and the time-sorted dmil files (brtool output)
    findjson = glob.glob('{}/{}_*.json'.format(info_dir, pt_id))
    if findjson:
        voi_lps_path = findjson[0]
    else:
        print 'VOI JSON file was NOT found!'
        continue

    CENTER, HALFLENGTHS, time_sorted_IMNAMES = ImSeg.GetVOIInfo(voi_lps_path)

    # DMI MRI volumes (from brtool output)
    pt_series_dir = glob.glob('{}/{:0>3d}/*/*'.format(imdir, pt_id))[0]
    pt_dmi_list = ['{}/dce/{}'.format(pt_series_dir, ff) for ff in time_sorted_IMNAMES]
    test = [pt_dmi_list[1]]
    for ff in test:
        dceSeries = dmi.DMI(ff)
        mri_itk_img = ITKImageHelper.generate_oriented_itkImage(dceSeries)

        # figure out the patient MRN
        tag_mrn = dicom.tag.Tag(0x0010,0x0020)
        pt_mrn = dceSeries._DS[tag_mrn].value

        # the generated tumor mask from shuang method
        mask_fname = '{}/TumorMask_shuang.nrrd'.format(pt_series_dir)
        mask_itk_img = ITKImageHelper.itkImage_read(mask_fname)

        # define texture feature list
        theimf = ImF.ImageFeature(mri_itk_img,feature_list,mask_itk_img)
        # theimf = ImF.ImageFeature(mri_itk_img,feature_list)
        Rneighbor = 1
        theimf._compute_texture_features(Rneighbor)
        print 'texture feature complete!'
        theimf._compute_first_order_stats()
        print 'first order stats complete!'
        theimf._compute_shape_size_features()
        print 'shape size feature complete!'

        tmp_dict = theimf.feature_output
        tmp_dict['pt_id'] = pt_id
        tmp_dict['pt_mrn'] = pt_mrn
        tmp_dict['dce_series_dmi_fname'] = ff

        pt_features_data = pt_features_data.append(tmp_dict,ignore_index=True)

        del theimf

pt_features_data['pt_id'] = pt_features_data['pt_id'].astype('int')
# print pt_features_data['dce_series_dmi_fname']
# print pt_features_data['pt_id']

# save the dataframe to .csv
dataoutfname = '{}/clare_work/Data/pt_features_data.csv'.format(rootdir)
pt_features_data.to_csv(dataoutfname)

