# -*- coding: utf-8 -*-
"""
Created on 7/26/16 4:57 PM

@author: shuang Shih-ying Huang
@goal: re-organize the PETCT image files so it's organize by pt_id instead of accession number

"""

import pandas as pd
import os
import DataTool as dbt
import subprocess
import glob
import re
import shutil


# # Get a list of pt_id that MRI features have run on..
# rootdir =  '/data/francgrp1/clare_work/Data'
# alljsons = glob.glob('{}/her2_ImageFeatures/MRI_features_data*.json'.format(rootdir))
# df_all = pd.DataFrame()
# for jj in alljsons:
#     df = pd.read_json(jj)
#     pd_series = df.ix[0,['pt_id','pt_mrn','pt_accession_num']].astype('str')
#     pd_series['dce_fname'] = df['dce_series_dmi_fname'].unique().tolist()
#     df_all = df_all.append(pd_series,ignore_index=True)
#
#
# outfname = '{}/her2_ClinicalData/her2_MRI_ImgFeat_info.json'.format(rootdir)
# df_all.to_json(outfname)

# go through all PETCT images on lassen:/data/francgrp1/incoming
img_fdir = '/data/francgrp1/incoming'

# get pt_id vs anon_accession_num df
data_fdir = '/data/francgrp1/clare_work/her2_TallyLink/PtInfo_PETMRI'
fname = '{}/MRI_PET_data.json'.format(data_fdir)

# make sure DTYPE is specified (astype('str') didn't cast the accession number correctly (off by 1)
mri_pet_df = pd.read_json(fname,dtype={'pt_id':str,'Anon_Accession_num_PET':str,'Accession_num_PET':str,'anon_acc_num':str,'exam_dir_num':str,'pt_id':str,'pt_mrn':str})
print mri_pet_df[mri_pet_df['pt_id'] == '96']

# get the best_Ct and best_PET series number
db_host = 'lassen.radiology.ucsf.edu'
db_auth_dir = '/data/francgrp1/packages/MONGO/AUTH'
db_user_name, db_pw = dbt.get_db_auth_info(db_auth_dir)
db_name = 'her2'
db_collection = 'petct_studies'

mri_pet_anon_acc_num = mri_pet_df['anon_acc_num'].tolist()
qr = {'AccessionNumber': {'$in': mri_pet_anon_acc_num}}
qr_var_display = {'MRN': 1, 'AccessionNumber': 1, 'best_pet_series': 1, 'best_ct_series': 1}
qr_df = dbt.read_mongo(db=db_name,collection=db_collection,query=qr,query_var_select=qr_var_display,host=db_host,username=db_user_name,password=db_pw)

# # join the qr_df with mri_pet_df to save for later
# jdf = pd.merge(mri_pet_df,qr_df,how='left',left_on='Anon_Accession_num_PET',right_on='AccessionNumber')
#
# # save the updated BEST petct series number to json files
# outfname = '{}/MRI_PET_bestCTPET_data.json'.format(data_fdir)
# jdf.to_json(outfname)
#
#
# go through all the entry in the mri_pet_df to re-organize the PETCT image data
root_petctdir = '/data/francgrp1/breast_radiomics/her2/PETCT'
the_pt_id = mri_pet_df['pt_id'].unique().tolist()
# the_pt_id = [120,123,137,139,10,30,44,79]
for pid in the_pt_id:
    print 'pt_id: {}'.format(pid)
    exam_dir_num = mri_pet_df.ix[mri_pet_df['pt_id'] == pid,'exam_dir_num'].iloc[0]
    anon_acc_num = mri_pet_df.ix[mri_pet_df['pt_id'] == pid,'Anon_Accession_num_PET'].iloc[0]
    exam_dir_name = '{}/{}'.format(root_petctdir,exam_dir_num)
    accnum_dirname_dict = dbt.get_examdir_acc_nums(exam_dir_name)
    print exam_dir_name
    print accnum_dirname_dict

    if anon_acc_num in accnum_dirname_dict.keys():
        src_dirname = '{}/{}'.format(exam_dir_name,accnum_dirname_dict[anon_acc_num])
        print src_dirname
    else:
        print 'anon acc #: {} is not found in the exam dir {}'.format(anon_acc_num,exam_dir_name)
        print 'accnum_dirname_dict: {}'.format(accnum_dirname_dict)
        continue

    # # remove des_dirname if it exists
    # des_dirname = '{}/{:0>3s}'.format(root_petctdir, pid)
    # if os.path.exists(des_dirname):
    #     subprocess.call(['rm','-r',des_dirname]) #instead of using shutil.rmtree
    # else:
    #     print '{} DOES NOT EXIST!'.format(des_dirname)
    #
    # # try subprocess to run linux command
    # subprocess.call(['mkdir','-p',des_dirname])
    # subprocess.call(['cp','-r',src_dirname,des_dirname])
    #
    # # src_dirnamehutil.copytree didnt work with symlinks, instead it copies files in the linked directory...
    # shutil.copytree(src_dirname,des_dirname,symlinks=True)

# # deteremine the PETCT exam closet to MRI exam date (MRI_PET_data.json is missing the below MRNs)
# fname = '/data/francgrp1/breast_radiomics/her2/TABLES/closest_pets.json'
# df = pd.read_json(fname,dtype={'PatientID':str})
#
# test_mrn = ['08075719', '47454751', '32140564', '46540182', '42704318', '49520689']
#
# for tt in test_mrn:
#     print df.loc[['exam_dir'],int(tt)].values
#     print df.loc['mrn',int(tt)]

# # go through all BEST PET Data Link
# img_fdir = '/data/francgrp1/breast_radiomics/her2'
# best_pet_dir = '{}/BEST_PETS_LINKS'.format(img_fdir)
# anon_accession_num = []
# best_pet_series_num = []
# img_dir_name = []
# for root, dirs, files in os.walk(best_pet_dir):
#     check = re.search(r'(.+)group_(\d+)',root)
#     if check:
#         for ss in dirs:
#             strfind = re.search(r'(\d+)-(\d+)',ss)
#             if strfind:
#                 anon_accession_num.append(strfind.group(1))
#                 best_pet_series_num.append(strfind.group(2))
#                 img_dir_name.append('{}/{}'.format(root,ss))
#             else:
#                 print 'cannot determine accession number or series id!'
#
# best_pet_df = pd.DataFrame({'anon_accession_num': anon_accession_num,'best_pet_series_num': best_pet_series_num,'img_dir_name': img_dir_name})
#
# # JOIN/merge to determine anon_access_num vs pt_id
# jdf = pd.merge(mri_pet_df,best_pet_df,how='left',left_on='Anon_Accession_num_PET',right_on='anon_accession_num')
# print jdf.loc[pd.isnull(jdf['anon_accession_num']),['pt_id']]


