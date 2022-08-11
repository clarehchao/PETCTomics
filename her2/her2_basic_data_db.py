# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 21:32:26 2016

@author: rharnish
"""

#%%

import os
import sys

import pandas as pd

from glob import glob
from pprint import pprint

import pymongo
import dicom

import machine_specific_paths


#%%

acc_num_lut_dir   = '/data/francgrp1/SpyderWorkspace/acc_num_lut'
sys.path.append(acc_num_lut_dir)
import acc_num_lut

#%%

paths = machine_specific_paths.MachineSpecificPaths()
db_dir   = paths.db_dir

tabular_data_dir    = '/data/francgrp1/breast_radiomics/her2/TABLES'
#marjan_table_path   = os.path.join(tabular_data_dir,'TumorMRIDCEvsPETSUV_Marjan.xlsx')
#marjan_table_path   = os.path.join(tabular_data_dir,'FirstPETData_Marjan.xls')
mrn_acc_table_path  = os.path.join(tabular_data_dir,'clarerequestPT2_11april.csv')

#%% Get connection to DB and "her2" Collection

def mongodb_connect(client_uri):
    try:
        client = pymongo.MongoClient(defaultURI,serverSelectionTimeoutMS=2000)
        # force connection on a request as the
        # connect=True parameter of MongoClient seems
        # to be useless here
        client.server_info()
        print "Connected to server {}".format(client_uri)
        return client
    except pymongo.errors.ServerSelectionTimeoutError as err:
        # do whatever you need
        print "Failed to connect to server {}".format(client_uri)
        print(err)
        return None

#defaultURI = 'mongodb://localhost:27017/'
defaultURI = 'mongodb://lassen.radiology.ucsf.edu:27017/'

'''
Try to launch mongodb -- if already running 
this shouldn't hurt anything
'''
#db_launch  = 'mongod --auth --dbpath ' + db_dir
db_launch  = 'mongod --dbpath ' + db_dir
os.system(db_launch)

client = mongodb_connect(defaultURI)
if client == None:
    print 'unable to launch or connect to DB'

# main database is her2    
db = client.her2

# get db credentials
import getpass
db_auth_dir  = '/data/francgrp1/packages/MONGO/AUTH'
db_user_name = getpass.getuser()
db_pw_path   = os.path.join(db_auth_dir,db_user_name+'.auth')
f = open(db_pw_path)
db_pw = f.readline().strip()
f.close()
db.authenticate(db_user_name, password=db_pw)

coll_exam_paths = db.exam_paths
c = coll_exam_paths.find()
recordCount = c.count()
print 'found {} records in db.exam_paths collection'.format(recordCount)

#%%

'''
TODO: package this stuff and put it somewhere it can be 
reused
'''

# DB keys
dbKey_AccessionNumber = 'AccessionNumber'
dbKey_PatientID       = 'PatientID'
dbKey_StudyDate       = 'StudyDate'

# DICOM tags
tag_PatientID       = dicom.tag.Tag(0x0010,0x0020)
tag_AccessionNumber = dicom.tag.Tag(0x0008,0x0050)
tag_StudyDate       = dicom.tag.Tag(0x0008,0X0020)


#%%

'''
Populate her2.mri_studies with some basic info
'''

# mri_studies collection
coll_mri_studies = db.mri_studies
c = coll_mri_studies.find()
recordCount = c.count()
print 'found {} records in db.mri_studies collection'.format(recordCount)

her2_dir = '/data/francgrp1/breast_radiomics/her2'
mri_dmis = glob(her2_dir + '/MRI/*/*/*/dce/*sagL_tp1.dmi')

for dmi in mri_dmis:
    
    # check if there is already a record for this exam_dir
    dce_dir  = os.path.dirname(dmi)
    exam_dir = os.path.dirname(dce_dir)
    c = coll_exam_paths.find({"exam_dir":exam_dir})
    recordCount = c.count()
    if recordCount > 0:
        continue
    
    print 'adding data from {}'.format(dmi)
    dcm_info = dicom.read_file(dmi)
    pat_id     = dcm_info[tag_PatientID].value
    acc_num    = dcm_info[tag_AccessionNumber].value
    study_date = dcm_info[tag_StudyDate].value
    coll_mri_studies.update_one({dbKey_AccessionNumber:acc_num},{'$set':{'TEST':True}}, upsert=True)
    coll_mri_studies.update_one({dbKey_AccessionNumber:acc_num},{'$set':{dbKey_PatientID:pat_id}}, upsert=True)
    coll_mri_studies.update_one({dbKey_AccessionNumber:acc_num},{'$set':{dbKey_StudyDate:study_date}}, upsert=True)
    
    coll_exam_paths.update_one({dbKey_AccessionNumber:acc_num},{'$set':{'exam_dir':exam_dir}}, upsert=True)    
    coll_exam_paths.update_one({dbKey_AccessionNumber:acc_num},{'$set':{'sampled_dcm':dmi}}, upsert=True)
    
#%%
    
#pet_ct_1s = glob(her2_dir + '/PETCT/???/*/*/1.dcm')
pet_ct_1s = glob(her2_dir + '/PETCT/OLD_LINKED_DIRS/*/*/*/1.dcm')

#%%

# petct_studies collection
coll_petct_studies = db.petct_studies
c = coll_petct_studies.find()
recordCount = c.count()
print 'found {} records in db.petct_studies collection'.format(recordCount)

visited_exam_dirs = {}
for dcm in pet_ct_1s:
    
    # check if there is already a record for this exam_dir
    series_dir  = os.path.dirname(dcm)
    exam_dir    = os.path.dirname(series_dir)
    if visited_exam_dirs.has_key(exam_dir):
        continue

    print 'adding data from {}'.format(dcm)
    dcm_info = dicom.read_file(dcm)
    pat_id     = dcm_info[tag_PatientID].value
    acc_num    = dcm_info[tag_AccessionNumber].value
    study_date = dcm_info[tag_StudyDate].value
#    dict_test = {'version':1,'b':0.98,'c':[-1,0,1,2]}
#    coll_petct_studies.update_one({dbKey_AccessionNumber:acc_num},{'$set':{'TEST':True}}, upsert=True)
#    coll_petct_studies.update_one({dbKey_AccessionNumber:acc_num},{'$set':{'dict_test':dict_test}}, upsert=True)
    coll_petct_studies.update_one({dbKey_AccessionNumber:acc_num},{'$set':{dbKey_PatientID:pat_id}}, upsert=True)
    coll_petct_studies.update_one({dbKey_AccessionNumber:acc_num},{'$set':{dbKey_StudyDate:study_date}}, upsert=True)
    coll_petct_studies.update_one({dbKey_AccessionNumber:acc_num},{'$set':{'exam_dir':exam_dir}}, upsert=True)

    coll_exam_paths.update_one({dbKey_AccessionNumber:acc_num},{'$set':{'exam_dir':exam_dir}}, upsert=True)    
    coll_exam_paths.update_one({dbKey_AccessionNumber:acc_num},{'$set':{'sampled_dcm':dcm}}, upsert=True)

    visited_exam_dirs[exam_dir] = 'visited'
    
#%%

# get anon acc num lut map
fmap,rmap = acc_num_lut.parseLUTFile(os.path.join(acc_num_lut_dir,'accNumLUT.properties'))


#%%
   
import datetime
def dicom_date_to_datetime(dicom_date_str):
    datetime_date = datetime.datetime.strptime(dicom_date_str, '%Y%m%d')
    datetime_str  = datetime_date.isoformat()
    return datetime_date, datetime_str

cutoff_date, cutoff_date_str = dicom_date_to_datetime('19971101')
s = coll_petct_studies.find({'StudyDate':{'$lte': cutoff_date_str}})
for rec in s:
    print '\n\n____________________________________'
    pprint(rec)
    
#%%

# load data frames from mongodb and csv file
df_petct   = pd.DataFrame(list(coll_petct_studies.find()))
df_mri     = pd.DataFrame(list(coll_mri_studies.find()))
df_mrn_acc = pd.read_csv(mrn_acc_table_path)

# get original acc nums into df_petct and change StudyDate to show it has been anonymized
df_petct.rename(columns={'AccessionNumber':'anon_acc_num'}, inplace=True)
df_petct.rename(columns={'StudyDate':'anon_study_date'},    inplace=True)
df_petct['acc_num'] = df_petct['anon_acc_num'].apply(lambda x: rmap[x])
df_petct['acc_num'] = df_petct['acc_num'].astype(str)

# rename column to acc_num and ensure string type
df_mrn_acc.rename(columns={'accession':'acc_num'}, inplace=True)
df_mrn_acc['acc_num'] = df_mrn_acc['acc_num'].astype(str)
df_mrn_acc['mrn'] = df_mrn_acc['mrn'].astype(str)

# rename columns to acc_num and mrn and ensure string type
df_mri.rename(columns={'AccessionNumber':'acc_num'}, inplace=True)
df_mri.rename(columns={'PatientID':'mrn'}, inplace=True)
df_mri['acc_num'] = df_mri['acc_num'].astype(str)

# use df_mrn_acc to add mrn column to df_petct
df_merged        = pd.merge(df_mrn_acc,df_petct,on='acc_num',how='right')
df_merged['mrn'] = df_merged['mrn'].apply(lambda x: '{0:0>8}'.format(x))


#%%

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

years    = mdates.YearLocator()   # every year
months   = mdates.MonthLocator()  # every month
yearsFmt = mdates.DateFormatter('%Y')

#%%

import time

date_incr = 3210

non_matching = []

dict_mri = df_mri.to_dict('records')
dict_pet = df_merged.to_dict('records')

closest_pets = {}
for rec in dict_mri:
    
    exam_dates = []
    exam_types = []
    
    print '\n\n{}'.format('--------------'*6)
    mrn                = rec['mrn']
    if pd.isnull(mrn):
        print 'NO MRN!!!'
        continue
    mri_sd_dcm         = rec['StudyDate']
    mri_sd, mri_sd_str = dicom_date_to_datetime(mri_sd_dcm)
    exam_dates.append(mri_sd)
    exam_types.append(1)
    print 'MRI mrn: {} study date: {}'.format(mrn,mri_sd_str)
    

    dict_matching_pet = df_merged.loc[df_merged['mrn']==mrn].to_dict('records')
    smallest_abs_date_diff = 365*100 # 100 yrs.
    if len(dict_matching_pet) > 0:
        for d in dict_matching_pet:
            print
            acc_num      = d['acc_num']
            asd_dcm      = d['anon_study_date']
            asd, asd_str = dicom_date_to_datetime(asd_dcm)
            pet_sd           = asd + datetime.timedelta(days=date_incr)
            print 'acc_num {}:  {} --> {}'.format(acc_num,asd,pet_sd)
            petb4mr = pet_sd <= mri_sd
            date_diff = (pet_sd - mri_sd).days
            if date_diff < smallest_abs_date_diff:
                smallest_abs_date_diff = date_diff
                closest_pet_dict = d
                closest_pet_dict['study_date'] = pet_sd
            if abs(date_diff)<30:
                print 'REAL CLOSE DATES!!!'
                
            print 'pet before mr?: {}'.format(petb4mr)
            print 'date offset: {}'.format(date_diff)
            exam_dates.append(pet_sd)
            exam_types.append(2)
        closest_pets[mrn] = closest_pet_dict
        found_match = True
    else:
        print 'no matching pets'
        found_match = False
        
    fig, ax = plt.subplots(figsize=(12, 2))
    ax.scatter(exam_dates, exam_types)
    
    
    # format the ticks
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(yearsFmt)
    ax.xaxis.set_minor_locator(months)
    
    datemin = datetime.date( min(exam_dates).year - 1 ,  1 ,  1 )
    datemax = datetime.date( max(exam_dates).year + 1 , 12 , 31 )
#    datemin = datetime.date(2000, 1, 1)
#    datemax = datetime.date(2016, 1, 1)
    ax.set_xlim(datemin, datemax)
    
    
    # format the coords message box
    def price(x):
        return '$%1.2f' % x
    ax.format_xdata = mdates.DateFormatter('%Y-%m-%d')
    ax.format_ydata = price
    ax.grid(True)
    
    # rotates and right aligns the x labels, and moves the bottom of the
    # axes up to make room for them
    fig.autofmt_xdate()
    
    plt.show()  
    
    if found_match == False:
        print '**************'
        print mrn
        non_matching.append(mrn)
        time.sleep(10)

print 'non matching: {}'.format(non_matching)
        
#%%
        
missing = ['08075719', 
           '47454751', 
           '32140564', 
           '46540182', 
           '42704318', 
           '49520689']
        
for m in missing:
    print '\n\n\n {}'.format(m)
    tmp = df_merged.loc[df_merged['mrn']==m].to_dict('records')
    pprint(tmp)
    
#%%

the_good = []
the_bad  = []

for k in closest_pets.keys():
    study_dir       = None
    best_pet_series = None
    print '\n\n------------------------------------------'
    print k
    closest_pet = closest_pets[k]
    anon_acc_num = closest_pet['anon_acc_num']

    c = coll_exam_paths.find({'AccessionNumber':anon_acc_num})
    for rec in c:
        print '\nexam_paths'
        pprint(rec)
        if rec.has_key('exam_dir'):
            study_dir = rec['exam_dir']
            print 'STUDY_DIR from exam_paths: {}'.format(study_dir)  
            
    d = coll_petct_studies.find({'AccessionNumber':anon_acc_num})    
    for rec in d:
        print '\npetct_studies'
        pprint(rec)
        if rec.has_key('best_pet_series'):
            best_pet_series = rec['best_pet_series']
            print 'BEST_PET_SERIES from petct_studies: {}'.format(best_pet_series)
    
    if not study_dir:   
        print '\n\n\n\nNO STUDY DIR!!!!!!!!!'
        the_bad.append(closest_pet)
        continue
    
    if not best_pet_series:   
        print '\n\n\n\nNO BEST PET SERIES!!!!!!!!!'
        the_bad.append(closest_pet)
        continue
    
    the_good.append(closest_pet)
    

#%%
'''
save closest_pets dict to json
get rid of _id key since it was 
causing problems 
'''

closest_pets_path = '/data/francgrp1/breast_radiomics/her2/TABLES/closest_pets.json'
for k in closest_pets.keys():
    closest_pets[k].pop('_id', None)

df_closest_pets = pd.DataFrame(closest_pets)
df_closest_pets.to_json(closest_pets_path)

#%%

##import process_launcher
#
#best_pets_links_dir = '/data/francgrp1/breast_radiomics/her2/BEST_PETS_LINKS'
#
#n_series_to_group = 5
#
#idx = -1
#group_incr = 0    
#for SD in the_good:
#    idx += 1
#    if idx % n_series_to_group == 0:
#        group_incr += 1
#        print group_incr
#
#    print '\n\n___________________________________'
#    study_dir = SD['study_dir']
#    if pd.isnull(study_dir):
#        print 'NO STUDY DIR!!! continuing....'
#        continue
#    best_pet_series_number = str(int(SD['best_pet_series']))
#    pet_dir_to_link = os.path.join(study_dir,best_pet_series_number)
#    print pet_dir_to_link
#    
#    acc_num = SD['anon_acc_num']
#    uid = os.path.basename(study_dir)
#    dest_dir = os.path.join(best_pets_links_dir,'group_'+str(group_incr))
#    dest_path = os.path.join(dest_dir,acc_num+'-'+best_pet_series_number) 
#    if not os.path.exists(dest_dir):
#        os.makedirs(dest_dir)
#
#    cmd = ['ln','-s',
#           pet_dir_to_link,
#           dest_path]
#    if not os.path.islink(dest_path):
#        print 'MAKE LINK'
#        print 'cmd: {}'.format(cmd)
##        process_launcher.run_command(cmd)
#    else:
#        print 'DO NOT MAKE LINK'
        
