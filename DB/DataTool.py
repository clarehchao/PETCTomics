# -*- coding: utf-8 -*-
"""
Created on 8/1/16 10:31 AM

@author: shuang Shih-ying Huang
@goal: a collection of database or image-info related tools
       e.g. python plug-in for mongodb and how to get data from a query via javascript
       e.g. extract dicom header info for further data analysis

"""

import pandas as pd
import pymongo
import dicom
import getpass
import os
import re

# DICOM tags
tag_PatientID = dicom.tag.Tag(0x0010, 0x0020)
tag_AccessionNumber = dicom.tag.Tag(0x0008, 0x0050)
tag_StudyDate = dicom.tag.Tag(0x0008, 0X0020)

# code source: http://stackoverflow.com/questions/16249736/how-to-import-data-from-mongodb-to-pandas
def _connect_mongo(host, port, username, password, db):
    """ A util for making a connection to mongo """

    if username and password:
        mongo_uri = 'mongodb://{}:{}@{}:{}/{}'.format(username, password, host, port, db)
        conn = pymongo.MongoClient(mongo_uri)
    else:
        conn = pymongo.MongoClient(host, port)

    return conn[db]

def get_db_auth_info(db_auth_dir):
    db_user_name = getpass.getuser()
    db_pw_path = os.path.join(db_auth_dir, db_user_name + '.auth')
    f = open(db_pw_path)
    db_pw = f.readline().strip()
    f.close()
    return db_user_name,db_pw


def read_mongo(db, collection, query={}, query_var_select={},host='localhost', port=27017, username=None, password=None, no_id=True):
    """ Read from Mongo and Store into DataFrame """

    # Connect to MongoDB
    db = _connect_mongo(host=host, port=port, username=username, password=password, db=db)

    # Make a query to the specific DB and Collection
    if query or query_var_select:
        cursor = db[collection].find(query,query_var_select)
    else:
        cursor = db[collection].find()

    # Expand the cursor and construct the DataFrame
    df = pd.DataFrame(list(cursor))

    # Delete the _id
    if no_id:
        df = df.drop('_id',axis=1)

    return df

def get_examdir_acc_nums(the_exam_dir):
    output_dict = {}
    subdirs1 = os.listdir(the_exam_dir)
    for sb in subdirs1:
        subdirs2 = os.listdir('{}/{}'.format(the_exam_dir,sb))
        dcm_fname = '{}/{}/{}/1.dcm'.format(the_exam_dir, sb, subdirs2[0])
        dcm_info = dicom.read_file(dcm_fname)
        acc_num = str(dcm_info[tag_AccessionNumber].value)
        output_dict[acc_num] = sb
        # find_pet_dirs = [re.search(r'1(\d*)',ss).group() for ss in subdirs2 if re.search(r'1(\d*)',ss)]
        # if find_pet_dirs:
        #     dcm_fname = '{}/{}/{}/1.dcm'.format(the_exam_dir,sb,find_pet_dirs[0])
        #     dcm_info = dicom.read_file(dcm_fname)
        #     acc_num = str(dcm_info[tag_AccessionNumber].value)
        #     output_dict[acc_num] = sb
        # else:
        #     print 'CANNOT FIND any direction prefixed with 1xxx'

    return output_dict