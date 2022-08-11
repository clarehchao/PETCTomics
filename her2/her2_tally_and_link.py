# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 18:23:32 2016

@author: rharnish
"""

#%%

her2_pet_link_dir = '/data/francgrp1/breast_radiomics/her2/PETCT'

#%%

acc_num_lut_dir   = '/data/francgrp1/SpyderWorkspace/acc_num_lut'
import sys
sys.path.append(acc_num_lut_dir)
import acc_num_lut

import numpy  as np
import pandas as pd
import os
import process_launcher
from glob import glob
import datetime

#%%

# get lut map
fmap,rmap = acc_num_lut.parseLUTFile(os.path.join(acc_num_lut_dir,'accNumLUT.properties'))

# get orig and anon acc_nums
her2_orig_acc_nums_path = os.path.join(acc_num_lut_dir,'her2_pet_acc_nums.txt')
her2_orig_acc_nums = np.loadtxt(her2_orig_acc_nums_path,dtype=str)
her2_anon_acc_nums = [fmap[an] for an in her2_orig_acc_nums]


#%% look for previous tally

tally_dir        = '/data/francgrp1/anon_acc_tallies'
previous_tallies = glob('{}/anon_acc_tally*csv'.format(tally_dir))
have_previous_tally = False
if len(previous_tallies) > 0:
    have_previous_tally = True
    newest_tally = previous_tallies[0]
    newest_time = datetime.datetime.fromtimestamp(os.stat(newest_tally).st_mtime)
    for t in previous_tallies:
        mod_time = datetime.datetime.fromtimestamp(os.stat(t).st_mtime) # This is a datetime.datetime object!  
        print mod_time
        if mod_time > newest_time:
            newest_time  = mod_time
            newest_tally = t
    print 'newest tally: {}'.format(newest_tally)

#%%

# load tally file as dict
#tally_file = '/data/francgrp1/anon_acc_tally_20160418.csv'
tally_file = newest_tally
cols = ['anon_acc_num', 'exam_dir']
tally_dict_records = pd.read_csv(tally_file, names = cols).to_dict('records')
tally_dict = {}
for r in tally_dict_records:
    tally_dict[str(r['anon_acc_num'])] = str(r['exam_dir'])

#%%

nums_we_have = set(her2_anon_acc_nums).intersection(set(tally_dict.keys()))

nums_we_want = set(her2_anon_acc_nums).difference(set(tally_dict.keys()))

#%%

for n in nums_we_have:
    filestore_exam_dir = tally_dict[n]
#    print n, filestore_exam_dir
    exam_dir_name = os.path.basename(filestore_exam_dir)
    pat_dir       = os.path.dirname(filestore_exam_dir)
    pat_dir_name  = os.path.basename(pat_dir)
    dest_dir      = os.path.join(her2_pet_link_dir,pat_dir_name)
    dest_path     = os.path.join(dest_dir,exam_dir_name)
    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)
    cmd = ['ln','-s',
           filestore_exam_dir,
           dest_path]
    if not os.path.islink(dest_path):
        print 'make link: {}'.format(cmd)
#        process_launcher.run_command(cmd)
    else:
        two = 1+1
