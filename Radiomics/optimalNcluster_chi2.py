# -*- coding: utf-8 -*-
"""
Created on 4/11/17 11:15 AM

@author: shuang Shih-ying Huang
@goal: determine the optimal N of clusters and the corresponding medican cluster consensus

"""

import pandas as pd
import numpy as np
import itertools as itt


# PET image radiomics
rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'
pet_dir = '{}/her2_Analysis/PET/IsoVoxel_IMGBIN128_GLCMBIN64'.format(rootdir)
fname = '{}/chi2_relevancetest_all.csv'.format(pet_dir)
chi2_all_df = pd.read_csv(fname)
fdf = chi2_all_df[chi2_all_df['N_mincluster'] > 5]
idx = fdf.groupby(['Ncluster','v2']).apply(lambda df: df.pval.argmin())

# determine the number of outcome variables
the_outcome_vars_list = chi2_all_df['v2'].unique().tolist()
print the_outcome_vars_list

for ov in the_outcome_vars_list:
    idx = pd.DataFrame(idx)
    df_tmp = idx.loc[(slice(None), [ov]), :]
    idx2 = df_tmp[0]
    tmp_df = fdf.ix[idx2,['cluster_method','cluster_linkage','dist_method','N_sample','Ncluster','v2','pval']]
    for ii in range(len(tmp_df.index)):
        if tmp_df.ix[tmp_df.index[ii],'cluster_linkage'] is np.nan:
            cc_dir = '{}/ConsensusCluster_{}_{}'.format(pet_dir,tmp_df.ix[tmp_df.index[ii],'cluster_method'],tmp_df.ix[tmp_df.index[ii],'dist_method'])
            # print 'NAN is found! {}'.format(cc_dir)
        else:
            cc_dir = '{}/ConsensusCluster_{}_{}_{}'.format(pet_dir, tmp_df.ix[tmp_df.index[ii], 'cluster_method'],tmp_df.ix[tmp_df.index[ii], 'cluster_linkage'],tmp_df.ix[tmp_df.index[ii], 'dist_method'])
            # print 'NAN is NOT found! {}'.format(cc_dir)
        cc_fname = '{}/ClusterConsensus.csv'.format(cc_dir)
        cc_df = pd.read_csv(cc_fname)
        avg_cc_ss = cc_df.groupby('k')['clusterConsensus'].median()
        tmp_df.ix[tmp_df.index[ii],'median_clusterconsensus'] = avg_cc_ss[int(tmp_df.ix[tmp_df.index[ii],'Ncluster'])]
    print tmp_df
#
# # combine N-stage and overall stage for PET data
# the_cm = 'hc'
# the_cl = 'average'
# the_dd = 'spearman'
# the_v2 = 'Overall_stage'
# the_ncluster = 3
# print fdf.ix[(fdf['cluster_method'] == the_cm) & (fdf['cluster_linkage'] == the_cl) & (fdf['dist_method'] == the_dd) & (fdf['v2'] == the_v2) & (fdf['Ncluster'] == the_ncluster)]
#
# cc_dir = '{}/ConsensusCluster_{}_{}_{}'.format(pet_dir, the_cm, the_cl, the_dd)
# cc_fname = '{}/ClusterConsensus.csv'.format(cc_dir)
# cc_df = pd.read_csv(cc_fname)
# avg_cc_ss = cc_df.groupby('k')['clusterConsensus'].median()
# print 'median cluster consensus [N_cluster = {}] = {}'.format(the_ncluster,avg_cc_ss[the_ncluster])


# MRI image radiomics -- take 1
# theTP = [1, 2, 3]
# theglcmBin = [64, 128, 256]
# img_bin = 256
# rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'
#
# chi2_all_df = pd.DataFrame()
# for tt, bb in itt.product(theTP, theglcmBin):
#     # print 'TP: {}, glcm_bin: {}'.format(tt,bb)
#     mri_dir = '{}/her2_Analysis/MRI/IsoVoxel_TP{}_GLCMBIN{}'.format(rootdir,tt,bb)
#     fname = '{}/chi2_relevancetest_all.csv'.format(mri_dir)
#     chi2_df = pd.read_csv(fname)
#     # print chi2_df.shape
#     chi2_df['TP'] = tt
#     chi2_df['glcm_bin'] = bb
#     chi2_df['img_bin'] = img_bin
#     chi2_all_df = chi2_all_df.append(chi2_df,ignore_index=True)
#
# print chi2_all_df
# fdf = chi2_all_df[(chi2_all_df['N_mincluster'] > 5) & (chi2_all_df['TP'] == 2) & (chi2_all_df['glcm_bin'] == 128)]

# # MRI image radiomics -- take 2
# mri_tp = 2
# mri_glcmbin = 128
# rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'
# mri_dir = '{}/her2_Analysis/MRI/IsoVoxel_TP{}_GLCMBIN{}'.format(rootdir, mri_tp, mri_glcmbin)
# fname = '{}/chi2_relevancetest_all.csv'.format(mri_dir)
# chi2_all_df = pd.read_csv(fname)
# N_cluster_min = 5
# fdf = chi2_all_df[(chi2_all_df['N_mincluster'] > N_cluster_min)]
#
# # check = chi2_all_df[(chi2_all_df['N_mincluster'] > N_cluster_min) & (chi2_all_df['pval'] < 0.05) & (chi2_all_df['v2'] == 'T_stage')]
# # print check
#
# idx = fdf.groupby(['Ncluster','v2']).apply(lambda df: df.pval.argmin())
#
# # determine the number of outcome variables
# the_outcome_vars_list = chi2_all_df['v2'].unique().tolist()
#
# for ov in the_outcome_vars_list:
#     idx = pd.DataFrame(idx)
#     df_tmp = idx.loc[(slice(None), [ov]), :]
#     idx2 = df_tmp[0]
#     # tmp_df = fdf.ix[idx2,['TP','glcm_bin','cluster_method','cluster_linkage','dist_method','Ncluster','v2','pval']]
#     tmp_df = fdf.ix[idx2, ['cluster_method', 'cluster_linkage', 'dist_method', 'Nsample','Ncluster', 'v2', 'pval']]
#     for ii in range(len(tmp_df.index)):
#         # mri_dir = '{}/her2_Analysis/MRI/IsoVoxel_TP{}_GLCMBIN{}'.format(rootdir, tmp_df.ix[tmp_df.index[ii], 'TP'],tmp_df.ix[tmp_df.index[ii], 'glcm_bin'])
#         mri_dir = '{}/her2_Analysis/MRI/IsoVoxel_TP{}_GLCMBIN{}'.format(rootdir, mri_tp, mri_glcmbin)
#         if tmp_df.ix[tmp_df.index[ii],'cluster_linkage'] is np.nan:
#             cc_dir = '{}/ConsensusCluster_{}_{}'.format(mri_dir,tmp_df.ix[tmp_df.index[ii],'cluster_method'],tmp_df.ix[tmp_df.index[ii],'dist_method'])
#             # print 'NAN is found! {}'.format(cc_dir)
#         else:
#             cc_dir = '{}/ConsensusCluster_{}_{}_{}'.format(mri_dir, tmp_df.ix[tmp_df.index[ii], 'cluster_method'],tmp_df.ix[tmp_df.index[ii], 'cluster_linkage'],tmp_df.ix[tmp_df.index[ii], 'dist_method'])
#             # print 'NAN is NOT found! {}'.format(cc_dir)
#         cc_fname = '{}/ClusterConsensus.csv'.format(cc_dir)
#         cc_df = pd.read_csv(cc_fname)
#         avg_cc_ss = cc_df.groupby('k')['clusterConsensus'].median()
#         tmp_df.ix[tmp_df.index[ii],'median_clusterconsensus'] = avg_cc_ss[int(tmp_df.ix[tmp_df.index[ii],'Ncluster'])]
#     print tmp_df
