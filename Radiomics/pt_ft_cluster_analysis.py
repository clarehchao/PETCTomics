#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 6/6/17

@author: shuang
@goal: combine the result from patient cluster and feature cluster

"""

import pandas as pd
import DataHelper as dh
import VizPlot as vp
from scipy.stats import zscore
import numpy as np
import re

rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'

# PET and MRI image setting
the_mri_tp = 2
the_mri_glcmbin = 128
the_pet_imgnorm_bin = 128
the_pet_glcmbin = 64
the_Ncluster = 3

# PATIENT cluster: get the consensus cluster setting for a given outcome with the best Chi2 result
im_dir = '{}/her2_Analysis/PETMRI/PETgbin{}_imgbin{}_MRItp{}_gbin{}'.format(rootdir,the_pet_glcmbin,the_pet_imgnorm_bin,the_mri_tp,the_mri_glcmbin)
chi2_fname = '{}/BestChi2_outcome_Ncluster{}.csv'.format(im_dir, the_Ncluster)
df_pt_bestchi2 = pd.read_csv(chi2_fname)

# FEATURE cluster
df_medCC_all, df_oi = dh.Median_Cluster_Consensus(im_dir, 'featureConsensusCluster')
the_cm, the_cl, the_dm = df_oi.ix[df_oi['k'] == the_Ncluster,['cluster_method','cluster_linkage','dist_method']].as_matrix()[0]
print 'feature cluster setting: {}, {}, {}'.format(the_cm, the_cl, the_dm)
if isinstance(the_cl, str):
    cc_dir = '{}/featureConsensusCluster_{}_{}_{}'.format(im_dir, the_cm, the_cl, the_dm)
else:
    cc_dir = '{}/featureConsensusCluster_{}_{}'.format(im_dir, the_cm, the_dm)
cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir, the_Ncluster)
df_cs = pd.read_csv(cs_class_fname)
df_cs.columns = ['feature', 'cs_class']
df_cs = df_cs.sort_values('cs_class')

assoc_fname = '{}/assoc_corr_test_result_all.csv'.format(im_dir)
df_assoc = pd.read_csv(assoc_fname)
jdf = pd.merge(df_cs, df_assoc, on='feature')
jdf_colnames = jdf.columns.tolist()
jdf_colnames[jdf_colnames.index('cs_class')] = 'feature_cs_class'
jdf.columns = jdf_colnames
oh1 = jdf['outcome'].unique().tolist()
oh2 = oh1[:]
oh2[oh2.index('Sjoerd_Grade')] = 'Tumor Grade'
oh2[oh2.index('Marjan_Histology')] = 'Tumor Histology'
dict_oc_rename = dict(zip(oh1,oh2))
jdf['outcome'] = jdf['outcome'].apply(lambda x: dict_oc_rename[x])

# get all the feature and clinical data
dtall_fname = '{}/data_all.csv'.format(im_dir)
df_data_all = pd.read_csv(dtall_fname)
df_data_all.rename(columns={'Sjoerd_Grade':'Tumor Grade', 'Marjan_Histology': 'Tumor Histology'},inplace=True)

# compute zscore for all the feature data
feat_var_name = [ss for ss in df_data_all.columns.tolist() if any([ss.startswith('FOstats_'), ss.startswith('ShapeSize_'), ss.startswith('texture_'), ss.startswith('MRI_'), ss.startswith('PET_')])]
df_data_all[feat_var_name] = df_data_all[feat_var_name].apply(zscore)

# perform metric calculation only for the features with significant chi2 result from patient clustering
df_ptchi2_outcome = df_pt_bestchi2[df_pt_bestchi2.pval < 0.05]
for ii in df_ptchi2_outcome.index.tolist():
    ov, cm, cl, dm, pval = df_ptchi2_outcome.ix[ii, ['outcome_var','cluster_method','cluster_linkage','dist_method','pval']]
    print 'outcome var: {}, patient cluster chi2 pval: {}'.format(ov, pval)
    # get the patient cluster result
    if isinstance(cl, str):
        cc_dir = '{}/patientConsensusCluster_{}_{}_{}'.format(im_dir, cm, cl, dm)
    else:
        cc_dir = '{}/patientConsensusCluster_{}_{}'.format(im_dir, cm, dm)

    pt_cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir, the_Ncluster)
    df_pt_cs_class = pd.read_csv(pt_cs_class_fname)
    df_pt_cs_class.columns = ['ptid_side', 'patient_cs_class']

    # join df_data_all with pt CS clustering result
    jdf_pt_cs = pd.merge(df_data_all, df_pt_cs_class, on='ptid_side')

    df_tmp = jdf[jdf['outcome'] == ov]

    # set all corr_coeff < 0.3 to NaN
    df_tmp.ix[df_tmp.corr_coeff < 0.3, 'corr_coeff'] = np.nan

    # set all pval > 0.05 to NaN
    df_tmp.ix[df_tmp.pval > 0.05, 'corr_coeff'] = np.nan
    rank = df_tmp.groupby('feature_cs_class')['corr_coeff'].rank(na_option='keep')

    # append the rank column to df_tmp appropriately
    df_tmp.ix[rank.index.tolist(), 'feature_cs_rank'] = rank
    df_tmp['feature_cs_rank'].fillna(0, inplace=True)

    df_info = df_tmp.ix[:,['feature','feature_cs_class','feature_cs_rank']]
    jdf_pt_cs = jdf_pt_cs.merge(jdf_pt_cs.apply(dh.func_feat_rank, axis=1, args=(df_info,)),
                      left_index=True, right_index=True)

    # groupby patient CS clustering result for featclass metric
    # TODO: rank according to correlation, maybe have a corr coeff cut off for ranking
    # i.e.: set the rank to 0's if corr coeff is < 0.3 or etc. (so it doesn't include all
    df1 = jdf_pt_cs.groupby('patient_cs_class')[['featCSclass_1', 'featCSclass_2', 'featCSclass_3']].mean()
    df2 = jdf_pt_cs.groupby('patient_cs_class')[['featCSclass_1', 'featCSclass_2', 'featCSclass_3']].std()

    # make another table with 'mean +- stdev' so one can plot with heatmap and annotate mean+-stdev
    # unpivot the tables first to combine mean + stdev
    df1_unpiv = df1.unstack().reset_index(name='value')
    df1_unpiv.rename(columns={'level_0': 'feature_cs_class', 'value': 'rank_metric_avg'}, inplace=True)
    df1_unpiv['feature_cs_class'] = df1_unpiv['feature_cs_class'].apply(lambda x: re.search('featCSclass_(\d+)', x).group(1))

    df2_unpiv = df2.unstack().reset_index(name='value')
    df2_unpiv.rename(columns={'level_0': 'feature_cs_class', 'value': 'rank_metric_std'}, inplace=True)
    df2_unpiv['feature_cs_class'] = df2_unpiv['feature_cs_class'].apply(lambda x: re.search('featCSclass_(\d+)', x).group(1))

    jdf1df2 = pd.merge(df1_unpiv, df2_unpiv, on=['feature_cs_class','patient_cs_class'])
    # jdf1df2['avg_std_str'] = jdf1df2.apply(lambda x: u'{:.1f} \u00B1 {:.1f}'.format(x['rank_metric_avg'], x['rank_metric_std']).encode('utf-8'), axis=1)
    jdf1df2['avg_std_str'] = jdf1df2.apply(lambda x: r'{:.1f} $\pm$ {:.1f}'.format(x['rank_metric_avg'], x['rank_metric_std']), axis=1)

    data_tab = jdf1df2.pivot(index='feature_cs_class', columns='patient_cs_class', values='rank_metric_avg')
    annot_tab = jdf1df2.pivot(index='feature_cs_class', columns='patient_cs_class', values='avg_std_str').as_matrix()
    fig_fname = '{}/pt_feat_rank_metric_{}.pdf'.format(im_dir, ov)
    vp.heatmap_plot(data_tab, annot_labels=annot_tab, fig_name=fig_fname)