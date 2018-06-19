#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 6/5/18

@author: shuang
@goal: compare the radiomics of the primary tumors vs recurred tumors of the same patient

"""


import numpy as np
import pandas as pd
import glob
import re
import VizPlot as vp


if __name__ == '__main__':
    # gather all radiomic data from primary and recurred tumors
    rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'
    savedir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2/her2_Analysis/PET_Recur_Tumor'

    # get radiomics of all the primary tumor data
    fname1 = '{}/her2_Analysis/PETMRI/PETbinwidth0.1_MRItp2_binwidth5/data_all.csv'.format(rootdir)
    df_prim_all = pd.read_csv(fname1)

    # print(df_prim_all.columns.tolist())

    # find all PET radiomics
    pat = re.compile('_pet')
    feat_names = [ss for ss in df_prim_all.columns.tolist() if re.search('([\w.]+)_pet', ss)]
    new_feat_names = [re.search('([\w.]+)_pet', ss).group(1) for ss in df_prim_all.columns.tolist() if
                      re.search('([\w.]+)_pet', ss)]
    newer_feat_names = [re.search('([\w.]+)_avg', ss).group(1) if re.search('([\w.]+)_avg', ss) else ss for ss in
                        new_feat_names]

    the_col_names = feat_names + ['ptid_side']
    df_prim = df_prim_all.loc[:, the_col_names]

    # change feature name
    col_dict = dict(zip(feat_names, newer_feat_names))
    df_prim.rename(col_dict, axis='columns', inplace=True)
    df_prim['tumor_type'] = 'Primary'
    # print(df_prim.columns.tolist())

    json_dir = '{}/her2_ImageFeatures/IsoVoxelSize'.format(rootdir)
    all_jsons = glob.glob('{}/*.json'.format(json_dir))

    df_recur = pd.DataFrame()
    for jj in all_jsons:
        df_tmp = pd.read_json(jj)
        df_recur = df_recur.append(df_tmp, ignore_index=True)
    df_recur['FOstats_min'] = df_recur['FOstats_minmax'].apply(lambda x: x[0])
    df_recur['FOstats_max'] = df_recur['FOstats_minmax'].apply(lambda x: x[1])
    df_recur.drop(columns=['FOstats_minmax'], inplace=True)

    # get the average of texture features

    pat = re.compile('texture_')
    texture_cols = [ss for ss in df_recur.columns.tolist() if pat.match(ss)]
    for tc in texture_cols:
        df_recur[tc + '_avg'] = df_recur[tc].apply(np.mean)
        df_recur.drop(tc, axis=1, inplace=True)
    df_recur['tumor_type'] = df_recur['tumor_tag'].map(lambda x: '_'.join(['Recur', x]))
    df_recur['ptid_side'] = df_recur[['pt_id', 'breast_side']].apply(lambda x: '{}_{}'.format(x[0], x[1]), axis=1)
    newer_feat_names = [re.search('([\w.]+)_avg', ss).group(1) if re.search('([\w.]+)_avg', ss) else ss for ss in
                        df_recur.columns.tolist()]
    col_dict = dict(zip(df_recur.columns.tolist(), newer_feat_names))
    df_recur.rename(col_dict, axis='columns', inplace=True)

    col_of_interest = df_prim.columns.tolist()
    df_recur_oi = df_recur.loc[:, col_of_interest]
    df_prim_oi = df_prim.loc[:, col_of_interest]

    # combine primary and recur tumor DFs
    df_all = pd.concat([df_prim_oi, df_recur_oi], ignore_index=True)
    # print(df_all)

    # the data ready for bokeh plot
    ptid_sides = list(df_recur_oi.ptid_side.unique())
    # print(ptid_sides)

    # PLOT THE ANNULAR WEDGE PLOT for all patients
    vp.Plot_Annular(df_all, ptid_sides, savedir)

    # PLOT THE PEARSON CORRELATION COEFFICIENT PLOT for all patients
    fname ='{}/Radiomics_primaryVSrecur_PearsonCorr.pdf'.format(savedir)
    vp.Plot_Radiomics_Corr(df_all, ptid_sides, fig_name=fname)




