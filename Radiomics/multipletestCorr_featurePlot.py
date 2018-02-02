#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 9/18/17

@author: shuang
@goal: apply multiple comparison correction to the statistical test for association between each variable to a given outcome

"""

import pandas as pd
import StatsTool as st
import VizPlot as vp
import numpy as np


im_dir = '/Users/shuang/Documents/Proj_Radiomics/data/her2/her2_Analysis/PETMRI/PETbinwidth0.1_MRItp2_binwidth5'
info_dict = {'spearmancorr': ('Spearman\'s\nCorr', ['Tumor_Grade','T_stage','N_stage','Overall_stage'],['Tumor Grade', 'T stage', 'N stage', 'Overall stage']),
             'corr': ('Mutiple Reg\nProportion of Variance',['DF_1yr','DF_2yr','DF_3yr','DF_4yr','DF_5yr','BC_subtype'], ['Disease free at {} yr'.format(ii) for ii in range(1,6)] + ['Breast cancer subtype'])}
fnames = ['{}/assoc_{}_all_v2.csv'.format(im_dir, k) for k in info_dict.keys()]
corr_method = 'fdr_bh'

# the_dfs = [pd.read_csv(ff) for ff in fnames]
# for ft, df in zip(filetags, the_dfs):
#     st.apply_mult_corr1(df, im_dir, ft)

# use default multiple test correction setting, e.g. fdr_bh with alpha=0.5
for ff in fnames:
    st.apply_mult_corr2(ff, corr_method=corr_method)

drop_feat_name = ['MRI_VOLUME_BLUE','MRI_VOLUME_GREEN','MRI_VOLUME_PURPLE','MRI_VOLUME_RED',
                  'MRI_VOLUME_TUMOR', 'MRI_VOLUME_TUM_BLU','MRI_VOLUME_WHITE','PET_SUV_max',
                  'MRI_PE2ROI_SER_MEAN','MRI_SERROI_PE2_PEAK']

corr_thresh = 0.20
rmeg2_thresh = 0.04

for k, (vv, foi, ftnames) in info_dict.items():
    corr_fname = '{}/assoc_{}_all_v2_pvalcorr.csv'.format(im_dir,k)
    df_assoc = pd.read_csv(corr_fname)

    # heatmap and tables of KW H test p-val
    corr_pval_colname = 'pval_corr_{}'.format(corr_method)
    # the_tab1 = df_assoc.pivot('feature', 'outcome', corr_pval_colname).ix[:,foi]
    tmp_tab = df_assoc.pivot('feature','outcome','corr_coeff').ix[:,foi]
    tmp_idx = tmp_tab.index.tolist()
    the_feats = [ss for ss in tmp_idx if ss not in drop_feat_name]
    the_tab2 = tmp_tab.loc[the_feats]
    the_tab2.rename(columns=dict(zip(foi, ftnames)), inplace=True)

    # rename the index to the desired index, but I believe VizPlot clustermap converts the index already
    # idx = the_tab2.index.tolist()
    # idx_new = [vp.featurename_redef_petmr(ss) for ss in idx]
    # idx_lut = dict(zip(idx, idx_new))
    # the_tab2.rename(index=idx_lut, inplace=True)

    # use correlation coeff instead, since the pval with multiple test correction wipe out lots potential variables
    fig_name = '{}/featureVSoutcome_{}_withmask.pdf'.format(im_dir, k)
    if k == 'spearmancorr':
        tmp1 = the_tab2.abs()
        indx_title = 'feature / Spearman\'s rank coefficient'
        thresh = corr_thresh
    else:
        tmp1 = the_tab2
        indx_title = 'feature / Multiple regression proportion of variance'
        thresh = rmeg2_thresh

    the_mask = tmp1 < thresh
    tmp2 = the_tab2[tmp1 >= thresh]

    vminmax = [np.nanmin(tmp2.values), np.nanmax(tmp2.values)]

    the_tab3 = vp.clustermap_plot_simple(the_tab2, idvar_color_pal_mode=2, row_label_title='Radiomic Features', vminmax=vminmax,mask=the_mask, fig_name=fig_name, value_title=vv, annot=True, fmt='.2f')
    # vp.clustermap_plot_simple(the_tab1, idvar_color_pal_mode=2, row_label_title='Radiomic Features', vminmax=vminmax,mask=the_mask, fig_name=fig_name, value_title='Kruskal-Wallis Test\nP-value')
    # # vp.heatmap_plot(the_tab1, the_mask=the_mask, vminmax=[0,0.05], fig_name=fig_name)

    # save as tables for supplementary data for manuscript
    tab_name = '{}/table_FeatsVS{}.csv'.format(im_dir, k)
    the_tab3.index.name = indx_title
    the_tab3.to_csv(tab_name)

# # heatmap of corr coeff
# the_tab2 = df_assoc.pivot('feature', 'outcome', 'corr_coeff')
# the_tab2 = the_tab2.ix[idx]
# the_tab2.columns = [col_label_lut[ss] for ss in the_tab2.columns]
# the_tab2_oi = the_tab2.ix[:,imf_order]
# if corr_type == 1:
#     fname_tab2 = '{}/corrcoeff_regr_featuresVSoutcome.csv'.format(im_dir)
# else:
#     fname_tab2 = '{}/spearman_corr_featuresVSoutcome.csv'.format(im_dir)
# the_tab2_oi.to_csv(fname_tab2)
#
#
# # the mask: true, will mask data, false: will NOT mask data
# if corr_type == 1:
#     fig_name = '{}/featureVSoutcome_corrcoeff_Ncluster{}_withmask.pdf'.format(cc_dir, the_Ncluster)
# else:
#     fig_name = '{}/featureVSoutcome_spearmancorr_Ncluster{}_withmask.pdf'.format(cc_dir, the_Ncluster)
# # the_mask = the_tab2_oi < 0.3
# # the_mask = the_tab1_oi > 0.05
# cor_coeff_cutoff = 0.3
# the_mask = (the_tab1_oi > 0.05) | (abs(the_tab2_oi) < cor_coeff_cutoff)
# # vminmax = [df_assoc['corr_coeff'].min(), df_assoc['corr_coeff'].max()]
# vminmax = [cor_coeff_cutoff, df_assoc['corr_coeff'].max()]
#
# if corr_type == 1:
#     var_title = 'Multiple\nCorr. Coeff.\nfor Regression'
# else:
#     var_title = 'Spearman\'s Rank\nCorrelation'
# vp.clustermap_plot_simple(the_tab2_oi, row_labels, row_label_title=row_label_title, vminmax=vminmax, mask=the_mask, fig_name=fig_name, value_title=var_title)
#
#
