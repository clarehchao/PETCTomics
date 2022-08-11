#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 6/20/18

@author: shuang
@goal: univariate feature analysis of Radiomic features to Oncotype Dx score (obtained from Elissa Price, UCSF)

"""

import pandas as pd
import VizPlot as vp
import numpy as np

rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'

# PET and MRI image setting
the_mri_tp = 2
the_pet_binwidth = 0.1
the_mri_binwidth = 5
im_dir = '{}/her2_Analysis/PETMRI/PETbinwidth{}_MRItp{}_binwidth{}'.format(rootdir,the_pet_binwidth, the_mri_tp, the_mri_binwidth)

corr_thresh = 0.50
rmeg2_thresh = 0.04

info_dict = {'spearmancorr': ('Spearman\'s\nCorr', ['ODx_score'],['Oncotype Dx Score'])}
drop_feat_name = ['MRI_VOLUME_BLUE','MRI_VOLUME_GREEN','MRI_VOLUME_PURPLE','MRI_VOLUME_RED',
                  'MRI_VOLUME_TUMOR', 'MRI_VOLUME_TUM_BLU','MRI_VOLUME_WHITE','PET_SUV_max',
                  'MRI_PE2ROI_SER_MEAN','MRI_SERROI_PE2_PEAK']

for k, (vv, foi, ftnames) in info_dict.items():
    print(k, vv, foi, ftnames)
    corr_fname = '{}/assoc_{}_oncotype.csv'.format(im_dir,k)
    df_assoc = pd.read_csv(corr_fname)

    tmp_tab = df_assoc.pivot('feature','outcome','corr_coeff').ix[:,foi]
    tmp_idx = tmp_tab.index.tolist()
    the_feats = [ss for ss in tmp_idx if ss not in drop_feat_name]
    the_tab2 = tmp_tab.loc[the_feats]
    the_tab2.rename(columns=dict(zip(foi, ftnames)), inplace=True)

    fig_name = '{}/featureVSoutcome_spearmancorr_OncotypeDx_withmask.pdf'.format(im_dir)
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

    df_tmp = the_tab2.reset_index()
    print(df_tmp['Oncotype Dx Score'].min(), df_tmp['Oncotype Dx Score'].max())
    print(df_tmp.ix[df_tmp['Oncotype Dx Score'].abs() >= 0.6, :])

    vminmax = [np.nanmin(tmp2.values), np.nanmax(tmp2.values)]

    # the_tab3 = vp.clustermap_plot_simple(the_tab2, idvar_color_pal_mode=2, row_label_title='Radiomic Features', vminmax=vminmax,mask=the_mask, fig_name=fig_name, value_title=vv, annot=True, fmt='.2f')
    the_tab3 = vp.clustermap_plot_simple(the_tab2, idvar_color_pal_mode=2, row_label_title='Radiomic Features',
                                         vminmax=vminmax, mask=the_mask, fig_name=fig_name, value_title=vv)









