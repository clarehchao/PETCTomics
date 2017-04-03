# -*- coding: utf-8 -*-
"""
Created on 1/24/17 3:52 PM

@author: shuang Shih-ying Huang
@goal: visualize the cluster consensus result from R with a given data set

"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import itertools as itt
import numpy as np
import math


def featurename_redef(the_name):
    # check fo the string prefixed by '_'
    # return the_name if no match was found
    if re.search(r'texture_(.*)_avg', the_name):
        return 'TX_' + re.search(r'texture_(.*)_avg', the_name).group(1)
    elif re.search(r'FOstats_(.*)', the_name):
        return 'FO_' + re.search(r'FOstats_(.*)', the_name).group(1)
    elif re.search(r'ShapeSize_(.*)', the_name):
        return 'SS_' + re.search(r'ShapeSize_(.*)', the_name).group(1)
    else:
        the_name

def ClustermapPlot_TumorStatus(the_df,fig_title,fig_fname):
    # make sure there's no NANs
    the_df = the_df.dropna()

    texture_cols = ['texture_autocorrelation', 'texture_cluster_prominence', 'texture_cluster_shade',
                    'texture_cluster_tendency', 'texture_contrast', 'texture_correlation', 'texture_diff_entropy',
                    'texture_dissimilarity',
                    'texture_energy', 'texture_entropy', 'texture_homogeneity1', 'texture_homogeneity2', 'texture_idmn',
                    'texture_idn', 'texture_inv_var', 'texture_maxprob', 'texture_sum_avg', 'texture_sum_entropy',
                    'texture_sum_var']

    the_df = the_df.sort_values('cs_class')
    the_csclass_sorted_ptidside = the_df['ptid_side'].tolist()

    # re-format the radiomic features into columnes of feature name and feature values
    img_feature_names = ['FOstats_energy', 'FOstats_entropy', 'FOstats_kurtosis', 'FOstats_mean', 'FOstats_min',
                         'FOstats_max', 'FOstats_skewness', 'FOstats_uniformity', 'FOstats_variance',
                         'ShapeSize_compactness1', 'ShapeSize_compactness2', 'ShapeSize_max_euc_dis',
                         'ShapeSize_spherical_disproportion', 'ShapeSize_sphericity', 'ShapeSize_surf_area_cm2',
                         'ShapeSize_surface2volratio', 'ShapeSize_vol_cm3'] + [ss + '_avg' for ss in texture_cols]

    df_flatten_data = pd.melt(the_df, id_vars='ptid_side', value_vars=img_feature_names, var_name='Image Feature')

    # change to pivot_table instead of pivot (just need to make sure all data is unique so it's not doing aggregate!!
    df_pivot = df_flatten_data.pivot_table(index='Image Feature', columns='ptid_side', values='value')
    df_pivot_sorted = df_pivot[the_csclass_sorted_ptidside]

    TripleNeg_status = the_df['TripleNeg'].tolist()
    tumor_grade = the_df['Sjoerd_Grade'].tolist()
    histology = the_df['Marjan_Histology'].tolist()
    cs_class = the_df['cs_class'].tolist()

    # build in tumor status for each primary_ID
    mcol = pd.MultiIndex.from_tuples(zip(histology, tumor_grade, TripleNeg_status, cs_class),
                                     names=['Histology', 'Grade', 'TripleNeg', 'CS_class'])
    df_pivot_sorted.columns = mcol

    # cluster colorbar
    img_feature_vars = df_pivot_sorted.index.tolist()

    # categorical label for the image feature vars
    img_feature_labels = []
    for ss in img_feature_vars:
        if re.search(r'FOstats_(.+)', ss):
            img_feature_labels.append('First-order Stats')
        elif re.search(r'ShapeSize_(.+)', ss):
            img_feature_labels.append('Shape and Size')
        elif re.search(r'texture_(.+)', ss):
            img_feature_labels.append('GLCM texture features')

    # re-assign index names
    img_feature_vars_new = [featurename_redef(ss) for ss in img_feature_vars]
    df_pivot_sorted.index = img_feature_vars_new

    img_feature_pal = sns.light_palette('navy', n_colors=len(set(img_feature_labels)), reverse=True)
    img_feature_lut = dict(zip(map(str, list(set(img_feature_labels))), img_feature_pal))
    img_feature_colors = pd.Series(img_feature_labels, index=df_pivot_sorted.index).map(img_feature_lut)
    df_img_feature_colors = pd.DataFrame(dict(ImageFeature=img_feature_colors))
    df_img_feature_colors.columns = ['']

    # categorize the clinical data into sub-group for visualization (set col_colors for the clustermap
    clinical_data_colors = pd.DataFrame()
    df_pivot_col_unique = df_pivot_sorted.columns.unique()

    colors = ['windows blue', 'amber', 'greyish', 'brick']
    histology_pal = sns.xkcd_palette(colors)
    histology_labels = [df_pivot_col_unique[i][0] for i in range(len(df_pivot_col_unique))]
    unique_histology_labels = list(set(histology_labels))
    histology_lut = dict(zip(unique_histology_labels, histology_pal))
    clinical_data_colors['Histology'] = pd.Series(histology_labels, index=df_pivot_col_unique).map(histology_lut)

    light_pal = sns.light_palette((210, 90, 60), input="husl")
    grade_pal = light_pal[1:len(light_pal):2]
    grade_labels = [df_pivot_col_unique[i][1] for i in range(len(df_pivot_col_unique))]
    unique_grade_labels = list(set(grade_labels))
    grade_lut = dict(zip(unique_grade_labels, grade_pal))
    clinical_data_colors['Grade'] = pd.Series(grade_labels, index=df_pivot_col_unique).map(grade_lut)

    the_pal = sns.color_palette('Set2', 10)
    tripleNeg_colors_labels = [df_pivot_col_unique[i][2] for i in range(len(df_pivot_col_unique))]
    clinical_data_colors['TripleNeg'] = pd.Series(tripleNeg_colors_labels, index=df_pivot_col_unique).map({0: (1., 1., 1.), 1: the_pal[2]})

    ch_pal = sns.color_palette('cubehelix',10)
    csclass_labels = [df_pivot_col_unique[i][3] for i in range(len(df_pivot_col_unique))]
    unique_csclass_labels = list(set(csclass_labels))
    Ndiv = int(round(len(ch_pal)/float(len(unique_csclass_labels))))
    csclass_pal = ch_pal[::Ndiv]
    csclass_lut = dict(zip(unique_csclass_labels,csclass_pal))
    clinical_data_colors['CS_class'] = pd.Series(csclass_labels, index=df_pivot_col_unique).map(csclass_lut)

    # Create a custom colormap for the heatmap values
    cmap = sns.diverging_palette(h_neg=10, h_pos=220, s=80, n=7, as_cmap=True)

    # go through different methods and metrics for plotting patients separated into triple neg and non triple-negative
    sns.set()
    sns.set(font="monospace")

    g = sns.clustermap(df_pivot_sorted, col_cluster=False, row_cluster=False,
                       col_colors=clinical_data_colors, row_colors=df_img_feature_colors,
                       z_score=0, cmap=cmap, linewidths=0, xticklabels=False, yticklabels=True,
                       figsize=(15, 15))
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    # display the legend for the image feature colormap by adding a empty bar plot but display the legend
    for label in list(set(img_feature_labels)):
        g.ax_row_dendrogram.bar(0, 0, color=img_feature_lut[label], label=label, linewidth=0)
    g.ax_row_dendrogram.legend(bbox_to_anchor=(1.4, 1.2), loc='best', title='Image Features')

    for label in unique_histology_labels:
        g.ax_col_dendrogram.bar(0, 0, color=histology_lut[label], label=label, linewidth=0)
    g.ax_col_dendrogram.legend(bbox_to_anchor=(0.85,0.56), loc='center', title='Tumor histology')

    for label in unique_grade_labels:
        g.ax_heatmap.bar(0, 0, color=grade_lut[label], label=label, linewidth=0)
    g.ax_heatmap.legend(bbox_to_anchor=(0.5, 1.2), loc='best', title='Tumor grade')

    for label in unique_csclass_labels:
        g.ax_col_colors.bar(0, 0, color=csclass_lut[label], label=label, linewidth=0)
    g.ax_col_colors.legend(bbox_to_anchor=(0.3, 3.0), loc='best', title='Cluster Consensus Class')

    # position the heatmap colorbar appropriately
    g.cax.set_position([0.15, .2, .03, .45])
    g.cax.set_title('z-score')
    g.ax_heatmap.set(xlabel='', ylabel='')

    g.fig.suptitle(fig_title, fontweight='bold')

    # create the file-save directories
    g.savefig(fig_fname)

    plt.close()
    # plt.show()

def ClustermapPlot_PtOutcome1(the_df,fig_title,fig_fname):
    # make sure there's no NANs
    the_df = the_df.dropna()

    texture_cols = ['texture_autocorrelation', 'texture_cluster_prominence', 'texture_cluster_shade',
                    'texture_cluster_tendency', 'texture_contrast', 'texture_correlation', 'texture_diff_entropy',
                    'texture_dissimilarity',
                    'texture_energy', 'texture_entropy', 'texture_homogeneity1', 'texture_homogeneity2', 'texture_idmn',
                    'texture_idn', 'texture_inv_var', 'texture_maxprob', 'texture_sum_avg', 'texture_sum_entropy',
                    'texture_sum_var']

    the_df = the_df.sort_values('cs_class')
    the_csclass_sorted_ptidside = the_df['ptid_side'].tolist()

    # re-format the radiomic features into columnes of feature name and feature values
    img_feature_names = ['FOstats_energy', 'FOstats_entropy', 'FOstats_kurtosis', 'FOstats_mean', 'FOstats_min',
                         'FOstats_max', 'FOstats_skewness', 'FOstats_uniformity', 'FOstats_variance',
                         'ShapeSize_compactness1', 'ShapeSize_compactness2', 'ShapeSize_max_euc_dis',
                         'ShapeSize_spherical_disproportion', 'ShapeSize_sphericity', 'ShapeSize_surf_area_cm2',
                         'ShapeSize_surface2volratio', 'ShapeSize_vol_cm3'] + [ss + '_avg' for ss in texture_cols]

    df_flatten_data = pd.melt(the_df, id_vars='ptid_side', value_vars=img_feature_names, var_name='Image Feature')

    # change to pivot_table instead of pivot (just need to make sure all data is unique so it's not doing aggregate!!
    df_pivot = df_flatten_data.pivot_table(index='Image Feature', columns='ptid_side', values='value')
    df_pivot_sorted = df_pivot[the_csclass_sorted_ptidside]


    recur_status = the_df['Recurrence_CR'].tolist()
    cs_class = the_df['cs_class'].tolist()

    # build in tumor status for each primary_ID
    mcol = pd.MultiIndex.from_tuples(zip(recur_status, cs_class), names=['Recurrence', 'CS_class'])
    df_pivot_sorted.columns = mcol

    # cluster colorbar
    img_feature_vars = df_pivot_sorted.index.tolist()

    # categorical label for the image feature vars
    img_feature_labels = []
    for ss in img_feature_vars:
        if re.search(r'FOstats_(.+)', ss):
            img_feature_labels.append('First-order Stats')
        elif re.search(r'ShapeSize_(.+)', ss):
            img_feature_labels.append('Shape and Size')
        elif re.search(r'texture_(.+)', ss):
            img_feature_labels.append('GLCM texture features')

    # re-assign index names
    img_feature_vars_new = [featurename_redef(ss) for ss in img_feature_vars]
    df_pivot_sorted.index = img_feature_vars_new

    img_feature_pal = sns.light_palette('navy', n_colors=len(set(img_feature_labels)), reverse=True)
    img_feature_lut = dict(zip(map(str, list(set(img_feature_labels))), img_feature_pal))
    img_feature_colors = pd.Series(img_feature_labels, index=df_pivot_sorted.index).map(img_feature_lut)
    df_img_feature_colors = pd.DataFrame(dict(ImageFeature=img_feature_colors))
    df_img_feature_colors.columns = ['']

    # categorize the clinical data into sub-group for visualization (set col_colors for the clustermap
    clinical_data_colors = pd.DataFrame()
    df_pivot_col_unique = df_pivot_sorted.columns.unique()

    the_pal = sns.color_palette('Set2', 10)
    recur_labels = [df_pivot_col_unique[i][0] for i in range(len(df_pivot_col_unique))]
    unique_recur_labels = list(set(recur_labels))
    recur_pal = [(1., 1., 1.),the_pal[3],the_pal[7]]
    recur_lut = dict(zip(unique_recur_labels,recur_pal))
    clinical_data_colors['Recurrence'] = pd.Series(recur_labels, index=df_pivot_col_unique).map(recur_lut)

    ch_pal = sns.color_palette('cubehelix',10)
    csclass_labels = [df_pivot_col_unique[i][1] for i in range(len(df_pivot_col_unique))]
    unique_csclass_labels = list(set(csclass_labels))
    Ndiv = int(round(len(ch_pal)/float(len(unique_csclass_labels))))
    csclass_pal = ch_pal[::Ndiv]
    csclass_lut = dict(zip(unique_csclass_labels,csclass_pal))
    clinical_data_colors['CS_class'] = pd.Series(csclass_labels, index=df_pivot_col_unique).map(csclass_lut)

    # Create a custom colormap for the heatmap values
    cmap = sns.diverging_palette(h_neg=10, h_pos=220, s=80, n=7, as_cmap=True)

    # go through different methods and metrics for plotting patients separated into triple neg and non triple-negative
    sns.set()
    sns.set(font="monospace")

    g = sns.clustermap(df_pivot_sorted, col_cluster=False, row_cluster=False,
                       col_colors=clinical_data_colors, row_colors=df_img_feature_colors,
                       z_score=0, cmap=cmap, linewidths=0, xticklabels=False, yticklabels=True,
                       figsize=(15, 15))
    plt.setp(g.ax_heatmap.get_yticklabels(),rotation=0)

    # display the legend for the image feature colormap by adding a empty bar plot but display the legend
    for label in list(set(img_feature_labels)):
        g.ax_row_dendrogram.bar(0, 0, color=img_feature_lut[label], label=label, linewidth=0)
    g.ax_row_dendrogram.legend(bbox_to_anchor=(1.4, 1.2), loc='best', title='Image Features')

    for label in unique_csclass_labels:
        g.ax_col_dendrogram.bar(0, 0, color= csclass_lut[label], label=label, linewidth=0)
    g.ax_col_dendrogram.legend(bbox_to_anchor=(0.85,0.56), loc='center', title='Cluster Consensus Class')

    for label in unique_recur_labels:
        g.ax_heatmap.bar(0, 0, color=recur_lut[label], label=label, linewidth=0)
    g.ax_heatmap.legend(bbox_to_anchor=(0.5, 1.2), loc='best', title='Recurrence')


    # position the heatmap colorbar appropriately
    g.cax.set_position([0.15, .2, .03, .45])
    g.cax.set_title('z-score')
    g.ax_heatmap.set(xlabel='', ylabel='')

    g.fig.suptitle(fig_title, fontweight='bold')

    # create the file-save directories
    g.savefig(fig_fname)

    plt.close()
    # plt.show()


def ClustermapPlot_PtOutcome2(the_df,fig_title,fig_fname):
    # make sure there's no NANs
    the_df = the_df.dropna()

    texture_cols = ['texture_autocorrelation', 'texture_cluster_prominence', 'texture_cluster_shade',
                    'texture_cluster_tendency', 'texture_contrast', 'texture_correlation', 'texture_diff_entropy',
                    'texture_dissimilarity',
                    'texture_energy', 'texture_entropy', 'texture_homogeneity1', 'texture_homogeneity2', 'texture_idmn',
                    'texture_idn', 'texture_inv_var', 'texture_maxprob', 'texture_sum_avg', 'texture_sum_entropy',
                    'texture_sum_var']

    the_df = the_df.sort_values('cs_class')
    the_csclass_sorted_ptidside = the_df['ptid_side'].tolist()

    # re-format the radiomic features into columnes of feature name and feature values
    img_feature_names = ['FOstats_energy', 'FOstats_entropy', 'FOstats_kurtosis', 'FOstats_mean', 'FOstats_min',
                         'FOstats_max', 'FOstats_skewness', 'FOstats_uniformity', 'FOstats_variance',
                         'ShapeSize_compactness1', 'ShapeSize_compactness2', 'ShapeSize_max_euc_dis',
                         'ShapeSize_spherical_disproportion', 'ShapeSize_sphericity', 'ShapeSize_surf_area_cm2',
                         'ShapeSize_surface2volratio', 'ShapeSize_vol_cm3'] + [ss + '_avg' for ss in texture_cols]

    df_flatten_data = pd.melt(the_df, id_vars='ptid_side', value_vars=img_feature_names, var_name='Image Feature')

    # change to pivot_table instead of pivot (just need to make sure all data is unique so it's not doing aggregate!!
    df_pivot = df_flatten_data.pivot_table(index='Image Feature', columns='ptid_side', values='value')
    df_pivot_sorted = df_pivot[the_csclass_sorted_ptidside]

    recur_type = the_df['Recurrence_Type'].tolist()
    cs_class = the_df['cs_class'].tolist()

    # build in tumor status for each primary_ID
    mcol = pd.MultiIndex.from_tuples(zip(recur_type, cs_class), names=['Recurrence type', 'CS_class'])
    df_pivot_sorted.columns = mcol

    # cluster colorbar
    img_feature_vars = df_pivot_sorted.index.tolist()

    # categorical label for the image feature vars
    img_feature_labels = []
    for ss in img_feature_vars:
        if re.search(r'FOstats_(.+)', ss):
            img_feature_labels.append('First-order Stats')
        elif re.search(r'ShapeSize_(.+)', ss):
            img_feature_labels.append('Shape and Size')
        elif re.search(r'texture_(.+)', ss):
            img_feature_labels.append('GLCM texture features')

    # re-assign index names
    img_feature_vars_new = [featurename_redef(ss) for ss in img_feature_vars]
    df_pivot_sorted.index = img_feature_vars_new

    img_feature_pal = sns.light_palette('navy', n_colors=len(set(img_feature_labels)), reverse=True)
    img_feature_lut = dict(zip(map(str, list(set(img_feature_labels))), img_feature_pal))
    img_feature_colors = pd.Series(img_feature_labels, index=df_pivot_sorted.index).map(img_feature_lut)
    df_img_feature_colors = pd.DataFrame(dict(ImageFeature=img_feature_colors))
    df_img_feature_colors.columns = ['']

    # categorize the clinical data into sub-group for visualization (set col_colors for the clustermap)
    clinical_data_colors = pd.DataFrame()
    df_pivot_col_unique = df_pivot_sorted.columns.unique()

    colors = ["windows blue", "amber", "greyish", "faded green", "dusty purple"]
    recurtype_pal = sns.xkcd_palette(colors)
    recurtype_labels = [df_pivot_col_unique[i][0] for i in range(len(df_pivot_col_unique))]
    unique_recurtype_labels = list(set(recurtype_labels))
    recurtype_lut = dict(zip(unique_recurtype_labels, recurtype_pal))
    clinical_data_colors['Recurrence_Type'] = pd.Series(recurtype_labels, index=df_pivot_col_unique).map(recurtype_lut)

    ch_pal = sns.color_palette('cubehelix',10)
    csclass_labels = [df_pivot_col_unique[i][1] for i in range(len(df_pivot_col_unique))]
    unique_csclass_labels = list(set(csclass_labels))
    Ndiv = int(round(len(ch_pal)/float(len(unique_csclass_labels))))
    csclass_pal = ch_pal[::Ndiv]
    csclass_lut = dict(zip(unique_csclass_labels,csclass_pal))
    clinical_data_colors['CS_class'] = pd.Series(csclass_labels, index=df_pivot_col_unique).map(csclass_lut)

    # Create a custom colormap for the heatmap values
    cmap = sns.diverging_palette(h_neg=10, h_pos=220, s=80, n=7, as_cmap=True)

    # go through different methods and metrics for plotting patients separated into triple neg and non triple-negative
    sns.set()
    sns.set(font="monospace")

    g = sns.clustermap(df_pivot_sorted, col_cluster=False, row_cluster=False,
                       col_colors=clinical_data_colors, row_colors=df_img_feature_colors,
                       z_score=0, cmap=cmap, linewidths=0, xticklabels=False, yticklabels=True,
                       figsize=(15, 15))
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    # display the legend for the image feature colormap by adding a empty bar plot but display the legend
    for label in list(set(img_feature_labels)):
        g.ax_row_dendrogram.bar(0, 0, color=img_feature_lut[label], label=label, linewidth=0)
    g.ax_row_dendrogram.legend(bbox_to_anchor=(1.4, 1.2), loc='best', title='Image Features')

    for label in unique_csclass_labels:
        g.ax_col_dendrogram.bar(0, 0, color=csclass_lut[label], label=label, linewidth=0)
    g.ax_col_dendrogram.legend(bbox_to_anchor=(0.85,0.56), loc='center', title='Cluster Consensus Class')

    for label in unique_recurtype_labels:
        g.ax_heatmap.bar(0, 0, color=recurtype_lut[label], label=label, linewidth=0)
    g.ax_heatmap.legend(bbox_to_anchor=(0.5, 1.23), loc='best', title='Recurrence Type')


    # position the heatmap colorbar appropriately
    g.cax.set_position([0.15, .2, .03, .45])
    g.cax.set_title('z-score')
    g.ax_heatmap.set(xlabel='', ylabel='')

    g.fig.suptitle(fig_title, fontweight='bold')

    # create the file-save directories
    g.savefig(fig_fname)

    plt.close()
    # plt.show()

def ClustermapPlot(the_df,varname,fig_title,fig_fname):
    # make sure there's no NANs
    the_df = the_df.dropna()

    texture_cols = ['texture_autocorrelation', 'texture_cluster_prominence', 'texture_cluster_shade',
                    'texture_cluster_tendency', 'texture_contrast', 'texture_correlation', 'texture_diff_entropy',
                    'texture_dissimilarity',
                    'texture_energy', 'texture_entropy', 'texture_homogeneity1', 'texture_homogeneity2', 'texture_idmn',
                    'texture_idn', 'texture_inv_var', 'texture_maxprob', 'texture_sum_avg', 'texture_sum_entropy',
                    'texture_sum_var']

    the_df = the_df.sort_values('cs_class')
    the_csclass_sorted_ptidside = the_df['ptid_side'].tolist()

    # re-format the radiomic features into columnes of feature name and feature values
    img_feature_names = ['FOstats_energy', 'FOstats_entropy', 'FOstats_kurtosis', 'FOstats_mean', 'FOstats_min',
                         'FOstats_max', 'FOstats_skewness', 'FOstats_uniformity', 'FOstats_variance',
                         'ShapeSize_compactness1', 'ShapeSize_compactness2', 'ShapeSize_max_euc_dis',
                         'ShapeSize_spherical_disproportion', 'ShapeSize_sphericity', 'ShapeSize_surf_area_cm2',
                         'ShapeSize_surface2volratio', 'ShapeSize_vol_cm3'] + [ss + '_avg' for ss in texture_cols]

    df_flatten_data = pd.melt(the_df, id_vars='ptid_side', value_vars=img_feature_names, var_name='Image Feature')

    # change to pivot_table instead of pivot (just need to make sure all data is unique so it's not doing aggregate!!
    df_pivot = df_flatten_data.pivot_table(index='Image Feature', columns='ptid_side', values='value')
    df_pivot_sorted = df_pivot[the_csclass_sorted_ptidside]

    var_type = the_df[varname].tolist()
    cs_class = the_df['cs_class'].tolist()

    # build in tumor status for each primary_ID
    mcol = pd.MultiIndex.from_tuples(zip(var_type, cs_class), names=[varname, 'CS_class'])
    df_pivot_sorted.columns = mcol

    # cluster colorbar
    img_feature_vars = df_pivot_sorted.index.tolist()

    # categorical label for the image feature vars
    img_feature_labels = []
    for ss in img_feature_vars:
        if re.search(r'FOstats_(.+)', ss):
            img_feature_labels.append('First-order Stats')
        elif re.search(r'ShapeSize_(.+)', ss):
            img_feature_labels.append('Shape and Size')
        elif re.search(r'texture_(.+)', ss):
            img_feature_labels.append('GLCM texture features')

    # re-assign index names
    img_feature_vars_new = [featurename_redef(ss) for ss in img_feature_vars]
    df_pivot_sorted.index = img_feature_vars_new

    img_feature_pal = sns.light_palette('navy', n_colors=len(set(img_feature_labels)), reverse=True)
    img_feature_lut = dict(zip(map(str, list(set(img_feature_labels))), img_feature_pal))
    img_feature_colors = pd.Series(img_feature_labels, index=df_pivot_sorted.index).map(img_feature_lut)
    df_img_feature_colors = pd.DataFrame(dict(ImageFeature=img_feature_colors))
    df_img_feature_colors.columns = ['']

    # categorize the clinical data into sub-group for visualization (set col_colors for the clustermap)
    clinical_data_colors = pd.DataFrame()
    df_pivot_col_unique = df_pivot_sorted.columns.unique()

    # colors = ["windows blue", "amber", "greyish", "faded green", "dusty purple"]
    # var_pal = sns.xkcd_palette(colors)
    var_pal = sns.color_palette('hls',len(df_pivot_col_unique))
    var_labels = [df_pivot_col_unique[i][0] for i in range(len(df_pivot_col_unique))]
    unique_var_labels = list(set(var_labels))
    var_lut = dict(zip(unique_var_labels, var_pal))
    clinical_data_colors[varname] = pd.Series(var_labels, index=df_pivot_col_unique).map(var_lut)

    ch_pal = sns.color_palette('cubehelix',10)
    csclass_labels = [df_pivot_col_unique[i][1] for i in range(len(df_pivot_col_unique))]
    unique_csclass_labels = list(set(csclass_labels))
    Ndiv = int(round(len(ch_pal)/float(len(unique_csclass_labels))))
    csclass_pal = ch_pal[::Ndiv]
    csclass_lut = dict(zip(unique_csclass_labels,csclass_pal))
    clinical_data_colors['CS_class'] = pd.Series(csclass_labels, index=df_pivot_col_unique).map(csclass_lut)

    # Create a custom colormap for the heatmap values
    cmap = sns.diverging_palette(h_neg=10, h_pos=220, s=80, n=7, as_cmap=True)

    # go through different methods and metrics for plotting patients separated into triple neg and non triple-negative
    sns.set()
    sns.set(font="monospace")

    g = sns.clustermap(df_pivot_sorted, col_cluster=False, row_cluster=False,
                       col_colors=clinical_data_colors, row_colors=df_img_feature_colors,
                       z_score=0, cmap=cmap, linewidths=0, xticklabels=False, yticklabels=True,
                       figsize=(15, 15))
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    # display the legend for the image feature colormap by adding a empty bar plot but display the legend
    for label in list(set(img_feature_labels)):
        g.ax_row_dendrogram.bar(0, 0, color=img_feature_lut[label], label=label, linewidth=0)
    g.ax_row_dendrogram.legend(bbox_to_anchor=(1.4, 1.2), loc='best', title='Image Features')

    for label in unique_csclass_labels:
        g.ax_col_dendrogram.bar(0, 0, color=csclass_lut[label], label=label, linewidth=0)
    g.ax_col_dendrogram.legend(bbox_to_anchor=(0.85,0.56), loc='center', title='Cluster Consensus Class')

    for label in unique_var_labels:
        g.ax_heatmap.bar(0, 0, color=var_lut[label], label=label, linewidth=0)
    g.ax_heatmap.legend(bbox_to_anchor=(0.65, 1.23), loc='best', title=varname,ncol=int(math.ceil(len(unique_var_labels)*0.5)))


    # position the heatmap colorbar appropriately
    g.cax.set_position([0.15, .2, .03, .45])
    g.cax.set_title('z-score')
    g.ax_heatmap.set(xlabel='', ylabel='')

    g.fig.suptitle(fig_title, fontweight='bold')

    # create the file-save directories
    g.savefig(fig_fname)

    plt.close()
    # plt.show()

if __name__ == '__main__':
    # gather outcome data
    rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'
    outcome_fname = '{}/her2_ClinicalData/her2_outcome.csv'.format(rootdir)
    her2_outcome_df = pd.read_csv(outcome_fname,dtype={'MRN': str,'Recurrence Type Summary': str})

    def recur_func(x):
        if x == 'NEVER DISEASE FREE':
            return 2
        elif x == 'NONE/DISEASE FREE' or x is np.nan:
            return 0
        else:
            return 1
    # Recurrence_CR (recurrence info based on Cancer Registry field 'Recurrence Type Summary'): -1, always have cancer, 0, no recurrence, and 1, recur
    her2_outcome_df['Recurrence_CR'] = her2_outcome_df['Recurrence Type Summary'].map(recur_func)

    # Recurrence type: 0: assume NaN or NONE/DISEASE FREE are no recurrence
    # 1: bone recur, 2: local recur, 3: distant/systemic recur
    # -1: unknown recur type AND never disease free
    def recur_type_func(x):
        if x is np.nan or x == 'NONE/DISEASE FREE':
            return 0
        elif x == 'DIST RECUR, BONE':
            return 1
        elif re.search(r'LOCAL RECUR[\w\s.]*',x):
            return 2
        elif re.search(r'DIST RECUR[\w\s.]*',x) and x != 'DIST RECUR, BONE':
            return 3
        else:
            return 4
    her2_outcome_df['Recurrence_Type'] = her2_outcome_df['Recurrence Type Summary'].map(recur_type_func)

    # catergorize tumor and lymph node stages
    def Tstage_func(x):
        tstage_list = ['T0', 'Tis', 'T1mic', 'T1', 'T1a', 'T1b', 'T1c', 'T2', 'T3', 'T4', 'T4a', 'T4b', 'T4c', 'T4d',
                       'TX', 'calc']
        tstage_dict = dict(zip([ss.lower() for ss in tstage_list], range(len(tstage_list))))
        tstage_dict['t2yp'] = 7  # this is the exception for T2, not sure what T2yp is but assume it's the same as T2
        if x is np.nan:
            # print 'x is nan'
            return x
        elif x.lower().replace(' ','') in tstage_dict.keys():
            # print 'in dict,{}'.format(tstage_dict[x.lower().replace(' ','')])
            return tstage_dict[x.lower().replace(' ','')]

    her2_outcome_df['T-stage'] = her2_outcome_df['T-stage at surgery'].map(Tstage_func)

    def Nstage_func(x):
        nstage_list = ['N0', 'N1', 'N2a', 'N2b', 'N3a', 'N3b', 'N3c', 'NX', 'N3a', 'calc']
        nstage_dict = dict(zip([ss.lower() for ss in nstage_list], range(len(nstage_list))))

        if x is np.nan:
            return x
        elif x.lower().replace(' ','') in nstage_dict.keys():
            return nstage_dict[x.lower().replace(' ','')]

    her2_outcome_df['N-stage'] = her2_outcome_df['N-stage at surgery'].map(Nstage_func)

    # # Recurrence field does not seem reliable...
    # her2_outcome_df['Recurrence'] = her2_outcome_df['Recurrence'].apply(lambda x: int(x) if x != '#VALUE!' and x is not np.nan else np.nan)
    # # disease status = 0 when time interval since diagnosis
    # her2_outcome_df['disease status'] = her2_outcome_df['Time interval since diagnosis'].apply(lambda x: 0 if x < 5.0 else 1)

    the_outcome_df = her2_outcome_df.loc[:, ['MRN', 'Recurrence_CR','Recurrence_Type','T-stage','N-stage']]

    # # PET image feature data
    # dist_method_list = ['euclidean','spearman','pearson']
    # cluster_method_dict = {'pam': None, 'km': None, 'kmdist': None,'hc':['centroid','median','average','mcquitty','complete']}
    # Ncluster_list = [2, 3, 4, 5]
    # theimgbin = [32, 64, 128]
    # theglcmbin = [32, 64, 128]
    #
    # # # for testing purpose
    # # dist_method_list = ['pearson']
    # # cluster_method_dict = {'hc': ['average']}
    # # Ncluster_list = [4]
    # # theimgbin = [32]
    # # theglcmbin = [32]
    #
    # for tt, bb in itt.product(theimgbin, theglcmbin):
    #     pet_dir = '{}/her2_Analysis/PET/IsoVoxel_IMGBIN{}_GLCMBIN{}'.format(rootdir,tt,bb)
    #     df_fname = '{}/PETdataAll_glcmNbin{}_normNbin{}.csv'.format(pet_dir,bb,tt)
    #     data_df = pd.read_csv(df_fname, dtype={'MRN': str})
    #     jdf = pd.merge(data_df, the_outcome_df, on='MRN')
    #     print 'read the file {}'.format(df_fname)
    #     print data_df.shape
    #     print jdf.shape
    #
    #     for ncluster in Ncluster_list:
    #         for k, val in cluster_method_dict.items():
    #             if val:
    #                 for ll,dd in itt.product(val,dist_method_list):
    #                     cc_dir = '{}/ConsensusCluster_{}_{}_{}'.format(pet_dir,k,ll,dd)
    #                     cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir,ncluster)
    #                     print 'read the file {}'.format(cs_class_fname)
    #                     cs_class_df = pd.read_csv(cs_class_fname)
    #                     cs_class_df.columns = ['ptid_side', 'cs_class']
    #
    #                     # combine the cs_class to the_df
    #                     the_df = pd.merge(jdf, cs_class_df, how='left', on='ptid_side')
    #                     fig_title = 'cluster method: {}, linkage: {}, distance method: {}'.format(k, ll, dd)
    #                     fig_fname = '{}/clustermap_Tstage_kmax{}.pdf'.format(cc_dir,ncluster)
    #                     ClustermapPlot(the_df,'T-stage',fig_title,fig_fname)
    #
    #                     fig_fname = '{}/clustermap_Nstage_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     ClustermapPlot(the_df, 'N-stage', fig_title, fig_fname)
    #
    #                     fig_fname = '{}/clustermap_tumorstatus_kmax{}.pdf'.format(cc_dir,ncluster)
    #                     ClustermapPlot_TumorStatus(the_df,fig_title,fig_fname)
    #
    #                     fig_fname = '{}/clustermap_recurstatus_kmax{}.pdf'.format(cc_dir,ncluster)
    #                     ClustermapPlot_PtOutcome1(the_df, fig_title, fig_fname)
    #
    #                     fig_fname = '{}/clustermap_recurtype_kmax{}.pdf'.format(cc_dir,ncluster)
    #                     ClustermapPlot_PtOutcome2(the_df, fig_title, fig_fname)
    #             else:
    #                 for dd in dist_method_list:
    #                     cc_dir = '{}/ConsensusCluster_{}_{}'.format(pet_dir,k,dd)
    #                     cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir, ncluster)
    #                     cs_class_df = pd.read_csv(cs_class_fname)
    #                     cs_class_df.columns = ['ptid_side', 'cs_class']
    #                     print 'read the file {}'.format(cs_class_fname)
    #
    #                     # combine the cs_class to the_df
    #                     the_df = pd.merge(jdf, cs_class_df, how='left', on='ptid_side')
    #                     print the_df.shape
    #
    #                     fig_title = 'cluster method: {}, distance method: {}'.format(k, dd)
    #                     fig_fname = '{}/clustermap_Tstage_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     ClustermapPlot(the_df, 'T-stage', fig_title, fig_fname)
    #
    #                     fig_fname = '{}/clustermap_Nstage_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     ClustermapPlot(the_df, 'N-stage', fig_title, fig_fname)
    #
    #                     fig_fname = '{}/clustermap_tumorstatus_kmax{}.pdf'.format(cc_dir,ncluster)
    #                     ClustermapPlot_TumorStatus(the_df, fig_title, fig_fname)
    #
    #                     fig_fname = '{}/clustermap_recurstatus_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     ClustermapPlot_PtOutcome1(the_df, fig_title, fig_fname)
    #
    #                     fig_fname = '{}/clustermap_recurtype_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     ClustermapPlot_PtOutcome2(the_df, fig_title, fig_fname)



    # MRI image feature data
    theTP = [1, 2, 3]
    theglcmBin = [64, 128, 256]
    Ncluster_list = [2,3,4,5]
    dist_method_list = ['euclidean','spearman','pearson']
    cluster_method_dict = {'pam': None, 'km': None, 'kmdist': None,'hc':['centroid','median','average','mcquitty','complete']}

    # # for testing purpose
    # theTP = [1]
    # theglcmBin = [64]
    # Ncluster_list = [2]
    # dist_method_list = ['euclidean']
    # cluster_method_dict = {'km': None}

    for tt, bb in itt.product(theTP, theglcmBin):
        mri_dir = '{}/her2_Analysis/MRI/IsoVoxel_TP{}_GLCMBIN{}'.format(rootdir, tt, bb)
        df_fname = '{}/MRIdataAll_tp{}_Nbin{}.csv'.format(mri_dir, tt, bb)
        data_df = pd.read_csv(df_fname, dtype={'PRIMARY_ID': str, 'MRN': str})
        jdf = pd.merge(data_df, the_outcome_df, on='MRN')
        print 'read the file {}'.format(df_fname)

        for ncluster in Ncluster_list:
            for k, val in cluster_method_dict.items():
                if val:
                    for ll,dd in itt.product(val,dist_method_list):
                        cc_dir = '{}/ConsensusCluster_{}_{}_{}'.format(mri_dir,k,ll,dd)
                        cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir,ncluster)
                        print 'read file {}'.format(cs_class_fname)
                        cs_class_df = pd.read_csv(cs_class_fname)
                        cs_class_df.columns = ['ptid_side', 'cs_class']

                        # combine the cs_class to the_df
                        the_df = pd.merge(jdf, cs_class_df, how='left', on='ptid_side')
                        fig_title = 'cluster method: {}, linkage: {}, distance method: {}'.format(k, ll, dd)
                        fig_fname = '{}/clustermap_Tstage_kmax{}.pdf'.format(cc_dir, ncluster)
                        ClustermapPlot(the_df, 'T-stage', fig_title, fig_fname)

                        fig_fname = '{}/clustermap_Nstage_kmax{}.pdf'.format(cc_dir, ncluster)
                        ClustermapPlot(the_df, 'N-stage', fig_title, fig_fname)

                        fig_fname = '{}/clustermap_tumorstatus_kmax{}.pdf'.format(cc_dir, ncluster)
                        ClustermapPlot_TumorStatus(the_df, fig_title, fig_fname)

                        fig_fname = '{}/clustermap_recurstatus_kmax{}.pdf'.format(cc_dir, ncluster)
                        ClustermapPlot_PtOutcome1(the_df, fig_title, fig_fname)

                        fig_fname = '{}/clustermap_recurtype_kmax{}.pdf'.format(cc_dir, ncluster)
                        ClustermapPlot_PtOutcome2(the_df, fig_title, fig_fname)

                else:
                    for dd in dist_method_list:
                        cc_dir = '{}/ConsensusCluster_{}_{}'.format(mri_dir,k,dd)
                        cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir, ncluster)
                        print 'read file {}'.format(cs_class_fname)
                        cs_class_df = pd.read_csv(cs_class_fname)
                        cs_class_df.columns = ['ptid_side', 'cs_class']

                        # combine the cs_class to the_df
                        the_df = pd.merge(jdf, cs_class_df, how='left', on='ptid_side')
                        fig_title = 'cluster method: {}, linkage: {}, distance method: {}'.format(k, ll, dd)
                        fig_fname = '{}/clustermap_Tstage_kmax{}.pdf'.format(cc_dir, ncluster)
                        ClustermapPlot(the_df, 'T-stage', fig_title, fig_fname)

                        fig_fname = '{}/clustermap_Nstage_kmax{}.pdf'.format(cc_dir, ncluster)
                        ClustermapPlot(the_df, 'N-stage', fig_title, fig_fname)

                        fig_fname = '{}/clustermap_tumorstatus_kmax{}.pdf'.format(cc_dir, ncluster)
                        ClustermapPlot_TumorStatus(the_df, fig_title, fig_fname)

                        fig_fname = '{}/clustermap_recurstatus_kmax{}.pdf'.format(cc_dir, ncluster)
                        ClustermapPlot_PtOutcome1(the_df, fig_title, fig_fname)

                        fig_fname = '{}/clustermap_recurtype_kmax{}.pdf'.format(cc_dir, ncluster)
                        ClustermapPlot_PtOutcome2(the_df, fig_title, fig_fname)


