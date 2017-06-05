#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 5/22/17

@author: shuang
@goal: this package will include various methods for visualizing the data for a given study

"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import math

def featurelabel_redef(ss):
    if re.search(r'FOstats_(.+)', ss):
        return 'First-order Stats'
    elif re.search(r'ShapeSize_(.+)', ss):
        return 'Shape and Size'
    elif re.search(r'texture_(.+)', ss):
        return 'GLCM texture features'
    else:
        return ss

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

def featurename_redef_petmr(the_name):
    # check fo the string prefixed by '_'
    # return the_name if no match was found
    if re.search(r'texture_(.*)_avg_(.*)', the_name):
        return '_'.join(('TX', re.search(r'texture_(.*)_avg_(.*)', the_name).group(1),re.search(r'texture_(.*)_avg_(.*)', the_name).group(2)))
    elif re.search(r'FOstats_(.*)_(.*)', the_name):
        return '_'.join(('FOstats', re.search(r'FOstats_(.*)_(.*)', the_name).group(1),re.search(r'FOstats_(.*)_(.*)', the_name).group(2)))
    elif re.search(r'ShapeSize_(.*)_(.*)', the_name):
        return '_'.join(('SS',re.search(r'ShapeSize_(.*)_(.*)', the_name).group(1),re.search(r'ShapeSize_(.*)_(.*)', the_name).group(2)))
    else:
        the_name

def featurelabel_redef_petmr(ss):
    if re.search(r'FOstats_(.+)_(.+)', ss):
        return ' '.join(('First-order Stats', re.search(r'FOstats_(.+)_(.+)', ss).group(2).upper()))
    elif re.search(r'ShapeSize_(.+)_(.+)', ss):
        return ' '.join(('Shape and Size', re.search(r'ShapeSize_(.+)_(.+)', ss).group(2).upper()))
    elif re.search(r'texture_(.+)_(.+)', ss):
        return ' '.join(('GLCM texture features',re.search(r'texture_(.+)_(.+)', ss).group(2).upper()))
    else:
        return ss

def ClustermapPlot(the_df, varname, fig_fname, img_feature_names,featname_def_func,featlabel_def_func,idvar_color_pal_mode = 1, fig_title=None, var_title=None, var_color_pal=None):
    """

    :param the_df: the input dataframe
    :param varname: a list of variable name, list
    :param fig_title: figure title
    :param fig_fname: figure file name to save
    :param var_title: variable title that one would like to show on the figure (instead of varname)
    :param var_color_pal: a list of variable color palette for cluster plots, list (if more than one variables are defined)
    :return: n/a
    assume that 'cs_class' is the independent variable for cluster class and 'ptid_side' is the sample's unique IDs
    """

    # make sure there's no NANs
    df_sub = the_df.ix[:, ['cs_class', 'ptid_side'] + varname + img_feature_names]
    df_sub = df_sub.dropna()
    print 'ClustermapPlot: df_sub.shape = {}'.format(df_sub.shape)

    df_sub = df_sub.sort_values('cs_class')
    the_csclass_sorted_ptidside = df_sub['ptid_side'].tolist()

    # re-format the radiomic features into columnes of feature name and feature values
    df_flatten_data = pd.melt(df_sub, id_vars='ptid_side', value_vars=img_feature_names, var_name='Image Feature')

    # change to pivot_table instead of pivot (just need to make sure all data is unique so it's not doing aggregate!!
    df_pivot = df_flatten_data.pivot_table(index='Image Feature', columns='ptid_side', values='value')
    df_pivot_sorted = df_pivot[the_csclass_sorted_ptidside]

    all_var_names = varname + ['cs_class']
    all_var_vals = [df_sub[cc].tolist() for cc in all_var_names]
    mcol = pd.MultiIndex.from_tuples(zip(*all_var_vals), names=all_var_names)

    df_pivot_sorted.columns = mcol

    # re-arrange the index name by mri and pet
    df_final = df_pivot_sorted.loc[img_feature_names]

    # cluster colorbar
    img_feature_vars = df_final.index.tolist()

    # categorical label for the image feature vars
    img_feature_labels = [featlabel_def_func(ss) for ss in img_feature_vars]

    # re-assign index names
    img_feature_vars_new = [featname_def_func(ss) for ss in img_feature_vars]
    df_final.index = img_feature_vars_new

    if idvar_color_pal_mode == 1:
        # set up the label color palette for the radiomic features
        img_feature_pal = sns.light_palette('navy', n_colors=len(set(img_feature_labels)), reverse=True)
        img_feature_lut = dict(zip(map(str, list(set(img_feature_labels))), img_feature_pal))
        img_feature_colors = pd.Series(img_feature_labels, index=df_final.index).map(img_feature_lut)
        df_img_feature_colors = pd.DataFrame(dict(ImageFeature=img_feature_colors))
        df_img_feature_colors.columns = ['']
        img_feature_labels_order = img_feature_lut.keys()
    elif idvar_color_pal_mode == 2:
        img_label_unique = list(set(img_feature_labels))
        imf_label_modality = [re.search(r'(.+) (.+)', ss).group(2) for ss in img_label_unique if re.search(r'(.+) (.+)', ss)]
        df_tmp = pd.DataFrame({'imf_label': img_label_unique, 'modality': imf_label_modality})
        # find out the freq of the different type of image modalities
        modality_freq = df_tmp['modality'].value_counts().as_matrix()

        # colorbrew_lib = ['YlGnBu','YlOrBr','PuBu','PuRd','PuBuGn']
        colorbrew_lib = ['Greens','Oranges','Purples']
        img_feature_pal = []
        for ii in range(len(modality_freq)):
            imf_pal_tmp = sns.color_palette(colorbrew_lib[ii], modality_freq[ii])
            img_feature_pal = img_feature_pal + imf_pal_tmp
        df_tmp.sort_values('modality',inplace=True)
        idx_lst = df_tmp.index.tolist()
        a,b = idx_lst.index(1), idx_lst.index(5)
        idx_lst[b], idx_lst[a] = idx_lst[a], idx_lst[b]
        imf_label_order = df_tmp.loc[idx_lst]

        img_feature_lut = dict(zip(imf_label_order['imf_label'].tolist(), img_feature_pal))
        img_feature_colors = pd.Series(img_feature_labels, index=df_final.index).map(img_feature_lut)
        df_img_feature_colors = pd.DataFrame(dict(ImageFeature=img_feature_colors))
        df_img_feature_colors.columns = ['']
        img_feature_labels_order = imf_label_order['imf_label'].tolist()


    # categorize the clinical data into sub-group for visualization (set col_colors for the clustermap)
    clinical_data_colors = pd.DataFrame()
    df_pivot_col_unique = df_final.columns.unique()

    # use this dict to save the unique_var_labels for displaying legend
    var_labels_dict = {}
    var_lut_dict = {}
    for ii, vv in enumerate(varname):
        if var_color_pal:
            var_pal = var_color_pal[ii]
        else:
            var_pal = sns.color_palette('Blues', len(df_pivot_col_unique))

        var_labels = [df_pivot_col_unique[i][ii] for i in range(len(df_pivot_col_unique))]
        unique_var_labels = list(set(var_labels))
        var_lut = dict(zip(unique_var_labels, var_pal))

        # save var_lables and var_lut for legend display
        var_labels_dict[vv] = unique_var_labels
        var_lut_dict[vv] = var_lut

        if var_title:
            clinical_data_colors[var_title[ii]] = pd.Series(var_labels, index=df_pivot_col_unique).map(var_lut)
        else:
            clinical_data_colors[varname[ii]] = pd.Series(var_labels, index=df_pivot_col_unique).map(var_lut)

    csclass_labels = [df_pivot_col_unique[i][ii + 1] for i in range(len(df_pivot_col_unique))]
    unique_csclass_labels = list(set(csclass_labels))
    csclass_pal = sns.light_palette('gray', len(unique_csclass_labels))
    csclass_lut = dict(zip(unique_csclass_labels, csclass_pal))
    clinical_data_colors['CS_class'] = pd.Series(csclass_labels, index=df_pivot_col_unique).map(csclass_lut)

    # Create a custom colormap for the heatmap values
    cmap = sns.diverging_palette(h_neg=10, h_pos=220, s=80, n=7, as_cmap=True)

    # go through different methods and metrics for plotting patients separated into triple neg and non triple-negative
    sns.set()
    sns.set(font="monospace")

    g = sns.clustermap(df_final, col_cluster=False, row_cluster=False,
                       col_colors=clinical_data_colors, row_colors=df_img_feature_colors,
                       z_score=0, cmap=cmap, linewidths=0, xticklabels=False, yticklabels=True,
                       figsize=(15, 15), vmin=-5., vmax=5.)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    # display the legend for the image feature colormap by adding a empty bar plot but display the legend
    for label in img_feature_labels_order:
        g.ax_row_dendrogram.bar(0, 0, color=img_feature_lut[label], label=label, linewidth=0)
    g.ax_row_dendrogram.legend(bbox_to_anchor=(1.4, 1.2), loc='best', title='Image Features')

    for label in unique_csclass_labels:
        g.ax_col_dendrogram.bar(0, 0, color=csclass_lut[label], label=label, linewidth=0)
    g.ax_col_dendrogram.legend(bbox_to_anchor=(0.85, 0.56), loc='center', title='Cluster Consensus Class')

    box_xx = 0.65
    box_yy = 1.23
    for ii, kk in enumerate(var_labels_dict.keys()):
        unique_var_labels = var_labels_dict[kk]
        var_lut = var_lut_dict[kk]
        hh_all = []
        for label in unique_var_labels:
            hh = g.ax_heatmap.bar(0, 0, color=var_lut[label], label=label, linewidth=0)
            hh_all.append(hh)
        lg = g.ax_heatmap.legend(handles=hh_all, bbox_to_anchor=(box_xx - 0.3 * ii, box_yy), loc='best', title=kk,
                                 ncol=int(math.ceil(len(unique_var_labels) * 0.5)))
        g.ax_heatmap.add_artist(lg)

    # position the heatmap colorbar appropriately
    g.cax.set_position([0.15, .2, .03, .45])
    g.cax.set_title('z-score')
    g.ax_heatmap.set(xlabel='', ylabel='')

    if fig_title:
        g.fig.suptitle(fig_title, fontweight='bold')

    # create the file-save directories
    g.savefig(fig_fname)

    plt.close()
    # plt.show()


def heatmap_plot(the_tab, row_labels, fig_name =''):
    """
    plot the data in the_df in a heatmap for paper and publication
    :param the_df: 
    :return: 
    """

    # TODO: plot a side bar to display feature cluster class
    # ==> might have to do a clustermap instead of a heatmap
    # sns.set()
    # sns.set(font="monospace")
    # ax = sns.heatmap(the_tab, annot=True, fmt='1.2f')
    # plt.setp(ax.get_yticklabels(), rotation=0)
    # plt.setp(ax.get_xticklabels(), rotation=90)

    sns.set()
    sns.set(font="monospace")

    # define the row colormap
    row_labels_pal = sns.light_palette('navy', n_colors=len(set(row_labels)), reverse=True)
    row_labels_lut = dict(zip(map(str, list(set(row_labels))), row_labels_pal))
    row_labels_colors = pd.Series(row_labels, index= the_tab.index).map(row_labels_lut)
    df_row_labels_colors = pd.DataFrame(dict(feature=row_labels_colors))
    df_row_labels_colors.columns = ['']


    g = sns.clustermap(the_tab, col_cluster=False, row_cluster=False, row_colors = df_row_labels_colors,linewidths=0, xticklabels=True, yticklabels=True,figsize=(15, 15))
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)

    if fig_name:
        g.savefig(fig_name)
        print 'save heatmap plot to {}'.format(fig_name)
        plt.close()
    else:
        plt.show()

