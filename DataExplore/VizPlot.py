#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 5/22/17

@author: shuang
@goal: this package will include various methods for visualizing the data for a given study

"""

import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import re
import math
import numpy as np
from bokeh.plotting import figure, show, output_file

from bokeh.palettes import Spectral3, Spectral4, Spectral5, Spectral6, Spectral7
from collections import OrderedDict
from bokeh.io import export_png
# from bokeh.io import export_svgs





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
        # print 'no other match, {}'.format(the_name)
        return the_name

def featurelabel_redef_petmr(ss):
    petmr_pat = re.compile('PET_|MRI_')
    if re.search(r'FOstats_(.+)_(.+)', ss):
        return ' '.join(('First-order Stats', re.search(r'FOstats_(.+)_(.+)', ss).group(2).upper()))
    elif re.search(r'ShapeSize_(.+)_(.+)', ss):
        return ' '.join(('Shape and Size', re.search(r'ShapeSize_(.+)_(.+)', ss).group(2).upper()))
    elif re.search(r'texture_(.+)_(.+)', ss):
        return ' '.join(('GLCM texture features',re.search(r'texture_(.+)_(.+)', ss).group(2).upper()))
    elif re.match(petmr_pat, ss):
        return 'Modality-specific metrics'
    else:
        return ss

def img_modality_def(ss):
    if re.search('(.+)_pet', ss) or re.search('PET_(.+)', ss):
        return 'pet'
    elif re.search('(.+)_mri', ss) or re.search('MRI_(.+)', ss):
        return 'mri'

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
    print 'ClustermapPlot: varialbe name {}, df_sub.shape = {}'.format(varname, df_sub.shape)

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
    # print img_feature_labels

    # re-assign index names
    img_feature_vars_new = [featname_def_func(ss) for ss in img_feature_vars]
    # print img_feature_vars_new
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
        if 'Modality-specific metrics' in img_label_unique:
            img_label_unique.remove('Modality-specific metrics')
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
        a,b = idx_lst.index(1), idx_lst.index(len(idx_lst)-1)
        idx_lst[b], idx_lst[a] = idx_lst[a], idx_lst[b]
        imf_label_order = df_tmp.loc[idx_lst]
        df_tmp2 = pd.DataFrame([['Modality-specific metrics','other']], columns=['imf_label','modality'])
        imf_label_order = imf_label_order.append(df_tmp2,ignore_index=True)
        img_feature_pal = img_feature_pal + sns.color_palette(colorbrew_lib[ii+1], 1)

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
        var_labels = [df_pivot_col_unique[i][ii] for i in range(len(df_pivot_col_unique))]
        unique_var_labels = list(set(var_labels))

        if var_color_pal:
            var_pal = var_color_pal[ii]
        else:
            # make sure the default palette is made with the # of unique outcome variables
            var_pal = sns.color_palette('Blues', len(var_labels))
        var_lut = dict(zip(unique_var_labels, var_pal))


        # save var_lables and var_lut for legend display
        if var_title:
            var_labels_dict[var_title[ii]] = unique_var_labels
            var_lut_dict[var_title[ii]] = var_lut
            clinical_data_colors[var_title[ii]] = pd.Series(var_labels, index=df_pivot_col_unique).map(var_lut)
        else:
            var_labels_dict[vv] = unique_var_labels
            var_lut_dict[vv] = var_lut
            clinical_data_colors[vv] = pd.Series(var_labels, index=df_pivot_col_unique).map(var_lut)

    csclass_labels = [df_pivot_col_unique[i][ii + 1] for i in range(len(df_pivot_col_unique))]
    unique_csclass_labels = list(set(csclass_labels))
    csclass_pal = sns.light_palette('gray', len(unique_csclass_labels))
    csclass_lut = dict(zip(unique_csclass_labels, csclass_pal))
    clinical_data_colors['CS_class'] = pd.Series(csclass_labels, index=df_pivot_col_unique).map(csclass_lut)

    # Create a custom colormap for the heatmap values
    cmap = sns.diverging_palette(h_neg=10, h_pos=220, s=80, n=7, as_cmap=True)

    # go through different methods and metrics for plotting patients separated into triple neg and non triple-negative
    sns.set()
    sns.set(font='Arial')

    # g = sns.clustermap(df_final, col_cluster=False, row_cluster=False,
    #                    col_colors=clinical_data_colors, row_colors=df_img_feature_colors,
    #                    z_score=0, cmap=cmap, linewidths=0, xticklabels=False, yticklabels=True,
    #                    figsize=(15, 15), vmin=-5., vmax=5.)
    g = sns.clustermap(df_final, col_cluster=False, row_cluster=False,
                       col_colors=clinical_data_colors, row_colors=df_img_feature_colors,
                       z_score=0, cmap=cmap, linewidths=0, xticklabels=False, yticklabels=False,
                       figsize=(15, 15), vmin=-5., vmax=5.)

    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    # display the legend for the image feature colormap by adding a empty bar plot but display the legend
    for label in img_feature_labels_order:
        g.ax_row_dendrogram.bar(0, 0, width=0,color=img_feature_lut[label], label=label, linewidth=0)
    g.ax_row_dendrogram.legend(bbox_to_anchor=(1.4, 1.2), loc='best', title='Image Features')

    for label in unique_csclass_labels:
        g.ax_col_dendrogram.bar(0, 0, width=0, color=csclass_lut[label], label=label, linewidth=0)
    g.ax_col_dendrogram.legend(bbox_to_anchor=(0.80, 0.48), loc='center', title='Cluster Consensus Class', ncol = len(unique_csclass_labels))

    box_xx = 0.65
    box_yy = 1.18
    for ii, kk in enumerate(var_labels_dict.keys()):
        unique_var_labels = var_labels_dict[kk]
        var_lut = var_lut_dict[kk]
        hh_all = []
        for label in unique_var_labels:
            hh = g.ax_heatmap.bar(0, 0, width=0, color=var_lut[label], label=label, linewidth=0)
            hh_all.append(hh)
        g.ax_heatmap.legend(bbox_to_anchor=(box_xx - 0.3 * ii, box_yy), loc='best', title=kk,
                             ncol=int(math.ceil(len(unique_var_labels) * 0.5)))

    # position the heatmap colorbar appropriately
    # g.cax.set_position([0.15, .2, .03, .45])
    g.cax.set_position([0.92, .2, .03, .45])
    g.cax.set_title('z-score')
    g.ax_heatmap.set(xlabel='', ylabel='')

    if fig_title:
        g.fig.suptitle(fig_title, fontweight='bold')

    # create the file-save directories
    g.savefig(fig_fname)

    plt.close()
    # plt.show()


def clustermap_plot_simple(the_tab, idvar_color_pal_mode, annot=None, fmt='', row_labels=None, row_label_title='', vminmax=[], mask=None, fig_name ='', value_title=''):
    """
    plot the data in the_df in a clustermap plot for paper and publication for simple tables
    :param the_df: 
    :return: 
    """

    sns.set()
    sns.set(font='Arial')

    # define the row colormap
    if idvar_color_pal_mode == 1:
        row_labels_pal = sns.light_palette('navy', n_colors=len(set(row_labels)), reverse=True)
        row_labels_lut = dict(zip(list(set(row_labels)), row_labels_pal))
        row_labels_colors = pd.Series(row_labels, index= the_tab.index).map(row_labels_lut)
        df_row_labels_colors = pd.DataFrame(dict(feature=row_labels_colors))
        df_row_labels_colors.columns = ['']
        row_labels_order = row_labels_lut.keys()


    elif idvar_color_pal_mode == 2:
        row_vars = the_tab.index.tolist()

        # sort the row label first by image modality or other
        row_vars.sort(key = img_modality_def)
        the_tab2 = the_tab.loc[row_vars, :]
        the_mask2 = mask.loc[row_vars, :]

        # categorical label for the image feature vars
        row_labels = [featurelabel_redef_petmr(ss) for ss in row_vars]

        # re-assign index names
        row_vars_new = [featurename_redef_petmr(ss) for ss in row_vars]

        the_tab2.index = row_vars_new
        the_mask2.index = row_vars_new
        # print(the_tab2)


        row_labels_unique = list(set(row_labels))
        row_labels_unique.remove('Modality-specific metrics')
        row_labels_modality = [re.search(r'(.+) (.+)', ss).group(2) for ss in row_labels_unique if re.search(r'(.+) (.+)', ss)]
        df_tmp = pd.DataFrame({'imf_label': row_labels_unique, 'modality': row_labels_modality})
        # find out the freq of the different type of image modalities
        modality_freq = df_tmp['modality'].value_counts().as_matrix()

        colorbrew_lib = ['Greens','Oranges','Purples']
        row_labels_pal = []
        for ii in range(len(modality_freq)):
            imf_pal_tmp = sns.color_palette(colorbrew_lib[ii], modality_freq[ii])
            row_labels_pal = row_labels_pal + imf_pal_tmp
        df_tmp.sort_values('modality',inplace =True)
        idx_lst = df_tmp.index.tolist()
        a,b = idx_lst.index(1), idx_lst.index(5)
        idx_lst[b], idx_lst[a] = idx_lst[a], idx_lst[b]
        row_labels_order = df_tmp.loc[idx_lst]
        df_tmp2 = pd.DataFrame([['Modality-specific metrics','other']], columns=['imf_label','modality'])
        row_labels_order = row_labels_order.append(df_tmp2,ignore_index=True)
        row_labels_pal = row_labels_pal + sns.color_palette(colorbrew_lib[ii+1], 1)

        row_labels_lut = dict(zip(row_labels_order['imf_label'].tolist(), row_labels_pal))
        row_labels_colors = pd.Series(row_labels, index=the_tab2.index).map(row_labels_lut)
        df_row_labels_colors = pd.DataFrame(dict(ImageFeature=row_labels_colors))
        df_row_labels_colors.columns = ['']
        row_labels_order = row_labels_order['imf_label'].tolist()
        # print(df_row_labels_colors)

    # hm_cmap = sns.light_palette((210, 90, 60), input="husl", as_cmap=True)
    hm_cmap = sns.cubehelix_palette(8, start=.7, rot=-.9, as_cmap=True)

    if mask is not None and vminmax:
        g = sns.clustermap(the_tab2, col_cluster=False, row_cluster=False, row_colors=df_row_labels_colors, cmap=hm_cmap,
                           mask=the_mask2, vmin=vminmax[0], vmax=vminmax[1], linewidths=0, xticklabels=True, yticklabels=True, figsize=(15, 15),annot=annot, fmt=fmt)
    else:
        g = sns.clustermap(the_tab2, col_cluster=False, row_cluster=False, row_colors=df_row_labels_colors, cmap=hm_cmap,
                           linewidths=0, xticklabels=True, yticklabels=True,figsize=(15, 15), annot=annot, fmt=fmt)


    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
    # plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45)

    for label in row_labels_order:
        g.ax_row_dendrogram.bar(0, 0, width=0, color=row_labels_lut[label], label=label, linewidth=0)
    if row_label_title:
        g.ax_row_dendrogram.legend(bbox_to_anchor=(1.4, 1.2), loc='best', title=row_label_title)

    g.cax.set_position([0.15, .2, .03, .45])
    g.cax.set_title(value_title)
    g.ax_heatmap.set(xlabel='', ylabel='')

    if fig_name:
        g.savefig(fig_name)
        print 'save heatmap plot to {}'.format(fig_name)
        plt.close()
    else:
        plt.show()

    return the_tab2


def heatmap_plot(the_tab, annot_labels=None, the_mask=None, vminmax=None, fig_name='', title=''):
    # sns.set()
    sns.set_context("poster")
    sns.set(font="futura std")
    fig, ax = plt.subplots(figsize=(15, 15))

    the_arg = {'data': the_tab, 'cmap':sns.light_palette((210, 90, 60), input="husl", as_cmap=True), 'ax': ax}

    if isinstance(annot_labels, np.ndarray):
        the_arg['annot'] = annot_labels
        the_arg['fmt'] = ''
    else:
        the_arg['annot'] = True
        the_arg['fmt'] = '1.3f'

    if the_mask is not None:
        the_arg['mask'] = the_mask

    if vminmax is not None:
        the_arg['vmin'] = vminmax[0]
        the_arg['vmax'] = vminmax[1]

    ax = sns.heatmap(**the_arg)

    # if isinstance(annot_labels, np.ndarray):
    #     ax = sns.heatmap(the_tab, annot=annot_labels, fmt='', cmap=sns.light_palette((210, 90, 60), input="husl", as_cmap=True))
    # else:
    #     ax = sns.heatmap(the_tab, annot=True, fmt='1.2f',cmap=sns.light_palette((210, 90, 60), input="husl", as_cmap=True))

    if title:
        ax.set_title(title)
    plt.setp(ax.get_yticklabels(), rotation=0)
    # plt.setp(ax.get_xticklabels(), rotation=90)

    if fig_name:
        plt.savefig(fig_name)
        print 'save heatmap plot to {}'.format(fig_name)
        plt.close()
    else:
        plt.show()


def rad(mic, a, b):
    return a * np.sqrt(np.log(mic * 1E4)) + b

def Annular_wedge_plot(df1, df3, filetag):
    inner_radius = 90
    outer_radius = 300 - 10
    minr = math.sqrt(math.log(.001 * 1E4))
    maxr = math.sqrt(math.log(100000 * 1E4))
    a = (outer_radius - inner_radius) / (minr - maxr)
    b = inner_radius - a * maxr

    width = 800
    height = 800
    PLOT_OPTS = dict(
        plot_width=width, plot_height=height, title="",
        x_axis_type=None, y_axis_type=None,
        x_range=(-420, 420), y_range=(-420, 420),
        min_border=0, outline_line_color="black",
        background_fill_color="#f0e1d2")

    # the angle for each bacteria type + legend (0.001 to 100)
    big_angle = 2.0 * np.pi / (len(df3) + 1)

    # why divide by 7? if you look at the figure, in each 'big_angle',
    # it's divided into 7 space since there are 3 kinds of recurred tumors
    # can change depending the # of category to plot within each circle
    N_tumor_type = len(df1['tumor_type'].unique())
    ndiv_small_angle = N_tumor_type * 2 + 1
    small_angle = big_angle / ndiv_small_angle

    p = figure(**PLOT_OPTS)

    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None


    N_tumor_type = len(df1.tumor_type.unique().tolist())
    print(df1.tumor_type.unique().tolist())
    print('N_tumor_type: {}'.format(N_tumor_type))
    tp_pal_dict = dict(zip([2, 3, 4, 5, 6, 7], [['#c64737', 'black'], Spectral3, Spectral4, Spectral5, Spectral6, Spectral7]))
    # print(tp_pal_dict)

    tumor_colors = tp_pal_dict[N_tumor_type]
    tumor_color_dict = OrderedDict(zip(df1.tumor_type.unique().tolist(), tumor_colors))
    # print(tumor_color_dict)

    radiomics_colors = ['#deebf7', '#e69584', '#bdbdbd']
    radiomics_color_dict = OrderedDict(zip(df3.Radiomics_type.unique().tolist(), radiomics_colors))
    # print(radiomics_color_dict)

    # annular wedges
    angles = np.pi / 2 - big_angle / 2 - df3.index.to_series() * big_angle
    colors = [radiomics_color_dict[rt] for rt in df3.Radiomics_type]
    p.annular_wedge(0, 0, inner_radius, outer_radius, -big_angle + angles, angles, color=colors)

    # small wedges
    for ii, tt in zip(range(ndiv_small_angle - 1, 1, -2), df1.tumor_type.tolist()):
        p.annular_wedge(0, 0, inner_radius, rad(df3[tt],a,b),
                        -big_angle + angles + (ii - 1) * small_angle, -big_angle + angles + ii * small_angle,
                        color=tumor_color_dict[tt])

    # circular axes and lables
    labels = np.power(10.0, np.arange(-3, 6))
    radii = a * np.sqrt(np.log(labels * 1E4)) + b
    p.circle(0, 0, radius=radii, fill_color=None, line_color="white")
    p.text(0, radii[:-1], [str(r) for r in labels[:-1]],
           text_font_size="8pt", text_align="center", text_baseline="middle")

    # radial axes
    p.annular_wedge(0, 0, inner_radius - 10, outer_radius + 10,
                    -big_angle + angles, -big_angle + angles, color="black")

    p.circle([-40, -40, -40], [-340, -360, -380], color=list(radiomics_color_dict.values()), radius=5)
    p.text([-30, -30, -30], [-340, -360, -380], text=radiomics_color_dict.keys(),
           text_font_size="8pt", text_align="left", text_baseline="middle")

    # tumor type legend
    rect_xx = [-40] * N_tumor_type

    if N_tumor_type > 5:
        rect_yy0 = 52
    elif N_tumor_type > 3:
        rect_yy0 = 26
    else:
        rect_yy0 = 18

    print('N_tumor_type={}, rect_yy0={}'.format(N_tumor_type, rect_yy0))
    rect_yy = range(rect_yy0, rect_yy0 - (N_tumor_type + 1) * 18, -18)

    txt_xx = [-15] * N_tumor_type
    txt_yy = rect_yy
    p.rect(rect_xx, rect_yy, width=30, height=13,
           color=list(tumor_color_dict.values()))
    p.text(txt_xx, txt_yy, text=list(tumor_color_dict),
           text_font_size="9pt", text_align="left", text_baseline="middle")
    show(p)

    # output to .html file
    htmlfile = '{}.html'.format(filetag)
    output_file(filename=htmlfile, title="Radiomics_PrimVsRecurTumors_Viz.py")

    # export to .png
    pngfile = '{}.png'.format(filetag)
    export_png(p, filename=pngfile)

    # # export to .svg (svg didn't seem to work or save the image correctly..)
    # p.output_backend = 'svg'
    # svgfile = '{}.svg'.format(filetag)
    # export_svgs(p, filename=svgfile)

def Plot_Annular(df_all, ptid_sides, savedir):
    for pp in ptid_sides:
        df1 = df_all[df_all['ptid_side'] == pp]
        val_vars = set(df1.columns.tolist()).symmetric_difference(['ptid_side', 'tumor_type'])

        # make an appropriate table
        df2 = pd.melt(df1, id_vars=['ptid_side', 'tumor_type'], value_vars=val_vars, var_name='Radiomics')
        df3 = df2.pivot(index='Radiomics', columns='tumor_type', values='value')
        df3.reset_index(inplace=True)

        # make another column to categorize radiomic feature to FOstats, shape and size and texture
        df3['Radiomics_type'] = df3['Radiomics'].apply(lambda x: re.split('_+', x)[0] if re.split('_+', x) else np.nan)

        filetag = '{}/Radiomics_primaryVSrecur_{}'.format(savedir,pp)
        Annular_wedge_plot(df1, df3, filetag=filetag)

def Plot_Radiomics_Corr(df_all, ptid_sides, fig_name=''):
    # compute pearson correlation coefficients for all primary to recur tumor pairs
    df_corr_all = pd.DataFrame()
    for ps in ptid_sides:
        df1 = df_all[df_all['ptid_side'] == ps]
        val_vars = set(df1.columns.tolist()).symmetric_difference(['ptid_side', 'tumor_type'])
        df2 = pd.melt(df1, id_vars=['ptid_side', 'tumor_type'], value_vars=val_vars, var_name='Radiomics')
        df3 = df2.pivot(index='Radiomics', columns='tumor_type', values='value')
        df3.reset_index(inplace=True)
        corr_mat = df3.corr()
        tmp1 = corr_mat.loc['Primary', :]
        the_ind = [ss for ss in tmp1.index.tolist() if ss != 'Primary']
        tmp2 = tmp1[the_ind].reset_index()
        tmp2['prim_recur pair name'] = tmp2['tumor_type'].map(lambda x: 'Primary_{}'.format(x))
        tmp2.drop(columns=['tumor_type'], inplace=True)
        tmp2['Tumor Pair Index'] = range(1, len(tmp2) + 1)
        tmp2['pt_id'] = ps
        tmp2.rename(columns={'Primary': 'Pearson Correlation'}, inplace=True)
        df_corr_all = pd.concat([df_corr_all, tmp2], ignore_index=True)

    sns.set(style="ticks")
    g = sns.FacetGrid(df_corr_all, col="pt_id", hue="pt_id", col_wrap=4, size=3)
    g.map(plt.plot, 'Tumor Pair Index', 'Pearson Correlation', marker="o", ms=6, color='green')

    # Adjust the tick positions and labels
    g.set(xlim=(0, 6), ylim=(0, 1.2))

    # Adjust the arrangement of the plots
    g.fig.tight_layout(w_pad=1)

    if fig_name:
        g.savefig(fig_name)
        print 'save facetgrid line plot to {}'.format(fig_name)
        plt.close()
    else:
        plt.show()

def Plot_Learning_OutCome_Heatmap(the_df, xlabel_name='', xticklabel_names='', fig_name=''):
    # plot the auc in a heatmap
    the_tab = the_df.pivot('clf_name', 'Dep_var', 'AUC_hat')
    the_label = the_df.pivot('clf_name', 'Dep_var', 'label')
    sns.set(font='Arial')
    sns.set_context('paper', font_scale=1.2)

    plt.figure(figsize=(11, 8.5))
    g = sns.heatmap(the_tab, annot=the_label, fmt='', annot_kws={'size': 10}, linewidth=0.2, cmap='Blues',
                    cbar_kws=dict(use_gridspec=False, location="top", label='AUC', shrink=0.3, anchor=(0.0, -0.2)))
    # g.set(xlabel='Clinical outcome', ylabel='Classification algorithm')
    g.set_xlabel(xlabel_name, {'weight': 'bold', 'size': 13})
    g.set_ylabel('Classification algorithm', {'weight': 'bold', 'size': 13})
    g.set_xticklabels(xticklabel_names, {'weight': 'bold', 'size': 11})
    g.set_yticklabels(g.get_yticklabels(), {'weight': 'bold', 'size': 11})

    if fig_name:
        plt.savefig(fig_name, bbox_inches='tight')
        print('save the heatmap figure: {}'.format(fig_name))
        plt.close()
    else:
        plt.show()


