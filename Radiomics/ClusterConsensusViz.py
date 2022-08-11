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
import itertools as itt


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
                       figsize=(15, 15),vmin=-5.,vmax=5.)
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


def ClustermapPlot(the_df,varname,fig_fname,fig_title=None, var_title=None, var_color_pal=None):
    """
    
    :param the_df: the input dataframe
    :param varname: a list of variable name, list
    :param fig_title: figure title
    :param fig_fname: figure file name to save
    :param var_title: variable title that one would like to show on the figure (instead of varname)
    :param var_color_pal: a list of variable color palette for cluster plots, list (if more than one variables are defined)
    :return: n/a
    """

    # define all image features of interest
    texture_cols = ['texture_autocorrelation', 'texture_cluster_prominence', 'texture_cluster_shade',
                    'texture_cluster_tendency', 'texture_contrast', 'texture_correlation', 'texture_diff_entropy',
                    'texture_dissimilarity',
                    'texture_energy', 'texture_entropy', 'texture_homogeneity1', 'texture_homogeneity2', 'texture_idmn',
                    'texture_idn', 'texture_inv_var', 'texture_maxprob', 'texture_sum_avg', 'texture_sum_entropy',
                    'texture_sum_var']

    img_feature_names = ['FOstats_energy', 'FOstats_entropy', 'FOstats_kurtosis', 'FOstats_mean', 'FOstats_min',
                         'FOstats_max', 'FOstats_skewness', 'FOstats_uniformity', 'FOstats_variance',
                         'ShapeSize_compactness1', 'ShapeSize_compactness2', 'ShapeSize_max_euc_dis',
                         'ShapeSize_spherical_disproportion', 'ShapeSize_sphericity', 'ShapeSize_surf_area_cm2',
                         'ShapeSize_surface2volratio', 'ShapeSize_vol_cm3'] + [ss + '_avg' for ss in texture_cols]

    # make sure there's no NANs
    df_sub = the_df.ix[:,['cs_class','ptid_side'] + varname + img_feature_names]
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

    # this is for single variable case
    # var_type = df_sub[varname].tolist()
    # cs_class = df_sub['cs_class'].tolist()
    #
    # # build in tumor status for each primary_ID
    # mcol = pd.MultiIndex.from_tuples(zip(var_type, cs_class), names= varname + ['CS_class'])

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
    # var_pal = sns.color_palette('hls',len(df_pivot_col_unique))
    # var_pal = sns.color_palette('Blues',len(df_pivot_col_unique))

    # this is for single variable case
    # # determine the color palette for a given variable
    # if var_color_pal:
    #     var_pal = var_color_pal
    # else:
    #     var_pal = sns.color_palette('Blues', len(df_pivot_col_unique))
    #
    # var_labels = [df_pivot_col_unique[i][0] for i in range(len(df_pivot_col_unique))]
    # unique_var_labels = list(set(var_labels))
    # var_lut = dict(zip(unique_var_labels, var_pal))
    # if var_title:
    #     clinical_data_colors[var_title] = pd.Series(var_labels, index=df_pivot_col_unique).map(var_lut)
    # else:
    #     clinical_data_colors[varname] = pd.Series(var_labels, index=df_pivot_col_unique).map(var_lut)
    #
    # csclass_labels = [df_pivot_col_unique[i][1] for i in range(len(df_pivot_col_unique))]

    # use this dict to save the unique_var_labels for displaying legend
    var_labels_dict = {}
    var_lut_dict = {}
    for ii,vv in enumerate(varname):
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

    csclass_labels = [df_pivot_col_unique[i][ii+1] for i in range(len(df_pivot_col_unique))]
    unique_csclass_labels = list(set(csclass_labels))
    # ch_pal = sns.color_palette('cubehelix',10)
    # Ndiv = int(round(len(ch_pal)/float(len(unique_csclass_labels))))
    # csclass_pal = ch_pal[::Ndiv]
    csclass_pal = sns.light_palette('gray',len(unique_csclass_labels))
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
                       figsize=(15, 15),vmin=-5.,vmax=5.)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    # display the legend for the image feature colormap by adding a empty bar plot but display the legend
    for label in list(set(img_feature_labels)):
        g.ax_row_dendrogram.bar(0, 0, color=img_feature_lut[label], label=label, linewidth=0)
    g.ax_row_dendrogram.legend(bbox_to_anchor=(1.4, 1.2), loc='best', title='Image Features')

    for label in unique_csclass_labels:
        g.ax_col_dendrogram.bar(0, 0, color=csclass_lut[label], label=label, linewidth=0)
    g.ax_col_dendrogram.legend(bbox_to_anchor=(0.85,0.56), loc='center', title='Cluster Consensus Class')

    # # for single variable use
    # for label in unique_var_labels:
    #     g.ax_heatmap.bar(0, 0, color=var_lut[label], label=label, linewidth=0)
    # g.ax_heatmap.legend(bbox_to_anchor=(0.65, 1.23), loc='best', title=varname,ncol=int(math.ceil(len(unique_var_labels)*0.5)))

    box_xx = 0.65
    box_yy = 1.23
    for ii,kk in enumerate(var_labels_dict.keys()):
        unique_var_labels = var_labels_dict[kk]
        var_lut = var_lut_dict[kk]
        hh_all = []
        for label in unique_var_labels:
            hh = g.ax_heatmap.bar(0,0, color=var_lut[label], label=label, linewidth=0)
            hh_all.append(hh)
        lg = g.ax_heatmap.legend(handles=hh_all,bbox_to_anchor=(box_xx-0.3*ii, box_yy), loc='best', title=kk,ncol=int(math.ceil(len(unique_var_labels) * 0.5)))
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

def str_find_dict(the_str,the_dict):
    for kk in the_dict.keys():
        ss = the_str.lower().replace(' ','')  #lower-case and remove white space
        match_list = the_dict[kk][0]
        for mm in match_list:
            # print 'the_str:{}, match: {}'.format(ss,mm)
            if ss.find(mm.lower()) >= 0:
                return kk
    print '{}: cannot find a key in the dictionary that matches the string'.format(the_str)
    return 0

def recur_func2(x):
    """
    out of the patients who was disease free in his/her life, 0 if no recur, 1 if recur
    :param x: entry in the recurrent type summary data
    :return: return nan if N/A, 0 if disease free/no recur, 1 if some sort of recurrence or never disease free
    """
    if x is np.nan or x.lower() == 'unknown' or x.lower() == 'never disease free':
        return np.nan
    elif x.lower().find('recur') >= 0:
        return 1
    elif x == 'NONE/DISEASE FREE':
        return 0

def recur_func1(x):
    if x is np.nan or x.lower() == 'unknown':
        return np.nan
    if x.lower() == 'never disease free':
        return 2
    elif x == 'NONE/DISEASE FREE':
        return 0
    else:
        return 1

def recur_func0(x):
    if x == 'NEVER DISEASE FREE':
        return 2
    # elif x == 'NONE/DISEASE FREE' or x is np.nan or x.lower() == 'unknown':
    elif x == 'NONE/DISEASE FREE' or x is np.nan:
        return 0
    else:
        return 1

def diseasefree5yr_func(row):
    """
    :param x: a row in a pandas dataframe
    :return: return nan if N/A, 0 if disease free/no recur after 5 years, 1 if some sort of recurrence or never disease free
    """
    dT = pd.to_datetime(row['Last contact/death']) - pd.to_datetime(row['Date of Primary Dx'])
    if isinstance(dT,pd.tslib.NaTType):
        dT_yr = -1
    else:
        dT_yr = dT.days/365.

    if row['Recurrence Type Summary'] is np.nan or row['Recurrence Type Summary'].lower() == 'unknown':
        return np.nan
    elif row['Recurrence Type Summary'] == 'NONE/DISEASE FREE' and dT_yr >= 5.0:
        return 1
    else:
        return 0

# # Recurrence type: 0: assume NaN or NONE/DISEASE FREE are no recurrence
# # 1: bone recur, 2: local recur, 3: distant/systemic recur
# # -1: unknown recur type AND never disease free
# def recur_type_func(x):
#     if x is np.nan or x == 'NONE/DISEASE FREE':
#         return 0
#     elif x == 'DIST RECUR, BONE':
#         return 1
#     elif re.search(r'LOCAL RECUR[\w\s.]*',x):
#         return 2
#     elif re.search(r'DIST RECUR[\w\s.]*',x) and x != 'DIST RECUR, BONE':
#         return 3
#     else:
#         return 4

# 0: no recur, 1: local recur, 2: distant recur (bone and others),
# nan or 'unknown' or 'never disease free' or 'recurred, type unknow': N/A,
def recur_type_func(x):
    if x is np.nan or x.lower().find('unknown') >= 0 or x.lower() == 'never disease free':
        return np.nan
    elif x == 'NONE/DISEASE FREE':
        return 0
    elif re.search(r'LOCAL RECUR[\w\s.]*',x):
        return 1
    elif re.search(r'DIST RECUR[\w\s.]*',x):
        return 2

def Tstage_func(x):
    tstage_dict = {'T0':(['T0','TX','calc','Tis'],0), 'T1':(['T1'],1), 'T2':(['T2'],2),'T3':(['T3'],3),'T4':(['T4'],4)}
    # add the numerical values to each or combined category
    # tstage_k = tstage_dict.keys()
    # for i in range(len(tstage_k)):
    #     tstage_dict[tstage_k[i]] = (tstage_dict[tstage_k[i]],i)

    if x is np.nan:
        # print 'x is nan'
        return x
    else:
        kk = str_find_dict(x,tstage_dict)
        if kk != 0:
            return tstage_dict[kk][1]
        else:
            return np.nan

def Nstage_func(x):
    nstage_dict = {'N0': (['N0','NX','calc'],0), 'N1': (['N1'],1), 'N2': (['N2'],2), 'N3': (['N3'],3)}

    # # add the numerical values to each or combined category
    # nstage_k = nstage_dict.keys()
    # for i in range(len(nstage_k)):
    #     nstage_dict[nstage_k[i]] = (nstage_dict[nstage_k[i]], i)

    if x is np.nan:
        return x
    else:
        kk = str_find_dict(x, nstage_dict)
        if kk != 0:
            return nstage_dict[kk][1]
        else:
            return np.nan

def Overallstage_func(x):
    ostage_dict = {'0': (['0','calc'],0), '1': (['1'],1), '2': (['2'],2), '3': (['3'],3), '4':(['4'],4)}

    # # add the numerical values to each or combined category
    # ostage_k = ostage_dict.keys()
    # for i in range(len(ostage_k)):
    #     ostage_dict[ostage_k[i]] = (ostage_dict[ostage_k[i]], i)

    if x is np.nan:
        return x
    else:
        kk = str_find_dict(x, ostage_dict)
        if kk != 0:
            return ostage_dict[kk][1]
        else:
            return np.nan

def BC_subtype_func(row):
# split the tumors into tumor subtypes (suggested by E. Jones)
# 0: HR+/HER2-, 1: HR+/HER2+, 2: HR-/HER2+, 3: triple neg, HR+ is ER+ or PR+
    is_HR_pos = (row['Sjoerd_PR'] == 1) or (row['Sjoerd_ER'] == 1)
    if is_HR_pos and row['Sjoerd_HER2'] == 0:
        return 0
    elif is_HR_pos and row['Sjoerd_HER2'] == 1:
        return 1
    elif ~is_HR_pos and row['Sjoerd_HER2'] == 1:
        return 2
    elif ~is_HR_pos and row['Sjoerd_HER2'] == 0:
        return 3
    else:
        return np.nan

if __name__ == '__main__':
    # gather outcome data
    rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'
    outcome_fname = '{}/her2_ClinicalData/her2_outcome.csv'.format(rootdir)
    her2_outcome_df = pd.read_csv(outcome_fname,dtype={'MRN': str,'Recurrence Type Summary': str})

    # Recurrence_CR (recurrence info based on Cancer Registry field 'Recurrence Type Summary'): -1, always have cancer, 0, no recurrence, and 1, recur
    # her2_outcome_df['Recurrence_CR_0'] = her2_outcome_df['Recurrence Type Summary'].map(recur_func0)
    her2_outcome_df['Recurrence_CR_1'] = her2_outcome_df['Recurrence Type Summary'].map(recur_func1)
    her2_outcome_df['Recurrence_CR_2'] = her2_outcome_df['Recurrence Type Summary'].map(recur_func2)
    her2_outcome_df['Diseasefree_5yr'] = her2_outcome_df.apply(diseasefree5yr_func, axis=1)
    her2_outcome_df['Recurrence_Type'] = her2_outcome_df['Recurrence Type Summary'].map(recur_type_func)

    # catergorize tumor and lymph node stages
    her2_outcome_df['T_stage'] = her2_outcome_df['T-stage at surgery'].map(Tstage_func)
    her2_outcome_df['N_stage'] = her2_outcome_df['N-stage at surgery'].map(Nstage_func)
    her2_outcome_df['Overall_stage'] = her2_outcome_df['Overall stage'].map(Overallstage_func)

    # # Recurrence field does not seem reliable...
    # her2_outcome_df['Recurrence'] = her2_outcome_df['Recurrence'].apply(lambda x: int(x) if x != '#VALUE!' and x is not np.nan else np.nan)
    # # disease status = 0 when time interval since diagnosis
    # her2_outcome_df['disease status'] = her2_outcome_df['Time interval since diagnosis'].apply(lambda x: 0 if x < 5.0 else 1)

    the_outcome_name = ['MRN','Diseasefree_5yr','Recurrence_CR_1', 'Recurrence_CR_2', 'Recurrence_Type','T_stage','N_stage','Overall_stage','T-stage at surgery']
    the_outcome_df = her2_outcome_df.loc[:, the_outcome_name]

    # # PET image feature data
    # # dist_method_list = ['euclidean','spearman','pearson']
    # # cluster_method_dict = {'pam': None, 'km': None, 'kmdist': None,'hc':['centroid','median','average','mcquitty','complete']}
    # # Ncluster_list = [2, 3, 4, 5]
    # # # theimgbin = [32, 64, 128]
    # # # theglcmbin = [32, 64, 128]
    # # theimgbin = [128]
    # # theglcmbin = [64]
    #
    # # for testing purpose
    # dist_method_list = ['spearman']
    # cluster_method_dict = {'hc': ['mcquitty']}
    # Ncluster_list = [3]
    # theimgbin = [128]
    # theglcmbin = [64]
    #
    # for tt, bb in itt.product(theimgbin, theglcmbin):
    #     pet_dir = '{}/her2_Analysis/PET/IsoVoxel_IMGBIN{}_GLCMBIN{}'.format(rootdir,tt,bb)
    #     df_fname = '{}/PETdataAll_glcmNbin{}_normNbin{}.csv'.format(pet_dir,bb,tt)
    #     data_df = pd.read_csv(df_fname, dtype={'pt_mrn': str, 'PRIMARY_ID': str, 'MRN': str, 'Anon_Accession_Num': str})
    #     data_df['BC_subtype'] = data_df.apply(BC_subtype_func, axis=1)
    #     jdf = pd.merge(data_df, the_outcome_df, on='MRN')
    #     print 'read the file {}'.format(df_fname)
    #     # print data_df.shape
    #     # print jdf.shape
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
    #
    #                     # # fig_title = 'cluster method: {}, linkage: {}, distance method: {}'.format(k, ll, dd)
    #
    #                     fig_fname = '{}/clustermap_BCsubtype_kmax{}.pdf'.format(cc_dir,ncluster)
    #                     ClustermapPlot(the_df,['BC_subtype'],fig_fname,var_title=['Breast cancer subtype'])
    #
    #                     # fig_fname = '{}/clustermap_Tstage_kmax{}.pdf'.format(cc_dir,ncluster)
    #                     # ClustermapPlot(the_df,['T_stage'],fig_fname,var_title=['T-stage'])
    #                     #
    #                     # fig_fname = '{}/clustermap_N_Overall_stage_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # varname = ['N_stage','Overall_stage']
    #                     # vartitle = ['N-stage','Overall-stage']
    #                     # N_unique_class = [len(the_df[varname[ii]].unique()) for ii in range(len(varname))]
    #                     # color_name = ['royalblue','forestgreen']
    #                     # var_pal_list = [sns.light_palette(cc, nn) for nn, cc in zip(N_unique_class, color_name)]
    #                     # ClustermapPlot(the_df,varname, fig_fname,var_title=vartitle,var_color_pal=var_pal_list)
    #                     #
    #                     # # fig_fname = '{}/clustermap_Nstage_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # # ClustermapPlot(the_df, 'N_stage', fig_title, fig_fname,var_title='N-stage')
    #                     # #
    #                     # # fig_fname = '{}/clustermap_Overallstage_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # # ClustermapPlot(the_df, 'Overall_stage', fig_title, fig_fname,var_title='Overall stage')
    #                     #
    #                     # fig_fname = '{}/clustermap_TripleNeg_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # tripneg_pal = [(1., 1., 1.),sns.color_palette('Set2', 10)[2]]
    #                     # ClustermapPlot(the_df, ['TripleNeg'], fig_fname, var_title=['Triple Neg'], var_color_pal=[tripneg_pal])
    #                     #
    #                     # fig_fname = '{}/clustermap_Grade_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # light_pal = sns.light_palette((210, 90, 60), input="husl")
    #                     # grade_pal = light_pal[1:len(light_pal):2]
    #                     # ClustermapPlot(the_df, ['Sjoerd_Grade'], fig_fname, var_title=['Grade'], var_color_pal=[grade_pal])
    #                     #
    #                     # fig_fname = '{}/clustermap_Histology_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # colors = ['windows blue', 'amber', 'greyish', 'brick']
    #                     # histology_pal = sns.xkcd_palette(colors)
    #                     # ClustermapPlot(the_df, ['Marjan_Histology'], fig_fname, var_title=['Histology'],var_color_pal=[histology_pal])
    #
    #                     # fig_fname = '{}/clustermap_recurstatus0_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # recur_pal0= [(1., 1., 1.), sns.color_palette('Set2', 10)[3], sns.color_palette('Set2', 10)[7]]
    #                     # ClustermapPlot(the_df, ['Recurrence_CR_0'], fig_fname, var_title=['Recurrence'],var_color_pal=[recur_pal0])
    #                     #
    #                     # fig_fname = '{}/clustermap_recurstatus1_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # recur_pal1 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3],sns.color_palette('Set2', 10)[7]]
    #                     # ClustermapPlot(the_df, ['Recurrence_CR_1'], fig_fname, var_title=['Recurrence'],var_color_pal=[recur_pal1])
    #                     #
    #                     # fig_fname = '{}/clustermap_recurstatus2_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # recur_pal2 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3]]
    #                     # ClustermapPlot(the_df, ['Recurrence_CR_2'], fig_fname, var_title=['Recurrence'],var_color_pal=[recur_pal2])
    #
    #                     # fig_fname = '{}/clustermap_df_5yr_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # df5yr_pal = [(1., 1., 1.), sns.color_palette('Set2', 10)[3]]
    #                     # ClustermapPlot(the_df, ['Diseasefree_5yr'], fig_fname, var_title=['5-yr Disease Free'],
    #                     #                var_color_pal=[df5yr_pal])
    #
    #                     # fig_fname = '{}/clustermap_recurtype_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # recurtype_pal = sns.light_palette('navy', n_colors=len(the_df['Recurrence_Type'].unique()), reverse=True)
    #                     # ClustermapPlot(the_df, ['Recurrence_Type'], fig_fname, var_title=['Recurrence Type'],var_color_pal=[recurtype_pal])
    #
    #                     # fig_fname = '{}/clustermap_df5yr_TN_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # varname = ['Diseasefree_5yr','TripleNeg']
    #                     # vartitle = ['5-yr Disease free','Triple Negative']
    #                     # pal1  = [(1., 1., 1.), sns.color_palette('Set2', 10)[3]]
    #                     # pal2 = [(1., 1., 1.), sns.color_palette('Set2', 10)[8]]
    #                     # var_pal_list = [pal1,pal2]
    #                     # ClustermapPlot(the_df,varname, fig_fname,var_title=vartitle,var_color_pal=var_pal_list)
    #
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
    #                     # print the_df.shape
    #
    #                     # # fig_title = 'cluster method: {}, distance method: {}'.format(k, dd)
    #                     # fig_fname = '{}/clustermap_Tstage_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # ClustermapPlot(the_df, ['T_stage'], fig_fname, var_title=['T-stage'])
    #                     #
    #                     # fig_fname = '{}/clustermap_N_Overall_stage_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # varname = ['N_stage', 'Overall_stage']
    #                     # vartitle = ['N-stage', 'Overall-stage']
    #                     # N_unique_class = [len(the_df[varname[ii]].unique()) for ii in range(len(varname))]
    #                     # color_name = ['royalblue', 'forestgreen']
    #                     # var_pal_list = [sns.light_palette(cc, nn) for nn, cc in zip(N_unique_class, color_name)]
    #                     # ClustermapPlot(the_df, varname, fig_fname, var_title=vartitle, var_color_pal=var_pal_list)
    #
    #                     # fig_fname = '{}/clustermap_Nstage_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # ClustermapPlot(the_df, 'N_stage', fig_title, fig_fname,var_title='N-stage')
    #                     #
    #                     # fig_fname = '{}/clustermap_Overallstage_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # ClustermapPlot(the_df, 'Overall_stage', fig_title, fig_fname,var_title='Overall stage')
    #
    #                     # fig_fname = '{}/clustermap_TripleNeg_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # tripneg_pal = [(1., 1., 1.), sns.color_palette('Set2', 10)[2]]
    #                     # ClustermapPlot(the_df, ['TripleNeg'], fig_fname, var_title=['Triple Neg'],
    #                     #                var_color_pal=[tripneg_pal])
    #                     #
    #                     # fig_fname = '{}/clustermap_Grade_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # light_pal = sns.light_palette((210, 90, 60), input="husl")
    #                     # grade_pal = light_pal[1:len(light_pal):2]
    #                     # ClustermapPlot(the_df, ['Sjoerd_Grade'], fig_fname, var_title=['Grade'],
    #                     #                var_color_pal=[grade_pal])
    #                     #
    #                     # fig_fname = '{}/clustermap_Histology_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # colors = ['windows blue', 'amber', 'greyish', 'brick']
    #                     # histology_pal = sns.xkcd_palette(colors)
    #                     # ClustermapPlot(the_df, ['Marjan_Histology'], fig_fname, var_title=['Histology'],
    #                     #                var_color_pal=[histology_pal])
    #
    #                     # fig_fname = '{}/clustermap_recurstatus0_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # recur_pal0 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3], sns.color_palette('Set2', 10)[7]]
    #                     # ClustermapPlot(the_df, ['Recurrence_CR_0'], fig_fname, var_title=['Recurrence'],
    #                     #                var_color_pal=[recur_pal0])
    #                     #
    #                     # fig_fname = '{}/clustermap_recurstatus1_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # recur_pal1 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3], sns.color_palette('Set2', 10)[7]]
    #                     # ClustermapPlot(the_df, ['Recurrence_CR_1'], fig_fname, var_title=['Recurrence'],
    #                     #                var_color_pal=[recur_pal1])
    #                     #
    #                     # fig_fname = '{}/clustermap_recurstatus2_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # recur_pal2 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3]]
    #                     # ClustermapPlot(the_df, ['Recurrence_CR_2'], fig_fname, var_title=['Recurrence'],
    #                     #                var_color_pal=[recur_pal2])
    #
    #                     # fig_fname = '{}/clustermap_df_5yr_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # df5yr_pal = [(1., 1., 1.), sns.color_palette('Set2', 10)[3]]
    #                     # ClustermapPlot(the_df, ['Diseasefree_5yr'], fig_fname, var_title=['5-yr Disease Free'],
    #                     #                var_color_pal=[df5yr_pal])
    #
    #                     # fig_fname = '{}/clustermap_recurtype_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # recurtype_pal = sns.light_palette('navy', n_colors=len(the_df['Recurrence_Type'].unique()),
    #                     #                                   reverse=True)
    #                     # ClustermapPlot(the_df, ['Recurrence_Type'], fig_fname, var_title=['Recurrence Type'],
    #                     #                var_color_pal=[recurtype_pal])
    #
    #                     # fig_fname = '{}/clustermap_df5yr_TN_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     # varname = ['Diseasefree_5yr', 'TripleNeg']
    #                     # vartitle = ['5-yr Disease free', 'Triple Negative']
    #                     # pal1 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3]]
    #                     # pal2 = [(1., 1., 1.), sns.color_palette('Set2', 10)[8]]
    #                     # var_pal_list = [pal1, pal2]
    #                     # ClustermapPlot(the_df, varname, fig_fname, var_title=vartitle, var_color_pal=var_pal_list)

    # MRI image feature data
    # theTP = [1, 2, 3]
    # theglcmBin = [64, 128, 256]
    # theTP = [2]
    # theglcmBin = [128]
    # Ncluster_list = [2,3,4,5]
    # dist_method_list = ['euclidean','spearman','pearson']
    # cluster_method_dict = {'pam': None, 'km': None, 'kmdist': None,'hc':['centroid','median','average','mcquitty','complete']}

    # for testing purpose
    theTP = [2]
    theglcmBin = [128]
    Ncluster_list = [2,3,4]
    dist_method_list = ['euclidean']
    cluster_method_dict = {'pam': None}

    for tt, bb in itt.product(theTP, theglcmBin):
        mri_dir = '{}/her2_Analysis/MRI/IsoVoxel_TP{}_GLCMBIN{}'.format(rootdir, tt, bb)
        df_fname = '{}/MRIdataAll_tp{}_Nbin{}.csv'.format(mri_dir, tt, bb)
        data_df = pd.read_csv(df_fname, dtype={'PRIMARY_ID': str, 'MRN': str})
        data_df['BC_subtype'] = data_df.apply(BC_subtype_func, axis=1)
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

                        # fig_title = 'cluster method: {}, linkage: {}, distance method: {}'.format(k, ll, dd)

                        fig_fname = '{}/clustermap_BCsubtype_kmax{}.pdf'.format(cc_dir, ncluster)
                        ClustermapPlot(the_df, ['BC_subtype'], fig_fname, var_title=['Breast cancer subtype'])

                        # fig_fname = '{}/clustermap_Tstage_kmax{}.pdf'.format(cc_dir, ncluster)
                        # ClustermapPlot(the_df, ['T_stage'], fig_fname, var_title=['T-stage'])
                        #
                        # fig_fname = '{}/clustermap_Nstage_kmax{}.pdf'.format(cc_dir, ncluster)
                        # ClustermapPlot(the_df, ['N_stage'], fig_fname,var_title=['N-stage'])
                        #
                        # fig_fname = '{}/clustermap_Overallstage_kmax{}.pdf'.format(cc_dir, ncluster)
                        # ClustermapPlot(the_df, ['Overall_stage'], fig_fname,var_title=['Overall stage'])
                        #
                        # fig_fname = '{}/clustermap_TripleNeg_kmax{}.pdf'.format(cc_dir, ncluster)
                        # tripneg_pal = [(1., 1., 1.), sns.color_palette('Set2', 10)[2]]
                        # ClustermapPlot(the_df, ['TripleNeg'], fig_fname, var_title=['Triple Neg'],
                        #                var_color_pal=[tripneg_pal])
                        #
                        # fig_fname = '{}/clustermap_Grade_kmax{}.pdf'.format(cc_dir, ncluster)
                        # light_pal = sns.light_palette((210, 90, 60), input="husl")
                        # grade_pal = light_pal[1:len(light_pal):2]
                        # ClustermapPlot(the_df, ['Sjoerd_Grade'], fig_fname, var_title=['Grade'],
                        #                var_color_pal=[grade_pal])
                        #
                        # fig_fname = '{}/clustermap_Histology_kmax{}.pdf'.format(cc_dir, ncluster)
                        # colors = ['windows blue', 'amber', 'greyish', 'brick']
                        # histology_pal = sns.xkcd_palette(colors)
                        # ClustermapPlot(the_df, ['Marjan_Histology'], fig_fname, var_title=['Histology'],
                        #                var_color_pal=[histology_pal])
                        #
                        # # fig_fname = '{}/clustermap_recurstatus0_kmax{}.pdf'.format(cc_dir, ncluster)
                        # # recur_pal0 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3], sns.color_palette('Set2', 10)[7]]
                        # # ClustermapPlot(the_df, ['Recurrence_CR_0'], fig_fname, var_title=['Recurrence'],var_color_pal=[recur_pal0])
                        #
                        # fig_fname = '{}/clustermap_recurstatus1_kmax{}.pdf'.format(cc_dir, ncluster)
                        # recur_pal1 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3], sns.color_palette('Set2', 10)[7]]
                        # ClustermapPlot(the_df, ['Recurrence_CR_1'], fig_fname, var_title=['Recurrence'],var_color_pal=[recur_pal1])
                        #
                        # fig_fname = '{}/clustermap_recurstatus2_kmax{}.pdf'.format(cc_dir, ncluster)
                        # recur_pal2 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3]]
                        # ClustermapPlot(the_df, ['Recurrence_CR_2'], fig_fname, var_title=['Recurrence'],var_color_pal=[recur_pal2])
                        #
                        # fig_fname = '{}/clustermap_recurtype_kmax{}.pdf'.format(cc_dir, ncluster)
                        # recurtype_pal = sns.light_palette('navy', n_colors=len(the_df['Recurrence_Type'].unique()),reverse=True)
                        # ClustermapPlot(the_df, ['Recurrence_Type'], fig_fname, var_title=['Recurrence Type'],var_color_pal=[recurtype_pal])
                else:
                    for dd in dist_method_list:
                        cc_dir = '{}/ConsensusCluster_{}_{}'.format(mri_dir,k,dd)
                        cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir, ncluster)
                        print 'read file {}'.format(cs_class_fname)
                        cs_class_df = pd.read_csv(cs_class_fname)
                        cs_class_df.columns = ['ptid_side', 'cs_class']

                        # combine the cs_class to the_df
                        the_df = pd.merge(jdf, cs_class_df, how='left', on='ptid_side')

                        # fig_title = 'cluster method: {}, linkage: {}, distance method: {}'.format(k, ll, dd)

                        fig_fname = '{}/clustermap_BCsubtype_kmax{}.pdf'.format(cc_dir, ncluster)
                        ClustermapPlot(the_df, ['BC_subtype'], fig_fname, var_title=['Breast cancer subtype'])

                        # fig_fname = '{}/clustermap_Tstage_kmax{}.pdf'.format(cc_dir, ncluster)
                        # ClustermapPlot(the_df, ['T_stage'], fig_fname, var_title=['T-stage'])
                        #
                        # fig_fname = '{}/clustermap_Nstage_kmax{}.pdf'.format(cc_dir, ncluster)
                        # ClustermapPlot(the_df, ['N_stage'], fig_fname,var_title=['N-stage'])
                        #
                        # fig_fname = '{}/clustermap_Overallstage_kmax{}.pdf'.format(cc_dir, ncluster)
                        # ClustermapPlot(the_df, ['Overall_stage'], fig_fname,var_title=['Overall stage'])
                        #
                        # fig_fname = '{}/clustermap_TripleNeg_kmax{}.pdf'.format(cc_dir, ncluster)
                        # tripneg_pal = [(1., 1., 1.), sns.color_palette('Set2', 10)[2]]
                        # ClustermapPlot(the_df, ['TripleNeg'], fig_fname, var_title=['Triple Neg'],
                        #                var_color_pal=[tripneg_pal])
                        #
                        # fig_fname = '{}/clustermap_Grade_kmax{}.pdf'.format(cc_dir, ncluster)
                        # light_pal = sns.light_palette((210, 90, 60), input="husl")
                        # grade_pal = light_pal[1:len(light_pal):2]
                        # ClustermapPlot(the_df, ['Sjoerd_Grade'], fig_fname, var_title=['Grade'],
                        #                var_color_pal=[grade_pal])
                        #
                        # fig_fname = '{}/clustermap_Histology_kmax{}.pdf'.format(cc_dir, ncluster)
                        # colors = ['windows blue', 'amber', 'greyish', 'brick']
                        # histology_pal = sns.xkcd_palette(colors)
                        # ClustermapPlot(the_df, ['Marjan_Histology'], fig_fname, var_title=['Histology'],
                        #                var_color_pal=[histology_pal])
                        #
                        # # fig_fname = '{}/clustermap_recurstatus0_kmax{}.pdf'.format(cc_dir, ncluster)
                        # # recur_pal0 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3], sns.color_palette('Set2', 10)[7]]
                        # # ClustermapPlot(the_df, ['Recurrence_CR_0'], fig_fname, var_title=['Recurrence'],
                        # #                var_color_pal=[recur_pal0])
                        #
                        # fig_fname = '{}/clustermap_recurstatus1_kmax{}.pdf'.format(cc_dir, ncluster)
                        # recur_pal1 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3], sns.color_palette('Set2', 10)[7]]
                        # ClustermapPlot(the_df, ['Recurrence_CR_1'], fig_fname, var_title=['Recurrence'],
                        #                var_color_pal=[recur_pal1])
                        #
                        # fig_fname = '{}/clustermap_recurstatus2_kmax{}.pdf'.format(cc_dir, ncluster)
                        # recur_pal2 = [(1., 1., 1.), sns.color_palette('Set2', 10)[3]]
                        # ClustermapPlot(the_df, ['Recurrence_CR_2'], fig_fname, var_title=['Recurrence'],
                        #                var_color_pal=[recur_pal2])
                        #
                        # fig_fname = '{}/clustermap_recurtype_kmax{}.pdf'.format(cc_dir, ncluster)
                        # recurtype_pal = sns.light_palette('navy', n_colors=len(the_df['Recurrence_Type'].unique()),
                        #                                   reverse=True)
                        # ClustermapPlot(the_df, ['Recurrence_Type'], fig_fname, var_title=['Recurrence Type'],
                        #                var_color_pal=[recurtype_pal])

