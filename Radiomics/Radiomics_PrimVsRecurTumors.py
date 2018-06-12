#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 6/5/18

@author: shuang
@goal: compare the radiomics of the primary tumors vs recurred tumors of the same patient

"""

from bokeh.plotting import figure, show, output_file
import numpy as np
from bokeh.palettes import Spectral3, Spectral4, Spectral5, Spectral6, Spectral7
from collections import OrderedDict
import pandas as pd
import glob
import re
from math import log, sqrt
from bokeh.io import export_png
# from bokeh.io import export_svgs



def rad(mic, a, b):
    return a * np.sqrt(np.log(mic * 1E4)) + b

def Annular_wedge_plot(df1, df3, filetag):
    inner_radius = 90
    outer_radius = 300 - 10
    minr = sqrt(log(.001 * 1E4))
    maxr = sqrt(log(100000 * 1E4))
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

