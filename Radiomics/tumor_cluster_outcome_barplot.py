#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 9/26/17

@author: shuang
@goal: plot the proportion table of the tumor clusters vs outcome

"""

import pandas as pd
import StatsTool as st
import seaborn as sns
import matplotlib.pyplot as plt

rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'

the_mri_tp = 2
the_pet_binwidth = 0.1
the_mri_binwidth = 5
im_dir = '{}/her2_Analysis/PETMRI/PETbinwidth{}_MRItp{}_binwidth{}'.format(rootdir,the_pet_binwidth, the_mri_tp, the_mri_binwidth)


fname = '{}/data_all.csv'.format(im_dir)
df_data_all = pd.read_csv(fname)

# look at the data for clustering algm settings
the_N_cluster = 3
# outcome_name_list = ['Recurrence_Type','Recurrence','Tumor_Grade','Tumor_Histology', 'T_stage', 'N_stage', 'Overall_stage','BC_subtype']
outcome_name_list = ['Recurrence_Type','Recurrence','Tumor_Grade','Tumor_Histology_updated', 'T_stage', 'N_stage', 'Overall_stage','BC_subtype','Overall_stage_2']

# the_N_cluster = 2
# outcome_name_list = ['BoneMetsOrNot','TripleNeg'] + more_outcome

chi2_outcome_fname = '{}/BestChi2_outcome_Ncluster{}.csv'.format(im_dir, the_N_cluster)
df_chi2_outcome_all = pd.read_csv(chi2_outcome_fname)

cluster_name = 'cs_class'
the_outcome_vars_list = ['Tumor_Grade', 'T_stage', 'BC_subtype']
light_pal = sns.light_palette((210, 90, 60), input="husl")
the_color_pals = [light_pal[1:len(light_pal):2], sns.color_palette('Blues', len(df_data_all['T_stage'].unique())), sns.color_palette('YlOrRd', len(df_data_all['BC_subtype'].unique()))]
cc_prefix = 'patientConsensusCluster'

for ov, pal in zip(the_outcome_vars_list, the_color_pals):
    cm, cl, dm = df_chi2_outcome_all.ix[df_chi2_outcome_all['outcome_var'] == ov, ['cluster_method','cluster_linkage','dist_method']].values.flatten()
    if isinstance(cl, str):
        cc_dir = '{}/{}_{}_{}_{}'.format(im_dir, cc_prefix, cm,cl,dm)
    else:
        cc_dir = '{}/{}_{}_{}'.format(im_dir, cc_prefix, cm, dm)

    cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir, the_N_cluster)
    cs_class_df = pd.read_csv(cs_class_fname)
    cs_class_df.columns = ['ptid_side', 'cs_class']

    # combine the cs_class to the_df
    the_df = pd.merge(df_data_all, cs_class_df, on='ptid_side')

    # print proportion table for chi2 test
    t1, t2, t3 = st.proportion_table(the_df, 'cs_class', ov)

    # flatten the proportion table
    t2_f = t2.unstack().reset_index()
    t2_f.rename(columns={0: 'N_count'}, inplace=True)
    t2_f['N_count'] = t2_f['N_count'].apply(lambda x: x*100.)

    # plot the bar charts
    sns.set(style="whitegrid", font='futura std')

    # Draw a nested barplot to show survival for class and sex
    g = sns.factorplot(x='cs_class', hue=ov, data=the_df,
                       size=3, kind='count', aspect=1.8, palette=pal, legend=False, edgecolor='none')
    g.despine(left=True)
    if ov == 'BC_subtype':
        plt.legend(title=ov, loc='best')
    else:
        plt.legend(title=ov, loc='best', ncol=len(the_df[ov].unique()) + 1)
    g.set_xlabels('Tumor cluster class')
    g.set_ylabels('Number of count')

    figname = '{}/proportional_table_{}.pdf'.format(im_dir, ov)
    g.savefig(figname)
