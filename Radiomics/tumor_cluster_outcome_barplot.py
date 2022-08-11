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
the_N_cluster_rep = 10000
# outcome_name_list = ['Recurrence_Type','Recurrence','Tumor_Grade','Tumor_Histology', 'T_stage', 'N_stage', 'Overall_stage','BC_subtype']
outcome_name_list = ['Recurrence_Type','Recurrence','Tumor_Grade','Tumor_Histology_updated', 'T_stage', 'N_stage', 'Overall_stage','BC_subtype','Overall_stage_2']



# the_N_cluster = 2
# outcome_name_list = ['BoneMetsOrNot','TripleNeg'] + more_outcome

# chi2_outcome_fname = '{}/BestChi2_outcome_Ncluster{}.csv'.format(im_dir, the_N_cluster)
# chi2_outcome_fname = '{}/BestMedCCChi2_outcome_Ncluster{}.csv'.format(im_dir, the_N_cluster)
chi2_outcome_fname = '{}/BestCVChi2_outcome_Ncluster{}_Rep{}.csv'.format(im_dir, the_N_cluster, the_N_cluster_rep)
df_chi2_outcome_all = pd.read_csv(chi2_outcome_fname)

cluster_name = 'cs_class'
light_pal = sns.light_palette((210, 90, 60), input="husl")
# the_outcome_vars_list = ['Tumor_Grade', 'T_stage', 'BC_subtype']
# the_color_pals = [light_pal[1:len(light_pal):2], sns.color_palette('Blues', len(df_data_all['T_stage'].unique())), sns.color_palette('YlOrRd', len(df_data_all['BC_subtype'].unique()))]

the_outcome_vars_list = ['Tumor_Grade', 'Overall_stage', 'BC_subtype', 'Recurrence']
fig_xlabel = ['Tumor Grade', 'Overall stage', 'Breast cancer subtype', 'Disease recurrence']
leg_label = [['T1','T2','T3'], ['S0','S1','S2','S3','S4'],['HR+/HER2-','HR+/HER2+','HR-/HER2+','HR-/HER2-'],['No recur','Recur', 'Never disease-free']]
the_color_pals = [light_pal[1:len(light_pal):2],
                  sns.color_palette('Purples', len(df_data_all['Overall_stage'].unique())),
                  sns.color_palette('Blues', len(df_data_all['BC_subtype'].unique())),
                  sns.color_palette('YlOrRd', len(df_data_all['Recurrence'].unique()))]
cc_prefix = 'patientConsensusCluster'

for ov, str_xlabel, lg_label, pal in zip(the_outcome_vars_list, fig_xlabel,leg_label, the_color_pals):
    print(ov)
    cm, cl, dm = df_chi2_outcome_all.ix[df_chi2_outcome_all['outcome_var'] == ov, ['cluster_method','cluster_linkage','dist_method']].values.flatten()
    if isinstance(cl, str):
        cc_dir = '{}/{}_Rep{}_{}_{}_{}'.format(im_dir, cc_prefix, the_N_cluster_rep, cm,cl,dm)
    else:
        cc_dir = '{}/{}_Rep{}_{}_{}'.format(im_dir, cc_prefix, the_N_cluster_rep, cm, dm)

    cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir, the_N_cluster)
    cs_class_df = pd.read_csv(cs_class_fname)
    cs_class_df.columns = ['ptid_side', 'cs_class']

    # combine the cs_class to the_df
    the_df = pd.merge(df_data_all, cs_class_df, on='ptid_side')

    # print proportion table for chi2 test
    # t1, t2, t3 = st.proportion_table(the_df, 'cs_class', ov)
    #
    # # flatten the proportion table
    # t2_f = t2.unstack().reset_index()
    # t2_f.rename(columns={0: 'N_count'}, inplace=True)
    # t2_f['N_count'] = t2_f['N_count'].apply(lambda x: x*100.)

    # t1_f = t1.unstack().reset_index()
    # t1_f.rename(columns={0: 'N_count'}, inplace=True)

    #  plot the bar charts
    sns.set(style="whitegrid", font='Arial')
    # # Draw a factor plot
    # g = sns.factorplot(x='cs_class', hue=ov, data=the_df,
    #                    size=3, kind='count', aspect=1.8, palette=pal, legend=False, edgecolor='none')
    # g.despine(left=True)
    # g.set_xlabels('Tumor cluster class')
    # g.set_ylabels('Frequency (%)')

    # normalize each cluster group before plotting
    sns.set_context('paper', font_scale=1.2)

    # plot the two proportion table into bar charts
    # normalize to total # of each cluster class
    plt.figure()
    gp = the_df.groupby('cs_class', sort=False)
    cts = gp[ov].value_counts(normalize=True, sort=False)
    dict_tmp = [{'cs_class': cs, ov: outcome, 'percentage': perc * 100.} for (cs, outcome), perc in dict(cts).items()]
    df_tmp = pd.DataFrame(dict_tmp)
    g = sns.barplot(x='cs_class', y='percentage', hue=ov, data = df_tmp, palette = pal, edgecolor='none')
    g.legend_.remove()
    g.set_xlabel('Tumor cluster class', fontsize=13, weight='bold')
    g.set_ylabel('Frequency (%)', fontsize=13, weight='bold')
    sns.despine(left=True)
    plt.ylim(0,90)

    # have to set frameon=True since seaborn turn the frame off automatically
    leg = plt.legend(title=str_xlabel, loc='best', frameon=True)
    frame = leg.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('None')
    for ii in range(len(leg.get_texts())):
        leg.get_texts()[ii].set_text(lg_label[ii])
    plt.tight_layout()

    # if ov == 'BC_subtype':
    #     plt.legend(title=str_xlabel, loc='best')
    # else:
    #     plt.legend(title=str_xlabel, loc='best', ncol=len(the_df[ov].unique()) + 1)

    # figname = '{}/proportional_table_{}.pdf'.format(im_dir, ov)
    figname = '{}/bestCVChi2_proportional_table_{}_freqwrtCS.pdf'.format(im_dir, ov)
    # g.savefig(figname)
    plt.savefig(figname)
    plt.close()

    # normalize to the total # of each outcome variable
    plt.figure()
    gp = the_df.groupby(ov, sort=False)
    cts = gp['cs_class'].value_counts(normalize=True, sort=False)
    dict_tmp = [{'cs_class': cs, ov: outcome, 'percentage': perc * 100.} for (outcome, cs), perc in dict(cts).items()]
    df_tmp = pd.DataFrame(dict_tmp)
    g = sns.barplot(x=ov, y='percentage', hue='cs_class', data=df_tmp, palette=pal, edgecolor='none')
    # g.set(xlabel=str_xlabel, ylabel='Frequency (%)',xticklabels=lg_label)
    g.set_xlabel(str_xlabel, fontsize=13, weight='bold')
    g.set_ylabel('Frequency (%)', fontsize=13, weight='bold')
    g.set_xticklabels(lg_label)
    g.legend_.remove()
    sns.despine(left=True)
    plt.ylim(0, 90)

    # have to set frameon=True since seaborn turn the frame off automatically
    leg = plt.legend(title='Tumor cluster class', loc='best', frameon=True)
    frame = leg.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('None')
    plt.tight_layout()

    figname = '{}/bestCVChi2_proportional_table_{}_freqwrtOutcome.pdf'.format(im_dir, ov)
    plt.savefig(figname)
    plt.clf()



