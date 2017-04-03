# -*- coding: utf-8 -*-
"""
Created on 2/21/17 12:24 PM

@author: shuang Shih-ying Huang
@goal: perform chi2 test for the association of radiomics features/clusters with clinical features

"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import itertools as itt
import numpy as np
import scipy.stats as stats

def chi2test(the_df,var1,var2):
    the_observed_tab = pd.crosstab(the_df[var1], the_df[var2], margins=True)

    the_observed_tab.columns = the_df[var2].unique().tolist() + ['row_totals']
    the_observed_tab.index = the_df[var1].unique().tolist() + ['col_totals']

    observed = the_observed_tab.ix[0:-1, 0:-1]

    # manual way to get chi2 test result
    # expected = np.outer(the_observed_tab["row_totals"][0:-1],
    #                     the_observed_tab.ix["col_totals"][0:-1]) / float(the_observed_tab['row_totals'][-1])
    #
    # expected = pd.DataFrame(expected)
    # expected.columns = the_df[var2].unique().tolist()
    # expected.index = the_df[var1].unique().tolist()
    #
    # chi_squared_stat = (((observed - expected) ** 2) / expected).sum().sum()
    # dof = (observed.shape[0] - 1) * (observed.shape[1] - 1)
    # crit = stats.chi2.ppf(q=0.95,  # Find the critical value for 95% confidence*
    #                       df=dof)
    # p_value = 1 - stats.chi2.cdf(x=chi_squared_stat,  # Find the p-value
    #                              df=dof)


    chi2,pval,dof,expected = stats.chi2_contingency(observed=observed)
    # chi2_contingency may have different result as the above calcuation when dof=1
    # since the default setting for 'correction' is true ==> Yate's correction for continuity is applied
    # adjust each observed values by 0.5 towards the corresponding expected value

    return chi2, pval, dof, expected


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
        tstage_list = ['T0','Tis','T1mic','T1','T1a','T1b','T1c','T2','T3','T4','T4a','T4b','T4c','T4d','TX','calc']
        tstage_dict = dict(zip(tstage_list,range(len(tstage_list))))
        tstage_dict['T2yp'] = 7  # this is the exception for T2, not sure what T2yp is but assume it's the same as T2
        if x in tstage_dict.keys():
            return tstage_dict[x]
        elif x.lower() == 'calc':
            return tstage_dict['calc']
        else:
            return np.nan
    her2_outcome_df['Tstage'] = her2_outcome_df['T-stage at surgery'].map(Tstage_func)

    def Nstage_func(x):
        nstage_list = ['N0','N1','N2a','N2b','N3a','N3b','N3c','NX','N3a','calc']
        nstage_dict = dict(zip(nstage_list, range(len(nstage_list))))
        if x in nstage_dict.keys():
            return nstage_dict[x]
        elif x.lower() == 'calc':
            return nstage_dict['calc']
        else:
            return np.nan
    her2_outcome_df['Nstage'] = her2_outcome_df['N-stage at surgery'].map(Nstage_func)

    the_outcome_df = her2_outcome_df.loc[:, ['MRN', 'Recurrence_CR', 'Recurrence_Type','Nstage','Tstage']]

    # PET image feature data
    pet_dir = '{}/her2_Analysis/PET/IsoVoxel_IMGBIN128_GLCMBIN64'.format(rootdir)
    df_fname = '{}/PETdataAll_glcmNbin64_normNbin128.csv'.format(pet_dir)
    data_df = pd.read_csv(df_fname, dtype={'pt_mrn': str, 'PRIMARY_ID': str, 'MRN': str, 'Anon_Accession_Num': str})
    jdf = pd.merge(data_df,the_outcome_df,on='MRN')
    print data_df.shape, jdf.shape

    dist_method_list = ['euclidean','spearman','pearson']
    cluster_method_dict = {'pam': None, 'km': None, 'kmdist': None,'hc':['centroid','median','average','mcquitty','complete']}
    # dist_method_list = ['pearson']
    # cluster_method_dict = {'hc': ['complete']}
    theKmax_list = [2,3]
    cluster_name = 'cs_class'
    outcome_name_list = ['Recurrence_Type','Recurrence_CR','TripleNeg','Sjoerd_Grade','Marjan_Histology']

    chi2_all_df = pd.DataFrame(columns=('cluster_method','cluster_linkage','dist_method','Ncluster','N_mincluster','v1','v2','chi2','pval'))
    count = 0
    for ncluster in theKmax_list:
        for k, val in cluster_method_dict.items():
            if val:
                for ll,dd in itt.product(val,dist_method_list):
                    cc_dir = '{}/ConsensusCluster_{}_{}_{}'.format(pet_dir,k,ll,dd)
                    cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir,ncluster)
                    cs_class_df = pd.read_csv(cs_class_fname)
                    cs_class_df.columns = ['ptid_side', 'cs_class']

                    # combine the cs_class to the_df
                    the_df = pd.merge(jdf, cs_class_df, how='left', on='ptid_side')
                    print the_df.shape

                    # find the N of the smallest cluster in each cluster result for future analysis
                    N_mincluster = the_df['cs_class'].value_counts().min()

                    for v2 in outcome_name_list:
                        chi2,pval,dum1,dum2 = chi2test(the_df,cluster_name,v2)
                        chi2_all_df.loc[count] = [k,ll,dd,ncluster,N_mincluster,cluster_name,v2,chi2,pval]
                        count = count + 1

            else:
                for dd in dist_method_list:
                    cc_dir = '{}/ConsensusCluster_{}_{}'.format(pet_dir,k,dd)
                    cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir, ncluster)
                    cs_class_df = pd.read_csv(cs_class_fname)
                    cs_class_df.columns = ['ptid_side', 'cs_class']

                    # combine the cs_class to the_df
                    the_df = pd.merge(jdf, cs_class_df, how='left', on='ptid_side')
                    print the_df.shape

                    # find the N of the smallest cluster in each cluster result for future analysis
                    N_mincluster = the_df['cs_class'].value_counts().min()

                    for v2 in outcome_name_list:
                        chi2,pval,dum1,dum2 = chi2test(the_df,cluster_name,v2)
                        chi2_all_df.loc[count] = [k,None,dd,ncluster,N_mincluster,cluster_name,v2,chi2,pval]
                        count = count + 1

    # save the chi2 relevance test result to .csv
    fname = '{}/chi2_relevancetest_all.csv'.format(pet_dir)
    chi2_all_df.to_csv(fname)

    # print out the parameters with the smallest pval for each outcome variable

    # avoid the uneven classification result (like N = 1 or 2 in a class...)
    # find the N of the smallest cluster
    fdf = chi2_all_df[chi2_all_df['N_mincluster'] > 5]

    idx = fdf.groupby(['Ncluster','v2']).apply(lambda df: df.pval.argmin())
    print fdf.ix[idx,:]



    # # MRI image feature data
    # theTP = [1, 2, 3]
    # theglcmBin = [64, 128, 256]
    # theKmax = 5
    # dist_method_list = ['euclidean','spearman','pearson']
    # cluster_method_dict = {'pam': None, 'km': None, 'kmdist': None,'hc':['centroid','median','average','mcquitty','complete']}
    # # dist_method_list = ['euclidean']
    # # cluster_method_dict = {'km': None}
    #
    # for tt, bb in itt.product(theTP, theglcmBin):
    #     mri_dir = '{}/her2_Analysis/MRI/IsoVoxel_TP{}_GLCMBIN{}'.format(rootdir, tt, bb)
    #     df_fname = '{}/MRIdataAll_tp{}_Nbin{}.csv'.format(mri_dir, tt, bb)
    #     data_df = pd.read_csv(df_fname, dtype={'PRIMARY_ID': str, 'MRN': str})
    #     jdf = pd.merge(data_df, the_outcome_df, on='MRN')
    #
    #     for ncluster in range(2,theKmax+1):
    #         for k, val in cluster_method_dict.items():
    #             if val:
    #                 for ll,dd in itt.product(val,dist_method_list):
    #                     cc_dir = '{}/ConsensusCluster_{}_{}_{}'.format(mri_dir,k,ll,dd)
    #                     cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir,ncluster)
    #                     cs_class_df = pd.read_csv(cs_class_fname)
    #                     cs_class_df.columns = ['ptid_side', 'cs_class']
    #
    #                     # combine the cs_class to the_df
    #                     the_df = pd.merge(jdf, cs_class_df, how='left', on='ptid_side')
    #                     fig_title = 'cluster method: {}, linkage: {}, distance method: {}'.format(k,ll,dd)
    #                     fig_fname = '{}/clustermap_tumorstatus_kmax{}.pdf'.format(cc_dir,ncluster)
    #                     ClustermapPlot_TumorStatus(the_df, fig_title, fig_fname)
    #                     fig_fname = '{}/clustermap_recurstatus_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     ClustermapPlot_PtOutcome1(the_df, fig_title, fig_fname)
    #                     fig_fname = '{}/clustermap_recurtype_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     ClustermapPlot_PtOutcome2(the_df, fig_title, fig_fname)
    #             else:
    #                 for dd in dist_method_list:
    #                     cc_dir = '{}/ConsensusCluster_{}_{}'.format(mri_dir,k,dd)
    #                     cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir, ncluster)
    #                     cs_class_df = pd.read_csv(cs_class_fname)
    #                     cs_class_df.columns = ['ptid_side', 'cs_class']
    #
    #                     # combine the cs_class to the_df
    #                     the_df = pd.merge(jdf, cs_class_df, how='left', on='ptid_side')
    #                     fig_title = 'cluster method: {}, distance method: {}'.format(k, dd)
    #                     fig_fname = '{}/clustermap_tumorstatus_kmax{}.pdf'.format(cc_dir,ncluster)
    #                     ClustermapPlot_TumorStatus(the_df, fig_title, fig_fname)
    #                     fig_fname = '{}/clustermap_recurstatus_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     ClustermapPlot_PtOutcome1(the_df, fig_title, fig_fname)
    #                     fig_fname = '{}/clustermap_recurtype_kmax{}.pdf'.format(cc_dir, ncluster)
    #                     ClustermapPlot_PtOutcome2(the_df, fig_title, fig_fname)