# -*- coding: utf-8 -*-
"""
Created on 4/18/17 3:18 PM

@author: shuang Shih-ying Huang
@goal: display the proportion table to evaluate the chi2 test result

"""

import pandas as pd
import re
import numpy as np

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

# 0: no recur, 1: local recur, 2: distant recur (bone and others),
# nan or 'unknown' or 'never disease free': N/A,
def recur_type_func(x):
    if x is np.nan or x.lower().find('unknown') >= 0 or x.lower() == 'never disease free':
        return np.nan
    elif x == 'NONE/DISEASE FREE':
        return 0
    elif re.search(r'LOCAL RECUR[\w\s.]*',x):
        return 1
    elif re.search(r'DIST RECUR[\w\s.]*',x):
        return 2

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

def proportion_table(the_df,var1,var2):

    df_tmp = the_df.ix[:,[var1,var2]]
    df_tmp = df_tmp.dropna()

    the_observed_tab0 = pd.crosstab(df_tmp[var1], df_tmp[var2], margins=True)
    the_observed_tab0.columns = df_tmp[var2].unique().tolist() + ['row_totals']
    the_observed_tab0.index = df_tmp[var1].unique().tolist() + ['col_totals']
    print the_observed_tab0
    print '\n'

    the_observed_tab1 = pd.crosstab(df_tmp[var1], df_tmp[var2], margins=True,normalize='index')
    # the_observed_tab1.columns = df_tmp[var2].unique().tolist() + ['row_totals']
    # the_observed_tab1.index = df_tmp[var1].unique().tolist() + ['col_totals']
    print the_observed_tab1
    print '\n'

    the_observed_tab2 = pd.crosstab(df_tmp[var1], df_tmp[var2], margins=True, normalize='columns')
    # the_observed_tab2.columns = df_tmp[var2].unique().tolist() + ['row_totals']
    # the_observed_tab2.index = df_tmp[var1].unique().tolist() + ['col_totals']
    print the_observed_tab2
    print '\n'



if __name__ == '__main__':
    # gather outcome data
    rootdir = '/Users/shuang/Documents/Proj_Radiomics/Data/her2'
    outcome_fname = '{}/her2_ClinicalData/her2_outcome.csv'.format(rootdir)
    her2_outcome_df = pd.read_csv(outcome_fname,dtype={'MRN': str,'Recurrence Type Summary': str})

    # Recurrence_CR (recurrence info based on Cancer Registry field 'Recurrence Type Summary')
    her2_outcome_df['Recurrence_CR_1'] = her2_outcome_df['Recurrence Type Summary'].map(recur_func1)
    her2_outcome_df['Recurrence_CR_2'] = her2_outcome_df['Recurrence Type Summary'].map(recur_func2)

    her2_outcome_df['Diseasefree_5yr'] = her2_outcome_df.apply(diseasefree5yr_func,axis=1)
    her2_outcome_df['Recurrence_Type'] = her2_outcome_df['Recurrence Type Summary'].map(recur_type_func)

    # catergorize tumor and lymph node stages
    her2_outcome_df['T_stage'] = her2_outcome_df['T-stage at surgery'].map(Tstage_func)
    her2_outcome_df['N_stage'] = her2_outcome_df['N-stage at surgery'].map(Nstage_func)
    her2_outcome_df['Overall_stage'] = her2_outcome_df['Overall stage'].map(Overallstage_func)

    the_outcome_df = her2_outcome_df.loc[:, ['MRN','Diseasefree_5yr','Recurrence_CR_1', 'Recurrence_CR_2', 'Recurrence_Type','N_stage','T_stage','Overall_stage']]

    # the_outcome_of_interest = ['Recurrence_Type']
    # the_outcome_of_interest = ['TripleNeg']
    # the_cluster_of_interest = [2.0]
    #
    # the_outcome_of_interest = ['Recurrence_CR_1']
    # the_cluster_of_interest = [3.0]

    # the_outcome_of_interest = ['Diseasefree_5yr']
    # the_cluster_of_interest = [2.0]

    # the_outcome_of_interest = ['Sjoerd_Grade']
    # the_cluster_of_interest = [3.0]

    the_outcome_of_interest = ['T_stage']
    the_cluster_of_interest = [3.0]


    # PET image feature data
    print '\n'
    print 'PET IMAGE FEATURES .....'
    pet_imgbin = 128
    pet_glcmbin = 64
    pet_dir = '{}/her2_Analysis/PET/IsoVoxel_IMGBIN{}_GLCMBIN{}'.format(rootdir,pet_imgbin,pet_glcmbin)
    df_fname = '{}/PETdataAll_glcmNbin{}_normNbin{}.csv'.format(pet_dir, pet_glcmbin, pet_imgbin)
    data_df = pd.read_csv(df_fname, dtype={'pt_mrn': str, 'PRIMARY_ID': str, 'MRN': str, 'Anon_Accession_Num': str})
    jdf = pd.merge(data_df,the_outcome_df,on='MRN')

    chi2_fname = '{}/chi2_relevancetest_all.csv'.format(pet_dir)
    chi2_all_df = pd.read_csv(chi2_fname)

    # avoid the uneven classification result (like N = 1 or 2 in a class...)
    # find the N of the smallest cluster
    fdf = chi2_all_df[chi2_all_df['N_mincluster'] > 5]

    # print out the parameters with the smallest pval for each outcome variable
    idx = fdf.groupby(['Ncluster','v2']).apply(lambda df: df.pval.argmin())

    # # switch the chi2 output for Overall_stage since N_stage and Overstage can be combined
    # the_cm = 'kmdist'
    # the_cl = None
    # the_dd = 'pearson'
    # the_v2 = 'Recurrence_CR_1'
    # the_ncluster = 3

    # the_cm = 'hc'
    # the_cl = 'mcquitty'
    # the_dd = 'spearman'
    # the_v2 = 'Diseasefree_5yr'
    # the_ncluster = 2

    # if the_cl:
    #     df_oh = fdf.ix[(fdf['cluster_method'] == the_cm) & (fdf['cluster_linkage'] == the_cl) & (fdf['dist_method'] == the_dd) & (fdf['v2'] == the_v2) & (fdf['Ncluster'] == the_ncluster)]
    # else:
    #     df_oh = fdf.ix[(fdf['cluster_method'] == the_cm) & (fdf['dist_method'] == the_dd) & (fdf['v2'] == the_v2) & (fdf['Ncluster'] == the_ncluster)]
    # idx[the_ncluster, the_v2] = df_oh.index.get_values()[0]

    # look at the stats for T and N stage
    idx = pd.DataFrame(idx)
    # df_tmp = idx.loc[(slice(None), the_outcome_of_interest), :]
    df_tmp = idx.loc[(the_cluster_of_interest,the_outcome_of_interest), :]
    idx2 = df_tmp[0]
    cluster_name = 'cs_class'
    for ii in range(len(idx2)):
        cm,cl,dm,Nc,v2,chi2,pval = fdf.ix[idx2.iloc[ii],['cluster_method','cluster_linkage','dist_method','Ncluster','v2','chi2','pval']]
        print '{}, Nc = {}, pval = {}'.format(v2,Nc,pval)
        # print fdf.ix[idx2.iloc[ii],['cluster_method','cluster_linkage','dist_method','Ncluster']]
        if cl is not np.nan:
            print 'cluster: {}, {}, {}'.format(cm, cl, dm)
            cc_dir = '{}/ConsensusCluster_{}_{}_{}'.format(pet_dir,cm,cl,dm)
        else:
            print 'cluster: {}, {}'.format(cm, dm)
            cc_dir = '{}/ConsensusCluster_{}_{}'.format(pet_dir,cm,dm)
        cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir, int(Nc))
        cs_class_df = pd.read_csv(cs_class_fname)
        cs_class_df.columns = ['ptid_side', 'cs_class']

        # combine the cs_class to the_df
        the_df = pd.merge(jdf, cs_class_df, how='left', on='ptid_side')
        proportion_table(the_df,cluster_name,v2)

    # MRI image feature
    print '\n'
    print 'MRI IMAGE FEATURES .....'
    mri_tp = 2
    mri_glcmbin = 128
    mri_dir = '{}/her2_Analysis/MRI/IsoVoxel_TP{}_GLCMBIN{}'.format(rootdir, mri_tp, mri_glcmbin)
    df_fname = '{}/MRIdataAll_tp{}_Nbin{}.csv'.format(mri_dir, mri_tp, mri_glcmbin)
    data_df = pd.read_csv(df_fname, dtype={'pt_mrn': str, 'PRIMARY_ID': str, 'MRN': str, 'Anon_Accession_Num': str})
    jdf = pd.merge(data_df, the_outcome_df, on='MRN')

    chi2_fname = '{}/chi2_relevancetest_all.csv'.format(mri_dir)
    chi2_all_df = pd.read_csv(chi2_fname)

    # avoid the uneven classification result (like N = 1 or 2 in a class...)
    # find the N of the smallest cluster
    fdf = chi2_all_df[chi2_all_df['N_mincluster'] > 5]

    # print out the parameters with the smallest pval for each outcome variable
    idx = fdf.groupby(['Ncluster', 'v2']).apply(lambda df: df.pval.argmin())

    # look at the stats for T and N stage
    idx = pd.DataFrame(idx)

    # df_tmp = idx.loc[(slice(None), the_outcome_of_interest), :]
    df_tmp = idx.loc[(the_cluster_of_interest, the_outcome_of_interest), :]
    idx2 = df_tmp[0]
    cluster_name = 'cs_class'
    for ii in range(len(idx2)):
        cm, cl, dm, Nc, v2, chi2, pval = fdf.ix[idx2.iloc[ii], ['cluster_method', 'cluster_linkage', 'dist_method', 'Ncluster', 'v2', 'chi2', 'pval']]
        print '{}, Nc = {}, pval = {}'.format(v2, Nc, pval)
        # print fdf.ix[idx2.iloc[ii],['cluster_method','cluster_linkage','dist_method','Ncluster']]
        if cl is not np.nan:
            print 'cluster: {}, {}, {}'.format(cm, cl, dm)
            cc_dir = '{}/ConsensusCluster_{}_{}_{}'.format(mri_dir, cm, cl, dm)
        else:
            print 'cluster: {}, {}'.format(cm, dm)
            cc_dir = '{}/ConsensusCluster_{}_{}'.format(mri_dir, cm, dm)
        cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir, int(Nc))
        cs_class_df = pd.read_csv(cs_class_fname)
        cs_class_df.columns = ['ptid_side', 'cs_class']

        # combine the cs_class to the_df
        the_df = pd.merge(jdf, cs_class_df, how='left', on='ptid_side')
        proportion_table(the_df, cluster_name, v2)