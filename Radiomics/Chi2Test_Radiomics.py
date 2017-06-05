# -*- coding: utf-8 -*-
"""
Created on 2/21/17 12:24 PM

@author: shuang Shih-ying Huang
@goal: perform chi2 test for the association of radiomics features/clusters with clinical features

"""

import pandas as pd
import re
import itertools as itt
import numpy as np
import scipy.stats as stats

def chi2test(the_df,var1,var2):
    # print var1, var2
    df_tmp = the_df.ix[:,[var1,var2]]
    df_tmp = df_tmp.dropna()
    # print df_tmp
    N_sample = df_tmp.shape[0]
    the_observed_tab = pd.crosstab(df_tmp[var1], df_tmp[var2], margins=True)
    the_observed_tab.columns = df_tmp[var2].unique().tolist() + ['row_totals']
    the_observed_tab.index = df_tmp[var1].unique().tolist() + ['col_totals']

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
    # chi2_contingency may have different result as the above calculation when dof=1
    # since the default setting for 'correction' is true ==> Yate's correction for continuity is applied
    # adjust each observed values by 0.5 towards the corresponding expected value

    # TODO: need to double check the output with R or other stats software
    # # calculate the crammer's V for degree of association or a measure of correlation for categorical data
    # n = observed.sum()
    # phi2 = chi2/n
    # r,k = observed.shape
    # phi2corr = max(0,phi2 - ((k-1)*(r-1))/(n-1))
    # rcorr = r - ((r-1)**2)/(n-1)
    # kcorr = k - ((k-1)**2)/(n-1)
    # crammersV = np.sqrt(phi2corr / min( (kcorr - 1), (rcorr - 1)))

    # print chi2,pval,dof,expected
    # return chi2, pval, dof, expected,N_sample, crammersV
    return chi2, pval, dof, expected, N_sample

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

def recur_dT(row):
    dT = pd.to_datetime(row['Last contact/death']) - pd.to_datetime(row['Date of Primary Dx'])
    if isinstance(dT, pd.tslib.NaTType):
        dT_yr = -1
    else:
        dT_yr = dT.days / 365.

    return dT_yr

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

# def recur_type_func(x):
#     if x is np.nan or x.lower() == 'unknown' or x.lower() == 'never disease free':
#         return np.nan
#     elif x == 'NONE/DISEASE FREE':
#         return 0
#     elif x == 'DIST RECUR, BONE':
#         return 1
#     elif re.search(r'LOCAL RECUR[\w\s.]*',x):
#         return 2
#     elif re.search(r'DIST RECUR[\w\s.]*',x) and x != 'DIST RECUR, BONE':
#         return 3

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


    # Recurrence_CR (recurrence info based on Cancer Registry field 'Recurrence Type Summary')
    # her2_outcome_df['Recurrence_CR_0'] = her2_outcome_df['Recurrence Type Summary'].map(recur_func0)
    her2_outcome_df['Recurrence_CR_1'] = her2_outcome_df['Recurrence Type Summary'].map(recur_func1)
    her2_outcome_df['Recurrence_CR_2'] = her2_outcome_df['Recurrence Type Summary'].map(recur_func2)
    her2_outcome_df['Diseasefree_5yr'] = her2_outcome_df.apply(diseasefree5yr_func,axis=1)
    her2_outcome_df['Recurrence_Type'] = her2_outcome_df['Recurrence Type Summary'].map(recur_type_func)

    # catergorize tumor and lymph node stages
    her2_outcome_df['T_stage'] = her2_outcome_df['T-stage at surgery'].map(Tstage_func)
    her2_outcome_df['N_stage'] = her2_outcome_df['N-stage at surgery'].map(Nstage_func)
    her2_outcome_df['Overall_stage'] = her2_outcome_df['Overall stage'].map(Overallstage_func)

    outcome_name_list = ['Recurrence_Type', 'Diseasefree_5yr', 'Recurrence_CR_1', 'Recurrence_CR_2', 'TripleNeg',
                         'Sjoerd_Grade','Marjan_Histology', 'T_stage', 'N_stage', 'Overall_stage','BC_subtype']
    the_outcome_df = her2_outcome_df.loc[:, ['MRN','Diseasefree_5yr','Recurrence_CR_1', 'Recurrence_CR_2', 'Recurrence_Type','N_stage','T_stage','Overall_stage']]
    # the_outcome_df = her2_outcome_df.loc[:,['MRN'] + outcome_name_list]

    # # PET image feature data
    # pet_dir = '{}/her2_Analysis/PET/IsoVoxel_IMGBIN128_GLCMBIN64'.format(rootdir)
    # df_fname = '{}/PETdataAll_glcmNbin64_normNbin128.csv'.format(pet_dir)
    # data_df = pd.read_csv(df_fname, dtype={'pt_mrn': str, 'PRIMARY_ID': str, 'MRN': str, 'Anon_Accession_Num': str})
    # data_df['BC_subtype'] = data_df.apply(BC_subtype_func,axis=1)
    # jdf = pd.merge(data_df,the_outcome_df,on='MRN')
    # # print data_df.shape, jdf.shape
    #
    # dist_method_list = ['euclidean','spearman','pearson']
    # cluster_method_dict = {'pam': None, 'km': None, 'kmdist': None,'hc':['centroid','median','average','mcquitty','complete']}
    # theKmax_list = [2,3,4,5]
    # cluster_name = 'cs_class'
    #
    # # #for testing
    # # dist_method_list = ['pearson']
    # # cluster_method_dict = {'kmdist': None}
    # # theKmax_list = [3]
    # # cluster_name = 'cs_class'
    # # outcome_name_list = ['Recurrence_CR']
    #
    # chi2_all_df = pd.DataFrame(columns=('cluster_method','cluster_linkage','dist_method','Ncluster','N_mincluster','N_sample','v1','v2','chi2','pval'))
    # count = 0
    # for ncluster in theKmax_list:
    #     for k, val in cluster_method_dict.items():
    #         if val:
    #             for ll,dd in itt.product(val,dist_method_list):
    #                 cc_dir = '{}/ConsensusCluster_{}_{}_{}'.format(pet_dir,k,ll,dd)
    #                 cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir,ncluster)
    #                 cs_class_df = pd.read_csv(cs_class_fname)
    #                 cs_class_df.columns = ['ptid_side', 'cs_class']
    #
    #                 # combine the cs_class to the_df
    #                 the_df = pd.merge(jdf, cs_class_df, how='left', on='ptid_side')
    #                 # print the_df.shape
    #
    #                 # find the N of the smallest cluster in each cluster result for future analysis
    #                 N_mincluster = the_df['cs_class'].value_counts().min()
    #
    #                 for v2 in outcome_name_list:
    #                     # print 'v2 is {}'.format(v2)
    #                     chi2,pval,dum1,dum2,N_sample = chi2test(the_df,cluster_name,v2)
    #                     # print chi2,pval, N_sample
    #                     chi2_all_df.loc[count] = [k,ll,dd,ncluster,N_mincluster,N_sample,cluster_name,v2,chi2,pval]
    #                     count = count + 1
    #
    #         else:
    #             for dd in dist_method_list:
    #                 cc_dir = '{}/ConsensusCluster_{}_{}'.format(pet_dir,k,dd)
    #                 cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir, ncluster)
    #                 cs_class_df = pd.read_csv(cs_class_fname)
    #                 cs_class_df.columns = ['ptid_side', 'cs_class']
    #
    #                 # combine the cs_class to the_df
    #                 the_df = pd.merge(jdf, cs_class_df, how='left', on='ptid_side')
    #                 # print the_df.shape
    #
    #                 # find the N of the smallest cluster in each cluster result for future analysis
    #                 N_mincluster = the_df['cs_class'].value_counts().min()
    #
    #                 for v2 in outcome_name_list:
    #                     # print 'v2 is {}'.format(v2)
    #                     chi2,pval,dum1,dum2,N_sample = chi2test(the_df,cluster_name,v2)
    #                     chi2_all_df.loc[count] = [k,None,dd,ncluster,N_mincluster,N_sample,cluster_name,v2,chi2,pval]
    #                     count = count + 1
    #
    # # save the chi2 relevance test result to .csv
    # fname = '{}/chi2_relevancetest_all.csv'.format(pet_dir)
    # chi2_all_df.to_csv(fname)

    # MRI image feature data
    # theTP = [1, 2, 3]
    # theglcmBin = [64, 128, 256]
    theTP = [2]
    theglcmBin = [128]
    dist_method_list = ['euclidean','spearman','pearson']
    cluster_method_dict = {'pam': None, 'km': None, 'kmdist': None,'hc':['centroid','median','average','mcquitty','complete']}
    # dist_method_list = ['euclidean']
    # cluster_method_dict = {'km': None}

    theKmax_list = [2, 3, 4, 5]
    cluster_name = 'cs_class'
    for tt, bb in itt.product(theTP, theglcmBin):
        mri_dir = '{}/her2_Analysis/MRI/IsoVoxel_TP{}_GLCMBIN{}'.format(rootdir, tt, bb)
        df_fname = '{}/MRIdataAll_tp{}_Nbin{}.csv'.format(mri_dir, tt, bb)
        data_df = pd.read_csv(df_fname, dtype={'PRIMARY_ID': str, 'MRN': str})
        data_df['BC_subtype'] = data_df.apply(BC_subtype_func, axis=1)
        jdf = pd.merge(data_df, the_outcome_df, on='MRN')
        count = 0

        chi2_all_df = pd.DataFrame(columns=('cluster_method', 'cluster_linkage', 'dist_method', 'Ncluster', 'N_mincluster', 'Nsample', 'v1', 'v2', 'chi2','pval'))
        for ncluster in theKmax_list:
            for k, val in cluster_method_dict.items():
                if val:
                    for ll,dd in itt.product(val,dist_method_list):
                        cc_dir = '{}/ConsensusCluster_{}_{}_{}'.format(mri_dir, k, ll, dd)
                        cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir,ncluster)
                        cs_class_df = pd.read_csv(cs_class_fname)
                        cs_class_df.columns = ['ptid_side', 'cs_class']

                        # combine the cs_class to the_df
                        the_df = pd.merge(jdf, cs_class_df, how='left', on='ptid_side')
                        # print the_df.shape

                        # find the N of the smallest cluster in each cluster result for future analysis
                        N_mincluster = the_df['cs_class'].value_counts().min()

                        for v2 in outcome_name_list:
                            # print 'v2 is {}'.format(v2)
                            chi2, pval, dum1, dum2, N_sample = chi2test(the_df, cluster_name, v2)
                            chi2_all_df.loc[count] = [k,ll,dd,ncluster,N_mincluster,N_sample,cluster_name,v2,chi2,pval]
                            count = count + 1
                else:
                    for dd in dist_method_list:
                        cc_dir = '{}/ConsensusCluster_{}_{}'.format(mri_dir, k, dd)
                        cs_class_fname = '{}/ConsensusClass_kmax{}.csv'.format(cc_dir, ncluster)
                        cs_class_df = pd.read_csv(cs_class_fname)
                        cs_class_df.columns = ['ptid_side', 'cs_class']

                        # combine the cs_class to the_df
                        the_df = pd.merge(jdf, cs_class_df, how='left', on='ptid_side')
                        # print the_df.shape

                        # find the N of the smallest cluster in each cluster result for future analysis
                        N_mincluster = the_df['cs_class'].value_counts().min()

                        for v2 in outcome_name_list:
                            # print 'v2 is {}'.format(v2)
                            chi2, pval, dum1, dum2, N_sample = chi2test(the_df, cluster_name, v2)
                            chi2_all_df.loc[count] = [k,None,dd,ncluster,N_mincluster,N_sample,cluster_name,v2,chi2,pval]
                            count = count + 1

        # save the chi2 relevance test result to .csv
        fname = '{}/chi2_relevancetest_all.csv'.format(mri_dir)
        chi2_all_df.to_csv(fname)