#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 5/23/17

@author: shuang
@goal: a collection of statistical tools used for data explore and other analysis

"""

import pandas as pd
import scipy.stats as stats
import numpy as np
import statsmodels.stats.multitest as smm
import os


def proportion_table(the_df,var1,var2, table_margin=False):

    df_tmp = the_df.ix[:,[var1,var2]]
    df_tmp = df_tmp.dropna()
    print(df_tmp[var1].value_counts())

    the_observed_tab0 = pd.crosstab(df_tmp[var1], df_tmp[var2], margins=table_margin)
    # the_observed_tab0.columns = df_tmp[var2].unique().tolist() + ['row_totals']
    # the_observed_tab0.index = df_tmp[var1].unique().tolist() + ['col_totals']
    print the_observed_tab0
    print '\n'

    the_observed_tab1 = pd.crosstab(df_tmp[var1], df_tmp[var2], margins=table_margin,normalize='index')
    # the_observed_tab1.columns = df_tmp[var2].unique().tolist() + ['row_totals']
    # the_observed_tab1.index = df_tmp[var1].unique().tolist() + ['col_totals']
    print the_observed_tab1
    print '\n'

    the_observed_tab2 = pd.crosstab(df_tmp[var1], df_tmp[var2], margins=table_margin, normalize='columns')
    # the_observed_tab2.columns = df_tmp[var2].unique().tolist() + ['row_totals']
    # the_observed_tab2.index = df_tmp[var1].unique().tolist() + ['col_totals']
    print the_observed_tab2
    print '\n'
    return (the_observed_tab0, the_observed_tab1, the_observed_tab2)

def chi2test(the_df,var1,var2):
    # print var1, var2
    df_tmp = the_df.ix[:,[var1,var2]]
    df_tmp = df_tmp.dropna()
    # print df_tmp
    N_sample = df_tmp.shape[0]
    the_observed_tab = pd.crosstab(df_tmp[var1], df_tmp[var2], margins=True)
    the_observed_tab.columns = df_tmp[var2].unique().tolist() + ['row_totals']
    the_observed_tab.index = df_tmp[var1].unique().tolist() + ['col_totals']
    # print(the_observed_tab)

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
    # Effective size: calculate the crammer's V for degree of association or a measure of correlation for categorical data
    # source: https://stackoverflow.com/questions/20892799/using-pandas-calculate-cram%C3%A9rs-coefficient-matrix
    # rule of thumb for evaluating the strength of association
    # < .10 = weak, .11 – .30 = moderate, > .31 = strong


    n = observed.sum().sum()
    phi2 = chi2/n
    r,k = observed.shape
    phi2corr = max(0,phi2 - ((k-1)*(r-1))/(n-1))
    rcorr = r - ((r-1)**2)/(n-1)
    kcorr = k - ((k-1)**2)/(n-1)
    crammersV = np.sqrt(phi2corr / min( (kcorr - 1), (rcorr - 1)))

    # print chi2,pval,dof,expected
    return chi2, pval, dof, expected,N_sample, crammersV, the_observed_tab
    # return chi2, pval, dof, expected, N_sample

def apply_mult_corr1(the_df, save_dir, file_tag):
    alpha = 0.05
    method = ['b', 'sidak', 'holm', 'hommel', 'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky']

    the_ov = the_df['outcome'].unique()
    for vv in the_ov:
        pval_raw = the_df.ix[the_df['outcome'] == vv, 'pval'].as_matrix().flatten()
        df_corr_pval = the_df.ix[the_df['outcome'] == vv, ['feature', 'pval']]
        for mm in method:
            rej, pval_corr = smm.multipletests(pval_raw, alpha=alpha, method=mm)[:2]
            df_corr_pval[mm] = pval_corr
            df_corr_pval['{}_rej'.format(mm)] = rej

        # sort the df based on pval
        df_corr_pval.sort_values('pval', inplace=True)
        # save the corrected pvals
        new_outcome_name = '_'.join(vv.split(' '))
        corr_pval_fname = '{}/{}_multitestCorr_{}.csv'.format(save_dir, file_tag, new_outcome_name)
        df_corr_pval.to_csv(corr_pval_fname, index=False)

def apply_mult_corr2(fname, corr_method='fdr_bh', alpha=0.05):

    the_df = pd.read_csv(fname)
    the_ov = the_df['outcome'].unique()
    corr_pval_colname = 'pval_corr_{}'.format(corr_method)
    pval_rej_colname = '{}_rej'.format(corr_method)
    for vv in the_ov:
        pval_raw = the_df.ix[the_df['outcome'] == vv, 'pval'].as_matrix().flatten()
        rej, pval_corr = smm.multipletests(pval_raw, alpha=alpha, method=corr_method)[:2]

        the_df.loc[the_df['outcome'] == vv, corr_pval_colname] = pval_corr
        the_df.loc[the_df['outcome'] == vv, pval_rej_colname] = rej

    # save the corrected pvals
    hh, tt = os.path.splitext(fname)
    corr_pval_fname = '{}_pvalcorr{}'.format(hh, tt)
    the_df.to_csv(corr_pval_fname, index=False)

