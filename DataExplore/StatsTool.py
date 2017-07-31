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

def proportion_table(the_df,var1,var2):

    df_tmp = the_df.ix[:,[var1,var2]]
    df_tmp = df_tmp.dropna()

    the_observed_tab0 = pd.crosstab(df_tmp[var1], df_tmp[var2], margins=True)
    # the_observed_tab0.columns = df_tmp[var2].unique().tolist() + ['row_totals']
    # the_observed_tab0.index = df_tmp[var1].unique().tolist() + ['col_totals']
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
    return chi2, pval, dof, expected,N_sample, crammersV
    # return chi2, pval, dof, expected, N_sample