#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 9/25/17

@author: shuang
@goal: a collection of tools used for training and validating a model and other related tasks

"""

import operator
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, f1_score, auc, roc_curve
from itertools import product
import pandas as pd
import numpy as np
import seaborn as sn
import matplotlib.pyplot as plt
import glob


def nested_CV(X, y, clf, params_dict, n_fold, n_trials, feat_names, feat_tag, coef_thresh, im_dir, clf_name, outcome_name):

    items = sorted(params_dict.items())
    param_keys, param_values = zip(*items)
    list_param_values = list(product(*param_values))

    arry_feat_names = np.array(feat_names)
    mean_fpr = np.linspace(0, 1, 100)

    unbias_score_arry = np.zeros(n_trials)
    lst_data_all = []
    data_col_names = ['trial_n', 'fold_n', 'N_sample', 'best_params', 'auc', 'feat_name', 'feat_importance']

    for nn in range(n_trials):
        print('n_trial #: {}'.format(nn))
        skf_outer = StratifiedKFold(n_splits=n_fold, shuffle=True)
        skf_outer.get_n_splits(X, y)

        outer_scores = []
        tprs = []
        # plt.figure()
        for fold, (train_index_outer, test_index_outer) in enumerate(skf_outer.split(X, y)):
            X_train_outer, X_test_outer = X[train_index_outer], X[test_index_outer]
            y_train_outer, y_test_outer = y[train_index_outer], y[test_index_outer]

            inner_mean_scores = []
            for ii, vv in enumerate(list_param_values):
                the_param = dict(zip(param_keys, vv))
                inner_scores = []
                # inner cross-validation
                skf_inner = StratifiedKFold(n_splits=n_fold, shuffle=True)
                for fold_inner, (train_index_inner, test_index_inner) in enumerate(skf_inner.split(X_train_outer, y_train_outer)):
                    # split the training data of outer CV
                    X_train_inner, X_test_inner = X_train_outer[train_index_inner], X_train_outer[test_index_inner]
                    y_train_inner, y_test_inner = y_train_outer[train_index_inner], y_train_outer[test_index_inner]
                    clf.set_params(**the_param)
                    clf.fit(X_train_inner, y_train_inner)

                    # score = f1_score(y_test_inner, clf.predict(X_test_inner))
                    score = roc_auc_score(y_test_inner, clf.predict_proba(X_test_inner)[:, 1])
                    # score = clf.score(X_test_inner, y_test_inner)
                    inner_scores.append(score)
                    # print('inner fold: {}, index: {}, param: {}, score: {}'.format(fold_inner+1, ii, vv, score))

                # calculate mean score for inner folds
                inner_mean_scores.append(np.mean(inner_scores))

            index, value = max(enumerate(inner_mean_scores), key=operator.itemgetter(1))

            # determine the optimal params
            the_opt_param = dict(zip(param_keys, list_param_values[index]))
            print('Best parameter of {} fold: {}'.format(fold + 1, the_opt_param))
            clf.set_params(**the_opt_param)
            clf.fit(X_train_outer, y_train_outer)

            # feature_importance_ or coef_ from an estimator
            feat_importance = getattr(clf, 'feature_importances_', None)
            if feat_importance is None and hasattr(clf, 'coef_'):
                feat_importance = clf.coef_
            # print(feat_importance.shape)

            if fold == 0:
                arry_coeff = np.zeros((feat_importance.flatten().shape[0], n_fold))
            arry_coeff[:, fold] = feat_importance.flatten()

            # plot ROC curves and compute AUCs
            probas_ = clf.predict_proba(X_test_outer)
            fpr, tpr, thresholds = roc_curve(y_test_outer, probas_[:, 1])
            tprs.append(np.interp(mean_fpr, fpr, tpr))
            tprs[-1][0] = 0.0
            score = auc(fpr, tpr)
            # plt.plot(fpr, tpr, lw=1, alpha=0.3, label='ROC fold %d (AUC = %0.2f)' % (fold, score))
            # score = clf.score(X_test_outer, y_test_outer)
            tmp = dict(zip(data_col_names, [nn, fold, len(X_train_outer), the_opt_param, score, feat_names, feat_importance.flatten()]))
            lst_data_all.append(tmp)
            # score = f1_score(y_test_outer, clf.predict(X_test_outer))
            # score = roc_auc_score(y_test_outer, clf.predict_proba(X_test_outer)[:,1])
            outer_scores.append(score)

        # show the prediction error estimate produced by nested CV
        print('unbiased prediction score: {}'.format(outer_scores))

        # # average all feature coeff's
        # avg_coeff = np.mean(arry_coeff, axis=1)
        # std_coeff = np.std(arry_coeff, axis=1)
        # idx = np.abs(avg_coeff) > coef_thresh
        # print('Trial # {}:, features with |coeff| > {}: {}'.format(nn, coef_thresh, arry_feat_names[idx]))
        # print('mean coeff: {}, std coeff: {}'.format(avg_coeff[idx], std_coeff[idx]))
        # unbias_score_arry[nn] = np.mean(outer_scores)
        #
        # # plot mean ROCs
        # mean_tpr = np.mean(tprs, axis=0)
        # mean_tpr[-1] = 1.0
        # mean_auc = auc(mean_fpr, mean_tpr)
        # std_auc = np.std(outer_scores)
        # #     plt.plot(mean_fpr, mean_tpr, color='b',
        # #              label=r'Mean ROC (Trial #{}, AUC = {:.2f} $\pm$ {:.2f}'.format(nn, mean_auc, std_auc),
        # #              lw=2, alpha=.8)
        # plt.plot(mean_fpr, mean_tpr, color='b',
        #          label=r'Mean ROC (AUC = {:.2f} $\pm$ {:.2f}'.format(mean_auc, std_auc),
        #          lw=2, alpha=.8)
        #
        # std_tpr = np.std(tprs, axis=0)
        # tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
        # tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
        # #     plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2)
        # plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
        #                  label=r'$\pm$ 1 std. dev.')
        #
        # plt.xlim([-0.05, 1.05])
        # plt.ylim([-0.05, 1.05])
        # plt.xlabel('False Positive Rate')
        # plt.ylabel('True Positive Rate')
        # # plt.title(fig_title)
        # plt.legend(loc="lower right")
        # fig_name = '{}/Learner/{}_IDV{}_DV{}_Trial{}.pdf'.format(im_dir, clf_name, feat_tag, outcome_name, nn)
        # plt.savefig(fig_name)
        # plt.close()

    # print(unbias_score_arry)
    df_data_all = pd.DataFrame(lst_data_all)
    df_data_all = df_data_all.ix[:, data_col_names]
    fname = '{}/Learner/{}_IDV{}_DV{}_Trial{}_{}folds.json'.format(im_dir, clf_name, feat_tag, outcome_name, n_trials, n_fold)
    df_data_all.to_json(fname)

def CombineFiles(data_dir, clf_name, outcome_var, feat_name, n_trial, k_fold):
    """
    combine all Run's of CV trials and return a final dataframe with the adjusted trial_n (0, ..., n_run*n_trial)
    :param data_dir:
    :param clf_name:
    :param outcome_var:
    :param feat_name:
    :param n_trial:
    :param k_fold:
    :return: the final concatenated dataframe
    """

    lst_fnames = glob.glob('{}/{}_IDV{}_DV{}_Trial{}_Run*_{}folds.json'.format(data_dir, clf_name, feat_name, outcome_var, n_trial, k_fold))
    lst_dfs = []
    for ii in range(len(lst_fnames)):
        df_tmp = pd.read_json(lst_fnames[ii])
        if ii > 0:
            df_tmp['trial_n'] = df_tmp['trial_n'] + ii*n_trial
        lst_dfs = lst_dfs + [df_tmp]
    the_df = pd.concat(lst_dfs, ignore_index=True)

    return the_df
