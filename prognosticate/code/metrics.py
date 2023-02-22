import pandas as pd
import argparse
import numpy as np
from pycox.evaluation import EvalSurv
import glob as glob
import re
import seaborn as sn
import matplotlib.pyplot as plt
import scipy
import matplotlib as mpl
from sksurv.metrics import (concordance_index_censored, concordance_index_ipcw, cumulative_dynamic_auc)
from sksurv.util import Surv
from scipy.stats import wilcoxon
import itertools as it
import shap
import math
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

parser = argparse.ArgumentParser(description='Prediction 1.0 Multiple Myeloma workflow')
parser.add_argument('--in_data_file', type=str, default='/Users/axr2376/Desktop/pred_1_0_paper/data_in', help='input dataset path')
parser.add_argument('--path', type=str, default='/Users/axr2376/Desktop/pred_1_0_paper/data_out/expt_1', help='input dataset path')
parser.add_argument('--legend', type=str, default='/Users/axr2376/Desktop/pred_1_0_paper/data_in/legend_PMMM_2022.xlsx', help='input dataset path')

parser.add_argument('--outcome', type=str, default='efs', help='none in the case of first level comparisons, os, efs, pfs')
parser.add_argument('--get_metric', type=str, default='C_index_harrell', help='options: NBLL, C_index_antolini, C_index_harrell, C_index_uno, Int_Brier_Score, AUC')
parser.add_argument('--preds_time_in_days', type=int, default=1825, help='time in days')
parser.add_argument('--thresh_m_to_p2', type=int, default=365, help='days')

parser.add_argument('--group_feat_add', type=str, default='0~1~2~100', help='0~1~2~3~4~5~6~7~8~100 ID group feature addition')

parser.add_argument('--c_index_state', action='store_true', default=False, help='calculate state c-index of models if True')
parser.add_argument('--c_index_group_addition', action='store_true', default=False, help='calculate multistate c-index of best model if True')
parser.add_argument('--rank_folds', action='store_true', default=False, help='plot ranking across splits for all models')
parser.add_argument('--surv_model', type=str, default='neural_cox_non_prop', help='models')
parser.add_argument('--agg_perm_imp', action='store_true', default=False, help='aggregate and plot')
parser.add_argument('--get_top_gen_feats', action='store_true', default=False, help='get top gen feats across 4 groups')
parser.add_argument('--heatmap_bygroup', action='store_true', default=False, help='heatmap of gen feats across 4 groups')
parser.add_argument('--shap', action='store_true', default=False, help='shap plots')
parser.add_argument('--get_frac_plot', action='store_true', default=False, help='fraction of patients plot using loo preds')
parser.add_argument('--re_rank_feats', action='store_true', default=False, help='re-rank permutation imp feats based on univariate feats')
parser.add_argument('--get_risk_scores', action='store_true', default=False, help='risk scores')
parser.add_argument('--relative_imp', action='store_true', default=False, help='relative imp scores; rank_folds True; groups_0to8 True')
parser.add_argument('--ext_val', action='store_true', default=False, help='external test set')
parser.add_argument('--plot_sub_groups', action='store_true', default=False, help='c-index, kfolds_rank')
parser.add_argument('--test_graph', action='store_true', default=False, help='test graphs')


def data_filter(df, ext_val):
    df = df[df['duration'] >= 0]
    df = df.dropna(subset=['os_time', 'pfs_time', 'pfs_code', 'os_code'])

    df['time_SCT'].fillna(0, inplace=True)

    if ext_val is not True: # not for independent Heidelberg cohort
        df_r_iss = pd.read_csv(args.in_data_file + '/ISS_RISS_R2ISS_LDH_all_cohort.txt', sep='\t')
        df_r_iss['sample'] = np.where(df_r_iss['study'] == 'MRC_XI', 'MRC_' + df_r_iss['sample'], df_r_iss['sample'])
        df_r_iss['sample'] = np.where(df_r_iss['study'] == 'UAMS', 'UAMS_' + df_r_iss['sample'], df_r_iss['sample'])
        df = pd.merge(df, df_r_iss[['sample', 'R_ISS', 'R2_ISS']], on='sample', how='inner')

    if ext_val:
        df['ISS'] = np.select([(df['ISS'] == 'I'), (df['ISS'] == 'II'), (df['ISS'] == 'III')], ['ISS1', 'ISS2', 'ISS3'])

    # map clinical vars
    map_gender = {'male': 0, 'female': 1}
    x_df = df.replace({'gender': map_gender})
    map_study = {'MMRF': 0, 'MPG': 1, 'Moffit': 2, 'MSKCC_292': 3, 'UAMS': 4,  'HD6': 5}
    x_df = x_df.replace({'study': map_study})
    map_ecog = {'ecog<2': 0, 'ecog>=2': 1}
    x_df = x_df.replace({'ecog': map_ecog})
    map_iss = {'ISS1': 0, 'ISS2': 1, 'ISS3': 2}
    x_df = x_df.replace({'ISS': map_iss})
    map_r_iss = {'R-ISS1': 0, 'R-ISS2': 1, 'R-ISS3': 2}
    x_df = x_df.replace({'R_ISS': map_r_iss})
    map_r2_iss = {'I': 0, 'II': 1, 'III': 2, 'IV':3}
    x_df = x_df.replace({'R2_ISS': map_r2_iss})
    x_df["combo"] = np.select([x_df.combo <= 1, x_df.combo > 1], [0, 1])
    x_df["age"] = np.select([(x_df.age < 65), (x_df.age >= 65) & (x_df.age < 75), x_df.age >= 75], [0, 1, 2])

    map_ldh = {'Low': 0, 'High': 1, 'Normal': 2}
    x_df = x_df.replace({'LDH_level': map_ldh})

    map_race = {'WHITE': 0, 'BLACK OR AFRICAN AMERICAN': 1, 'OTHER': 2}
    x_df = x_df.replace({'Race': map_race})

    x_df['pfs_time'] = x_df['pfs_time'].astype(int)
    x_df['os_time'] = x_df['os_time'].astype(int)
    x_df['duration'] = x_df['duration'].astype(int)
    x_df['time_SCT'] = x_df['time_SCT'].astype(int)

    x_df['pfs_time'] = np.where(
        (x_df['duration'] > x_df['pfs_time']) & (x_df['pfs_code'] == 0) & (x_df['os_code'] == 0),
        x_df[['duration', 'pfs_time', 'os_time']].max().max(), x_df['pfs_time'])

    x_df['os_time'] = np.where(
        (x_df['duration'] > x_df['os_time']) & (x_df['pfs_code'] == 0) & (x_df['os_code'] == 0),
        x_df[['duration', 'pfs_time', 'os_time']].max().max(), x_df['os_time'])

    x_df['os_time'] = np.where(
        (x_df['pfs_time'] > x_df['os_time']) & (x_df['pfs_code'] == 1),
        x_df[['pfs_time', 'os_time']].max().max(), x_df['os_time'])

    threshold = args.thresh_m_to_p2
    x_df['duration'] = np.where(x_df['duration'] >= threshold, threshold, x_df['duration'])

    # y_df = x_df[['pfs_code', 'pfs_time', 'os_code', 'os_time']]

    if args.outcome == 'pfs':
        x_df['pfs_time'] = np.select([(x_df['pfs_code'] == 0) & (x_df['os_code'] == 0),
                                  (x_df['pfs_code'] == 0) & (x_df['os_code'] == 1),
                                  (x_df['pfs_code'] == 1) & (x_df['os_code'] == 0),
                                  (x_df['pfs_code'] == 1) & (x_df['os_code'] == 1)
                                  ], [x_df['os_time'], x_df['os_time'], x_df['pfs_time'], x_df['pfs_time']])

        x_df['pfs_code'] = np.where((x_df['pfs_code'] == 1) | (x_df['os_code'] == 1), 1, 0)

    x_df['phases'] = np.where((x_df['phase'] == 'Phase1-only'), 'phase_1', 'phase_2')

    x_df.drop('phase', axis=1, inplace=True)

    x_df['stages'] = np.select([x_df['phases'] == 'phase_1',
                                   (x_df['phases'] == 'phase_2') & (x_df['pfs_time'] <= 540) & ~((x_df['pfs_code']==0) & (x_df['os_code']==1)),
                                   (x_df['phases'] == 'phase_2') & ((x_df['pfs_code']==0) & (x_df['os_code']==1)) & (x_df['os_time'] <= 540),
                                   (x_df['phases'] == 'phase_2') & (x_df['pfs_time'] > 540) & ~((x_df['pfs_code']==0) & (x_df['os_code']==1)),
                                   (x_df['phases'] == 'phase_2') & ((x_df['pfs_code']==0) & (x_df['os_code']==1)) & (x_df['os_time'] > 540)],
                                    ['induction', 'remission', 'remission', 'post-18months', 'post-18months'])

    x_df['sub_stages'] = np.select([(x_df['stages'] == 'induction') & (x_df['pfs_code'] == 0) & (x_df['os_code'] == 0),
                                       (x_df['stages'] == 'induction') & (x_df['pfs_code'] == 0) & (x_df['os_code'] == 1),
                                       (x_df['stages'] == 'induction') & (x_df['pfs_code'] == 1) & (x_df['os_code'] == 0),
                                       (x_df['stages'] == 'induction') & (x_df['pfs_code'] == 1) & (x_df['os_code'] == 1),
                                       (x_df['stages'] == 'remission') & (x_df['pfs_code'] == 0) & (x_df['os_code'] == 0),
                                       (x_df['stages'] == 'remission') & (x_df['pfs_code'] == 0) & (x_df['os_code'] == 1),
                                       (x_df['stages'] == 'remission') & (x_df['pfs_code'] == 1) & (x_df['os_code'] == 0),
                                       (x_df['stages'] == 'remission') & (x_df['pfs_code'] == 1) & (x_df['os_code'] == 1),
                                        (x_df['stages'] == 'post-18months') & (x_df['pfs_code'] == 0) & (x_df['os_code'] == 0),
                                        (x_df['stages'] == 'post-18months') & (x_df['pfs_code'] == 0) & (x_df['os_code'] == 1),
                                       (x_df['stages'] == 'post-18months') & (x_df['pfs_code'] == 1) & (x_df['os_code'] == 0),
                                       (x_df['stages'] == 'post-18months') & (x_df['pfs_code'] == 1) & (x_df['os_code'] == 1)
                                       ],
                ['induction_00', 'induction_01', 'induction_10', 'induction_11',
                 'remission_00', 'remission_01', 'remission_10', 'remission_11',
                 'post-18months_00', 'post-18months_01', 'post-18months_10', 'post-18months_11'])

    df_x_entire = x_df.copy()

    return df_x_entire


def load_dataset(val=False):
    if val is True:
        df = pd.read_csv(args.in_data_file + '/heidelberg_matrix.txt', sep='\t')
    else:
        df = pd.read_csv(args.in_data_file + '/PMMM_matrix_12052022.txt', sep='\t')

    df = df.drop('SCT_line', axis=1)
    df = df[df['duration'] >= 0]

    df = data_filter(df, val)
    print(df.shape)
    return df


def prep_ground_truth(filename, df, df_metr_conc, group):
    print(filename)
    df_preds = pd.read_csv(filename, index_col=None, header=0)
    df_preds = df_preds.rename(columns=lambda x: re.sub('\.','-',x))
    df_test = df[df['sample'].isin(list(df_preds))].sort_index()

    df_preds = df_preds[list(df_test['sample'])].sort_index()

    df_train = df[~df['sample'].isin(df_test['sample'].to_list())].sort_index()

    df_label_preds = df_test[df_test['sample'].isin(list(df_preds))]

    c_ind = None #c_antolini(df_preds, df_label_preds, group)
    auc_score = None ; #auc(df_preds, df_label_preds, df, group)

    # print('ext val', args.ext_val, df_preds.shape, df_label_preds.shape)
    if args.ext_val is not True:
        # print('in')
        int_br_score = brier_score(df_preds, df_label_preds, group)
        nbll_score = bll(df_preds, df_label_preds, group)
        c_ind_uno = c_uno(df_preds, df_label_preds, df_train, group)
    else:
        int_br_score = None
        nbll_score = None
        c_ind_uno = None
    c_ind_harell = c_harrell(df_preds, df_label_preds, group)

    model = filename.split('/')[-1].split('-')[-1].split('.')[0]
    group_id = filename.split('/')[-1].split('-')[1]
    fold = filename.split('/')[-1].split('-')[3]

    df_metr = pd.DataFrame({'C_index_antolini': [c_ind], 'C_index_harrell': [c_ind_harell],
                            'C_index_uno': [c_ind_uno], 'Model': [model], 'Fold' : fold,
                            'Int_Brier_Score': [int_br_score], 'AUC': [auc_score], 'NBLL': nbll_score, 'Group_id': group_id})
    df_metr_conc = pd.concat([df_metr_conc, df_metr], ignore_index=True)

    df_metr_conc.drop_duplicates(inplace=True)
    return df_metr_conc


'''
metrics
'''
def c_antolini(preds, Y, group):
    if group == 'Progress (P1)' or group == 'Progress & deceased (P1)':
        preds = preds.head(args.thresh_m_to_p2-1)
    else:
        preds = preds.head(args.preds_time_in_days-1)

    if args.outcome == 'os' or group == 'Progress & deceased (P1)' or \
                            group == 'Progress & deceased (P2)' or group == 'Non-progress & deceased (P2)':
        get_target = lambda df: (Y['os_time'].values, Y['os_code'].values)
    elif args.outcome != 'os' or group == 'Progress (P1)' or group == 'Progress (P2)':
        get_target = lambda df: (Y['pfs_time'].values, Y['pfs_code'].values)
    durations_val, events_val = get_target(Y)
    print(durations_val.shape, events_val.shape, preds.shape, Y.shape)

    ev = EvalSurv(preds, durations_val, events_val, censor_surv='km')
    c_ind = ev.concordance_td('antolini')
    return c_ind


def c_harrell(preds, Y, group):
    preds = 1 - preds

    if group == 'Progress (P1)' or group == 'Progress & deceased (P1)':
        preds_at_t = preds.iloc[args.thresh_m_to_p2-1]
    else:
        preds_at_t = preds.iloc[args.preds_time_in_days - 1]

    if args.outcome == 'os' or group == 'Progress & deceased (P1)' or \
                            group == 'Progress & deceased (P2)' or group == 'Non-progress & deceased (P2)':
        Y['os_code'] = Y['os_code'].astype('bool')
        Y = Surv.from_dataframe('os_code', 'os_time', Y)
        c_index = concordance_index_censored(Y['os_code'], Y['os_time'], np.asarray(preds_at_t))

    elif args.outcome != 'os' or group == 'Progress (P1)' or group == 'Progress (P2)':
        Y['pfs_code'] = Y['pfs_code'].astype('bool')
        Y = Surv.from_dataframe('pfs_code', 'pfs_time', Y)
        c_index = concordance_index_censored(Y['pfs_code'], Y['pfs_time'], np.asarray(preds_at_t))

    return c_index[0]


def c_uno(preds, Y, df, group):
    preds = 1 - preds

    if group != 'Progress (P1)' and group != 'Progress & deceased (P1)':
        preds_at_t = preds.iloc[args.preds_time_in_days-1]
    else:
        preds_at_t = preds.iloc[args.thresh_m_to_p2-1]

    if args.outcome == 'os' or group == 'Progress & deceased (P1)' or \
            group == 'Progress & deceased (P2)' or group == 'Non-progress & deceased (P2)':
        cols_test = Y['sample'].to_list()
        cols_train = set(df['sample'].tolist()) - set(cols_test)

        Y_train = df[df['sample'].isin(cols_train)][['os_code', 'os_time']]
        Y_train['os_code'] = df['os_code'].astype('bool')
        Y_train = Surv.from_dataframe('os_code', 'os_time', Y_train)

        Y['os_code'] = Y['os_code'].astype('bool')
        Y = Surv.from_dataframe('os_code', 'os_time', Y)

        try:
            c_index = concordance_index_ipcw(Y_train[['os_code', 'os_time']], Y[['os_code', 'os_time']],
                                         np.asarray(preds_at_t))
        except Exception:
            c_index = [np.nan]

    elif args.outcome != 'os' or group == 'Progress (P1)' or group == 'Progress (P2)':
        cols_test = Y['sample'].to_list()
        cols_train = set(df['sample'].tolist()) - set(cols_test)

        Y_train = df[df['sample'].isin(cols_train)][['pfs_code', 'pfs_time']]
        Y_train['pfs_code'] = df['pfs_code'].astype('bool')
        Y_train = Surv.from_dataframe('pfs_code', 'pfs_time', Y_train)

        Y['pfs_code'] = Y['pfs_code'].astype('bool')
        Y = Surv.from_dataframe('pfs_code', 'pfs_time', Y)
        try:
            c_index = concordance_index_ipcw(Y_train[['pfs_code', 'pfs_time']], Y[['pfs_code', 'pfs_time']],
                                         np.asarray(preds_at_t))
        except Exception:
            c_index = [np.nan]

    return c_index[0]


def brier_score(preds, Y, group):
    if group == 'Progress (P1)' or group == 'Progress & deceased (P1)':
        preds = preds.head(args.thresh_m_to_p2-1)
    else:
        preds = preds.head(args.preds_time_in_days-1)

    if args.outcome == 'os':
        get_target = lambda df: (Y['os_time'].values, Y['os_code'].values)
        time_grid = np.linspace(Y['os_time'].min(), args.preds_time_in_days-1, 100)
    elif group == group == 'Progress (P1)':
        get_target = lambda df: (Y['pfs_time'].values, Y['pfs_code'].values)
        time_grid = np.linspace(Y['pfs_time'].min(), args.thresh_m_to_p2-1, 100)
    elif group == 'Progress & deceased (P1)':
        get_target = lambda df: (Y['os_time'].values, Y['os_code'].values)
        time_grid = np.linspace(Y['os_time'].min(), args.thresh_m_to_p2-1, 100)
    elif args.outcome != 'os' or group == 'Progress (P2)':
        get_target = lambda df: (Y['pfs_time'].values, Y['pfs_code'].values)
        time_grid = np.linspace(Y['pfs_time'].min(), args.preds_time_in_days-1, 100)
    elif group == 'Progress & deceased (P2)' or group == 'Non-progress & deceased (P2)':
        get_target = lambda df: (Y['os_time'].values, Y['os_code'].values)
        time_grid = np.linspace(Y['os_time'].min(), args.preds_time_in_days-1, 100)
    durations_val, events_val = get_target(Y)

    ev = EvalSurv(preds, durations_val, events_val, censor_surv='km')
    int_br_score = ev.integrated_brier_score(time_grid)
    # ev.brier_score(time_grid).plot()
    # plt.ylabel('Brier score')
    return int_br_score


def bll(preds, Y, group):
    if group == 'Progress (P1)' or group == 'Progress & deceased (P1)':
        preds = preds.head(args.thresh_m_to_p2-1)
    else:
        preds = preds.head(args.preds_time_in_days-1)

    if args.outcome == 'os' or group == 'Progress & deceased (P1)' or \
            group == 'Progress & deceased (P2)' or group == 'Non-progress & deceased (P2)':
        get_target = lambda df: (Y['os_time'].values, Y['os_code'].values)
        time_grid = np.linspace(Y['os_time'].min(), Y['os_time'].max(), 100)

    elif args.outcome != 'os' or group == 'Progress (P1)' or group == 'Progress (P2)':
        get_target = lambda df: (Y['pfs_time'].values, Y['pfs_code'].values)
        time_grid = np.linspace(Y['pfs_time'].min(), Y['pfs_time'].max(), 100)

    durations_val, events_val = get_target(Y)

    ev = EvalSurv(preds, durations_val, events_val, censor_surv='km')
    nbll = ev.integrated_nbll(time_grid)
    # ev.nbll(time_grid).plot()
    # plt.ylabel('NBLL')
    return nbll


def auc(preds, Y, df, group):
    if group == 'Progress (P1)' or group == 'Progress & deceased (P1)':
        preds = preds.head(args.thresh_m_to_p2-1)
    else:
        print('nothing')
        # preds = preds.head(args.preds_time_in_days-1)

    if args.outcome == 'os' or group == 'Progress & deceased (P1)' or \
            group == 'Progress & deceased (P2)' or group == 'Non-progress & deceased (P2)':
        Y_test = Surv.from_dataframe('os_code', 'os_time', Y)
        cols_test = list(preds)
        cols_train = set(df['sample'].tolist()) - set(cols_test)
        Y_train = df[df['sample'].isin(cols_train)][['os_code', 'os_time']]
        Y_train = Surv.from_dataframe('os_code', 'os_time', Y_train)
        times = np.percentile(Y_train["os_time"], np.linspace(5, 85, 15))

    elif args.outcome != 'os' or group == 'Progress (P1)' or group == 'Progress (P2)':
        Y_test = Surv.from_dataframe('pfs_code', 'pfs_time', Y)
        cols_test = list(preds)
        cols_train = set(df['sample'].tolist()) - set(cols_test)
        Y_train = df[df['sample'].isin(cols_train)][['pfs_code','pfs_time']]
        # print('Y_train', Y_train)
        # print('df', df)
        Y_train = Surv.from_dataframe('pfs_code', 'pfs_time', Y_train)
        times = np.arange(Y_test['pfs_time'].min(), Y_test['pfs_time'].max(), 1)

    # preds = preds.head(len(times))

    preds = 1 - preds

    try:
        auc_, mean_auc = cumulative_dynamic_auc(Y_train, Y_test, preds.T, times)
    except Exception:
        mean_auc = np.nan
    # plt.plot(times, auc_, marker="o")

    print('mean auc', mean_auc)
    return mean_auc


'''
plot utils
'''
def label_diff(i, j, text, X, Y, ax):
    # Custom function to draw the diff bars
    print(X)

    x = (X[i]+X[j])/2
    y = 1.1*max(Y[i], Y[j])
    dx = abs(X[i]-X[j])

    props = {'connectionstyle':'bar','arrowstyle':'-','shrinkA':20,'shrinkB':20,'linewidth':4}
    ax.annotate(text, xy=(x,y+0.2), zorder=10)
    ax.annotate('', xy=(X[i],y), xytext=(X[j],y), arrowprops=props)


def barplot_annotate_brackets(num1, num2, data, center, he_ind, height, ax, yerr=None, dh=.03, barh=.03, fs=None, maxasterix=None):
    if type(data) is str:
        text = data
    else:
        text = ''
        p = .05

        while data < p:
            text += '*'
            p /= 10.

            if maxasterix and len(text) == maxasterix:
                break

        if len(text) == 0:
            text = 'n. s.'

    lx, ly = center[num1], height[he_ind]
    rx, ry = center[num2], height[he_ind]

    if yerr:
        ly += yerr[num1]
        ry += yerr[num2]

    ax_y0, ax_y1 = plt.gca().get_ylim()
    dh *= (ax_y1 - ax_y0)
    barh *= (ax_y1 - ax_y0)

    y = max(ly, ry) + dh

    barx = [lx, lx, rx, rx]
    bary = [y, y+barh, y+barh, y]
    mid = ((lx+rx)/2, y+barh)

    plt.plot(barx, bary, c='maroon')

    kwargs = dict(ha='center', va='bottom')
    if fs is not None:
        kwargs['fontsize'] = fs

    ax.text(*mid, text, **kwargs, fontsize=8, fontname='Arial')


'''
First level
'''
def rank_and_plot_first_level(df_, group, c_index_type):
    print(group)
    df_ = df_.dropna(subset=[c_index_type])

    df_ = df_[df_['Group'] == group]
    df_sub = df_[[c_index_type, 'Model', 'Kfolds']]

    df_sub = pd.pivot(df_sub, index='Kfolds', columns='Model', values=c_index_type)

    list_of_models = df_sub.median().sort_values(ascending=False).index.to_list()

    df_sub = 1 - df_sub
    df_sub_rank = df_sub.T.rank()

    df_sub_rank = pd.DataFrame(df_sub_rank.T.value_counts().rename('count')).reset_index()

    concat_df = pd.DataFrame()
    for item in list_of_models:
        df_model = df_sub_rank[[item, 'count']]
        df_model = df_model.groupby([item]).agg({'count': 'sum'}).rename(columns={'count': item})
        concat_df = pd.concat([concat_df, df_model], axis=1)

    # rank: perc of folds
    concat_df = concat_df.fillna(0)
    concat_df = concat_df / concat_df.sum()[0] * 100
    # concat_df[['deep_hit_single', 'neural_cox_non_prop', 'cph', 'cph_time', 'rsf']].plot(kind='bar', stacked=True)

    names_dict = {'rsf': 'RSF', 'neural_cox_non_prop': 'Neural CNPH', 'cph': 'CPH-R'}

    names_df = pd.DataFrame(names_dict.items(), columns=['Model_old_name', 'Model_new'])

    names_df = names_df.set_index('Model_old_name')
    concat_df = pd.merge(concat_df.T, names_df[['Model_new']],
                            left_index=True, right_index=True, how='left')

    concat_df = concat_df.set_index('Model_new')
    concat_df[[1, 2, 3]].plot(kind='bar', stacked=True)
    plt.xticks(rotation=45)
    plt.xlabel('')
    plt.legend()
    plt.title(group)
    plt.tight_layout()
    plt.savefig(args.path + '/plots/' + c_index_type + '_rankFolds_' + group + '.png', tight=True)


def rank_folds_first_level():
    """
    rank folds for state models
    """
    output_files = glob.glob(args.path + "/kfold/metrics/model_feat_combo/surv_group_" + str(args.group) + "*.csv")

    df_concat = pd.DataFrame()
    for file_name in output_files:
        df = pd.read_csv(file_name)
        df_concat = pd.concat([df_concat, df], ignore_index=True)

    for group in ['Progress (P1)', 'Progress & deceased (P1)', 'Progress (P2)', 'Progress & deceased (P2)']:
        for c_ind_type in ['C_index_val']:
            rank_and_plot_first_level(df_concat, group, c_ind_type)

    return None


def plot_c_index_state(perf_type):
    """
    State
    """
    output_files = glob.glob(args.path + "/kfold/metrics/model_feat_combo/surv_group_*.csv")
    # print(args.path + "/kfold/metrics/model_feat_combo/surv_group_" + str(args.group) + "*.csv")
    # print(output_files)

    data_list = []
    for filename in output_files:
        single_df = pd.read_csv(filename, index_col=None, header=0)
        single_df['Group_id'] = filename.split('/')[-1].split('_')[2]
        data_list.append(single_df)

    df = pd.concat(data_list, axis=0, ignore_index=True)

    df.mask(df.eq('None')).dropna(inplace=True)

    # filter groups
    df = df[(df['Group'] != 'Move to phase 2') & (df['Group'] != 'Non-progress & deceased (P2)')]

    # edit model name and color
    names_dict = {'rsf': 'RSF', 'cph': 'CPH-R', 'neural_cox_non_prop': 'Neural CNPH'}
    names_df = pd.DataFrame(names_dict.items(), columns=['Model_old_name', 'Model_new'])
    df = pd.merge(df, names_df[['Model_old_name', 'Model_new']],
                            left_on='Model', right_on='Model_old_name', how='left')

    #### ncox only
    df_cp = df.copy()
    df_cp = df_cp[df_cp['Model_new'] == 'Neural CNPH']  # edit / uncomment if more models are needed
    with plt.style.context("ggplot"):
        sn.set_theme(style="whitegrid")
        fig_ncox, ax_ncox = plt.subplots(2,2)
        ax_ncox = ax_ncox.flatten()
        for i, item in enumerate(df['Group'].unique()):
            ax_sub = ax_ncox[i]
            sn.boxplot(x="Model_new", y="C_index" + "_" + perf_type, data=df_cp[df_cp['Group'] == item],
                       palette='ch:s=-.2,r=.6', ax=ax_sub, showfliers=True)
            ax_sub.set(xlabel=None)
            ax_sub.set_xticks([])
            ax_sub.set_ylim([0.48, 0.8])

            legend = ax_sub.legend(prop={'size': 0})
            legend.get_frame().set_alpha(0.9)
            legend.remove()
            ax_sub.set_title(item, fontsize=12, fontname="Arial")
            ax_sub.set_ylabel("C-index", fontsize=12, fontname="Arial")
            ax_sub.tick_params(axis='both', which='major', labelsize=12)
            ax_sub.grid(True)

        plt.tight_layout()
        fig_ncox.savefig(args.path + '/plots/c_index_kfold_model_' + perf_type + '_ncox_only.png')

    #### all models
    df['Model_new'] = df['Model_new'] + ' (Group=' + df['Group_id'] + ')'
    # c index
    with plt.style.context("ggplot"):
        sn.set_theme(style="whitegrid")

        for i, item in enumerate(df['Group'].unique()):
            new_df = pd.DataFrame()
            subset_cdf = df[df['Group'] == item]

            c_df = subset_cdf.groupby(['Group', 'Model_new'])['C_index' + '_' + perf_type].median().sort_values(ascending=False).reset_index()
            median_df = c_df.groupby(['Model_new'])['C_index' + '_' + perf_type].median().sort_values(ascending=False).reset_index()

            for model_item in c_df['Model_new']:
                new_df = new_df.append(subset_cdf[subset_cdf['Model_new'] == model_item])

            fig, ax = plt.subplots()
            sn.boxplot(x="Model_new", y="C_index" + "_" + perf_type, data=new_df, order=list(median_df['Model_new']),
                       palette='ch:s=-.2,r=.6', ax=ax, showfliers=True)

            ax.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
            ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
            ax.grid(b=True, which='major', color='w', linewidth=1.0)
            ax.grid(b=True, which='minor', color='w', linewidth=0.5)

            ax.set(xlabel=None)
            ax.set_ylim([0.43, 1])

            for ax in fig.axes:
                mpl.pyplot.sca(ax)
                plt.xticks(rotation=45)

                ## p-value
                heights = [0.7, 0.75, 0.8, 0.85, 0.92]
                bars = np.arange(6)  # number of models
                he_ind = 0

                # neural cnph with other models
                for a, b in it.combinations(median_df.Model_new.unique().tolist(), 2):
                    if 'Neural' in a or 'Neural' in b:
                        print(a, b)
                        ind_a = median_df.index[median_df['Model_new'] == a].tolist()[0]
                        ind_b = median_df.index[median_df['Model_new'] == b].tolist()[0]

                        stat, p_val = scipy.stats.ranksums(
                            new_df[new_df['Model_new'] == a]['C_index' + '_' + perf_type],
                            new_df[new_df['Model_new'] == b]['C_index' + '_' + perf_type])

                        if p_val < 0.005:
                            text = 'p<0.005'
                        elif p_val < 0.001:
                            text = 'p<0.001'
                        elif p_val < 0.01:
                            text = 'p<0.01'
                        elif p_val < 0.05:
                            text = 'p<0.05'
                        elif p_val > 0.05:
                            text = 'p=' + str(round(p_val, 4))

                        print('ok', ind_a, ind_b, text)
                        barplot_annotate_brackets(ind_a, ind_b, text, bars, he_ind, heights, ax)
                        he_ind += 1

            legend = ax.legend(prop={'size': 0})
            legend.get_frame().set_alpha(0.9)
            legend.remove()
            ax.set_title(item, fontsize=12, fontname="Arial")
            ax.set_ylabel("C-index" + '(' + perf_type + ')', fontsize=12, fontname="Arial")
            ax.tick_params(axis='both', which='major', labelsize=12)
            ax.grid(True)

            plt.tight_layout()
            fig.savefig(args.path + '/plots/c_index_kfold_model_' + perf_type
                        + '_' + item + '.png')


'''
Multistate
'''
def rank_and_plot_multistate(df_, outcom, group, c_index_type):
    """
    rank and plot folds for multistate only
    """
    print(outcom, group)
    fig, ax = plt.subplots()
    sn.set(style='whitegrid')
    df_copy = df_.copy()

    if group == '2~10~100' and not args.ext_val:
        df_ = df_copy[df_copy['X_ID'].isin([
            'M-Neural CNPH (Group=100)', 'M-RSF (Group=100)', 'M-CPH-R (Group=100)',
            'M-CPH-R (Group=10)', 'M-CPH-R (Group=2)'])
        ]
    elif group == '0~1~2~3~4~5~6~7~9~10~100' and not args.ext_val:
        df_ = df_copy[df_copy['X_ID'].isin([
            'M-CPH-R (Group=0)', 'M-CPH-R (Group=1)', 'M-CPH-R (Group=2)',
            'M-Neural CNPH (Group=3)', 'M-Neural CNPH (Group=4)',
            'M-Neural CNPH (Group=5)', 'M-Neural CNPH (Group=6)', 'M-Neural CNPH (Group=7)',
            'M-Neural CNPH (Group=9)', 'M-CPH-R (Group=10)', 'M-Neural CNPH (Group=100)'])
        ]

    elif group == '2~9' and args.ext_val:
        df_ = df_copy[df_copy['X_ID'].isin([
            'M-Neural CNPH (Group=9)', 'M-CPH-R (Group=2)'])
        ]

    elif group == '2~9~100' and args.ext_val:
        df_ = df_copy[df_copy['X_ID'].isin([
            'M-Neural CNPH (Group=100)', 'M-RSF (Group=100)', 'M-CPH-R (Group=100)',
            'M-Neural CNPH (Group=9)', 'M-CPH-R (Group=2)'])
        ]

    if group == '0~1~2~3~4~5~6~7~9~10~100' and not args.ext_val:
        order_ = ['M-CPH-R (Group=10)', 'M-CPH-R (Group=2)', 'M-CPH-R (Group=0)', 'M-CPH-R (Group=1)',
                  'M-Neural CNPH (Group=3)', 'M-Neural CNPH (Group=4)',
                  'M-Neural CNPH (Group=5)', 'M-Neural CNPH (Group=6)', 'M-Neural CNPH (Group=7)',
                  'M-Neural CNPH (Group=9)', 'M-Neural CNPH (Group=100)']
    elif group == '2~10~100' and not args.ext_val:
        order_ = ['M-CPH-R (Group=10)', 'M-CPH-R (Group=2)', 'M-CPH-R (Group=100)',
                  'M-RSF (Group=100)', 'M-Neural CNPH (Group=100)']
    elif group == '2~9' and args.ext_val:
        order_ = ['M-CPH-R (Group=2)', 'M-Neural CNPH (Group=9)']
    elif group == '2~9~100' and args.ext_val:
        order_ = ['M-CPH-R (Group=2)', 'M-CPH-R (Group=100)', 'M-RSF (Group=100)',
                  'M-Neural CNPH (Group=100)', 'M-Neural CNPH (Group=9)']

    rename_x_id = {
                    'M-CPH-R (Group=100)': 'M-CPH-R', 'M-RSF (Group=100)': 'M-RSF',
                    'M-Neural CNPH (Group=100)': 'All features',
                    'M-CPH-R (Group=0)': 'R-ISS', 'M-CPH-R (Group=1)': 'R2-ISS',
                    'M-CPH-R (Group=2)': 'ISS',
                    'M-Neural CNPH (Group=3)': 'Plus_CLIN',
                    'M-Neural CNPH (Group=4)': 'Plus_Rec. Transloc',
                    'M-Neural CNPH (Group=5)': 'Plus_SCT',
                    'M-Neural CNPH (Group=6)': 'Plus_Cont. treat',
                    'M-Neural CNPH (Group=7)': 'Plus_Therapy',
                    'M-Neural CNPH (Group=8)': 'p<0.01 (univ) only',
                    'M-Neural CNPH (Group=9)': 'Plus_Top genomics',
                    'M-CPH-R (Group=10)': 'GEP70',
                    'M-Neural CNPH (Group=11)': 'ISS+SCT',
                    'M-Neural CNPH (Group=12)': 'ISS+Cont. treat',
                    'M-Neural CNPH (Group=13)': 'ISS+Rec. Transloc',
                    'M-Neural CNPH (Group=14)': 'ISS+Therapy',
                    'M-Neural CNPH (Group=15)': 'ISS+TopGen'
                   }

    rename_x_id_df = pd.DataFrame.from_dict(rename_x_id, orient='index').reset_index().rename(columns={'index': 'old', 0: 'new_X_ID'})

    df_ = pd.merge(df_, rename_x_id_df, left_on='X_ID', right_on='old')

    df_ = df_.dropna(subset=[c_index_type])

    df_sub = df_[[c_index_type, 'Fold', 'new_X_ID', 'X_ID']]
    df_sub_copy = df_sub.copy()

    df_[[c_index_type, 'Fold', 'new_X_ID']].to_csv(args.path + '/plots/' + 'cleanLabels_' + str(args.outcome) + '_'
                                        + str(args.get_metric) + '_ext_val' + str(args.ext_val) + '_' + str(args.preds_time_in_days) + '_group_' + group + '.csv', index=False)

    df_sub = pd.pivot(df_sub, index='Fold', columns='new_X_ID', values=c_index_type)

    list_of_models = df_sub_copy.groupby(['X_ID', 'new_X_ID'])[args.get_metric].median().to_frame().reset_index().set_index('X_ID').reindex(order_)['new_X_ID'].to_list()

    df_sub_c_ind_vals = df_sub.copy()

    df_sub = 1 - df_sub
    df_sub_rank = df_sub.T.rank(method='min')

    df_sub_rank = pd.DataFrame(df_sub_rank.T.value_counts().rename('count')).reset_index()

    concat_df = pd.DataFrame() # percentage ranked order
    for item in list_of_models:
        df_model = df_sub_rank[[item, 'count']]
        df_model = df_model.groupby([item]).agg({'count': 'sum'}).rename(columns={'count': item})
        concat_df = pd.concat([concat_df, df_model], axis=1)

    concat_df = concat_df.fillna(0)
    concat_df = concat_df / concat_df.sum()[0] #* 100

    concat_df = concat_df.T
    concat_df[list(range(1,concat_df.shape[0]+1))].plot(kind='bar',
                                                        width=1, ax=ax, stacked=True,
                                                        colormap='tab20')

    #ax.yaxis.set_major_locator(MultipleLocator(2))
    #ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.tick_params(which='both', width=2)
    ax.tick_params(which='major', length=7)
    ax.tick_params(which='minor', length=4)

    ax.grid(visible=True, which='major', axis='x', color='w', linewidth=1.0)
    ax.grid(visible=True, which='minor', axis='x', color='w', linewidth=0.5)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_color('black')
    ax.spines['left'].set_linewidth(1)

    ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fancybox=True, shadow=True, prop={'family': 'Arial', 'size': 10})

    # plt.legend()
    # plt.legend().get_frame().set_alpha(0.35)
    plt.xticks(rotation=80, fontsize=10, fontname='Arial')
    plt.xlabel('')
    ax.set_ylabel('Fraction of folds', fontsize=14, fontname="Arial")
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(args.path + '/plots/' + 'rankFolds_' + c_index_type + '_ext_val_' + str(args.ext_val) + '_' +
                            outcom + '_group_' + group + '.png', tight=True)

    # concat_df.to_csv(args.path + '/plots/rankfolds_' + group + '_' + str(args.groups_0to8) + '.txt', sep='\t')

    return df_sub_c_ind_vals


def rank_folds_multistate():
    """
    rank folds for multistate models
    """
    df_sub_c_ind_vals = None

    if args.ext_val:
        group_list = ['2~9', '2~9~100']
    else:
        group_list = ['2~10~100', '0~1~2~3~4~5~6~7~9~10~100']

    if args.plot_sub_groups:
        for group in group_list:
            for outcom in [args.outcome]:
                if args.ext_val is not True:
                    df = pd.read_csv(args.path + '/plots/' + args.outcome + '_ext_val_' + str(args.ext_val) + '_'
                        + str(args.preds_time_in_days) + '_group_' + '0~1~2~3~4~5~6~7~8~9~10~11~12~13~14~15~100' + '.csv')
                else:
                    df = pd.read_csv(args.path + '/plots/' + args.outcome + '_ext_val_' + str(args.ext_val) + '_'
                                     + str(args.preds_time_in_days) + '_group_' + '2~9~100' + '.csv')

                for c_ind_type in [args.get_metric]: # , 'C_index_harrell', 'C_index_uno'
                    df_sub_c_ind_vals = rank_and_plot_multistate(df, outcom, group, c_ind_type)

    return df_sub_c_ind_vals


def plot_metrics_multistate_models(df_metr_conc):
    """
    Multistate
    """
    print(df_metr_conc)
    if args.c_index_group_addition:
        item = 'X_ID'
        # string_vals = args.group_feat_add.split('~')
        # int_vals = [int(i) for i in string_vals]
        save_fil = 'group_feat_add'

    df_metr_conc.drop_duplicates(subset=['Model', 'Fold', 'Group_id'], inplace=True)

    names_dict = {'rsf': 'M-RSF', 'neural_cox_prop': 'M-Neural CPH', 'neural_cox_non_prop': 'M-Neural CNPH',
                  'cph': 'M-CPH-R', 'cph_ISS': 'M-CPH-R (ISS)',
                  'cph_R_ISS': 'M-CPH-R (R-ISS)', 'cph_R2_ISS': 'M-CPH-R (R2-ISS)'}

    names_df = pd.DataFrame(names_dict.items(), columns=['Model_old_name', 'Model_new'])

    df_metr_conc = pd.merge(df_metr_conc, names_df[['Model_old_name', 'Model_new']],
                            left_on='Model', right_on='Model_old_name', how='left')

    df_metr_conc['X_ID'] = df_metr_conc['Model_new'] + ' (Group=' + df_metr_conc['Group_id'] + ')'

    df_metr_conc.to_csv(args.path + '/plots/' + args.outcome + '_ext_val_' + str(args.ext_val) + '_' + str(args.preds_time_in_days) + '_group_' + str(args.group_feat_add) + '.csv', index=False) # comment

    if args.plot_sub_groups:
        # comment
        # del df_metr_conc
        # df_metr_conc = pd.read_csv(args.path + '/plots/' + args.outcome + '_ext_val_' + str(args.ext_val) + '_' + str(args.preds_time_in_days) + '_group_' + '0~1~2~3~4~5~6~7~8~9~10~11~12~13~14~15~100' + '.csv') # comment

        if args.ext_val:
            group_list = ['2~9', '2~9~100']
        else:
            group_list = ['2~10~100', '0~1~2~3~4~5~6~7~9~10~100']

        df_metr_conc_cp = df_metr_conc.copy()
        for group in group_list:
            if group == '2~10~100' and not args.ext_val:
                comp_ncnph = 'M-Neural CNPH (Group=100)' # feature group to compare for p-val stats
                comp_rsf = 'M-RSF (Group=100)' ; comp_cph = 'M-CPH-R (Group=100)'
                comp_gep = 'M-CPH-R (Group=10)' ; comp_iss = 'M-CPH-R (Group=2)'

                # GEP, ISS, All features
                df_metr_conc = df_metr_conc_cp[df_metr_conc_cp['X_ID'].isin([
                                    'M-Neural CNPH (Group=100)', 'M-RSF (Group=100)', 'M-CPH-R (Group=100)',
                                    'M-CPH-R (Group=10)', 'M-CPH-R (Group=2)'])
                                        ]
                if args.outcome == 'os':
                    heights = [0.675, 0.7, 0.745, 0.77]
                elif args.outcome == 'efs':
                    heights = [0.62, 0.7, 0.72, 0.75]

            elif group == '0~1~2~3~4~5~6~7~9~10~100' and not args.ext_val:
                comp_ncnph_100 = 'M-Neural CNPH (Group=100)' # feature group to compare for p-val stats
                comp_ncnph_9 = 'M-Neural CNPH (Group=9)' ; comp_ncnph_7 = 'M-Neural CNPH (Group=7)' ;
                comp_ncnph_6 = 'M-Neural CNPH (Group=6)' ; comp_ncnph_5 = 'M-Neural CNPH (Group=5)' ;
                comp_ncnph_4 = 'M-Neural CNPH (Group=4)' ; comp_ncnph_3 = 'M-Neural CNPH (Group=3)' ;
                comp_riss = 'M-CPH-R (Group=0)' ; comp_r2iss = 'M-CPH-R (Group=1)'
                comp_gep = 'M-CPH-R (Group=10)' ; comp_iss = 'M-CPH-R (Group=2)'

                df_metr_conc = df_metr_conc_cp[df_metr_conc_cp['X_ID'].isin([
                                'M-CPH-R (Group=0)', 'M-CPH-R (Group=1)', 'M-CPH-R (Group=2)',
                                'M-Neural CNPH (Group=3)', 'M-Neural CNPH (Group=4)',
                                'M-Neural CNPH (Group=5)', 'M-Neural CNPH (Group=6)', 'M-Neural CNPH (Group=7)',
                                'M-Neural CNPH (Group=9)', 'M-CPH-R (Group=10)', 'M-Neural CNPH (Group=100)'])
                                ]
                if args.outcome == 'os':
                    heights = [0.68, 0.65, 0.68, 0.7, 0.73, 0.76, 0.78, 0.76, 0.78, 0.76]
                elif args.outcome == 'efs':
                    heights = [0.62, 0.65, 0.61, 0.625, 0.65, 0.7, 0.725, 0.74, 0.77, 0.74]

            elif group == '2~9' and args.ext_val:
                comp_ncnph = 'M-Neural CNPH (Group=9)' # feature group to compare for p-val stats
                comp_iss = 'M-CPH-R (Group=2)'

                # GEP, ISS, All features
                df_metr_conc = df_metr_conc_cp[df_metr_conc_cp['X_ID'].isin([
                                    'M-Neural CNPH (Group=9)', 'M-CPH-R (Group=2)'])
                                        ]
                if args.outcome == 'os':
                    heights = [0.675]
                elif args.outcome == 'efs':
                    heights = [0.62]

            elif group == '2~9~100' and args.ext_val:
                comp_ncnph = 'M-Neural CNPH (Group=100)' # feature group to compare for p-val stats
                comp_rsf = 'M-RSF (Group=100)' ; comp_cph = 'M-CPH-R (Group=100)'
                comp_9 = 'M-Neural CNPH (Group=9)' ; comp_iss = 'M-CPH-R (Group=2)'

                # ISS, All features
                df_metr_conc = df_metr_conc_cp[df_metr_conc_cp['X_ID'].isin([
                                    'M-Neural CNPH (Group=100)', 'M-RSF (Group=100)', 'M-CPH-R (Group=100)',
                                    'M-Neural CNPH (Group=9)', 'M-CPH-R (Group=2)'])
                                        ]
                if args.outcome == 'os':
                    heights = [0.675, 0.7, 0.745, 0.77, 0.725]
                elif args.outcome == 'efs':
                    heights = [0.62, 0.64, 0.66, 0.68, 0.63]

            for metric in [args.get_metric]:
                fig, ax = plt.subplots()
                sn.set(style='whitegrid')

                if metric == 'Int_Brier_Score' or metric == 'NBLL':
                    median_df = df_metr_conc.groupby([item])[metric].median().sort_values(ascending=True).reset_index()

                if args.c_index_group_addition:
                    median_df = df_metr_conc.groupby([item])[metric].median().sort_values(ascending=False).reset_index()

                if group == '0~1~2~3~4~5~6~7~9~10~100' and not args.ext_val:
                    order_ = ['M-CPH-R (Group=10)', 'M-CPH-R (Group=2)', 'M-CPH-R (Group=0)', 'M-CPH-R (Group=1)',
                                'M-Neural CNPH (Group=3)', 'M-Neural CNPH (Group=4)',
                                'M-Neural CNPH (Group=5)', 'M-Neural CNPH (Group=6)', 'M-Neural CNPH (Group=7)',
                                'M-Neural CNPH (Group=9)', 'M-Neural CNPH (Group=100)']
                elif group == '2~10~100' and not args.ext_val:
                    order_ = ['M-CPH-R (Group=10)', 'M-CPH-R (Group=2)', 'M-CPH-R (Group=100)',
                              'M-RSF (Group=100)', 'M-Neural CNPH (Group=100)']
                elif group == '2~9' and args.ext_val:
                    order_ = ['M-CPH-R (Group=2)', 'M-Neural CNPH (Group=9)']
                    # order_ = list(median_df[item])
                elif group == '2~9~100' and args.ext_val:
                    order_ = ['M-CPH-R (Group=2)', 'M-CPH-R (Group=100)', 'M-RSF (Group=100)',
                              'M-Neural CNPH (Group=100)', 'M-Neural CNPH (Group=9)']

                median_df = median_df.set_index('X_ID').reindex(order_).reset_index()
                bars = np.arange(median_df.shape[0]) # max. number of models

                sn.boxplot(x=item, y=metric, data=df_metr_conc, order=order_,
                           palette='YlOrBr', showfliers=False, showmeans=False, ax=ax)

                ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
                ax.grid(visible=True, which='major', axis='x', color='w', linewidth=1.0)
                ax.grid(visible=True, which='minor', axis='x', color='w', linewidth=0.5)
                ax.set(xlabel=None)

                if group == '2~10~100' and not args.ext_val:
                    if args.outcome == 'os':
                        ax.set_ylim([0.45, 0.8])
                    elif args.outcome == 'efs':
                        ax.set_ylim([0.45, 0.8])
                elif group == '0~1~2~3~4~5~6~7~9~10~100' and not args.ext_val:
                    if args.outcome == 'os':
                        ax.set_ylim([0.45, 0.84])
                    elif args.outcome == 'efs':
                        ax.set_ylim([0.45, 0.8])
                elif group == '2~9' and args.ext_val:
                    if args.outcome == 'os':
                        ax.set_ylim([0.45, 0.8])
                    elif args.outcome == 'efs':
                        ax.set_ylim([0.45, 0.8])
                if group == '2~9~100' and args.ext_val:
                    if args.outcome == 'os':
                        ax.set_ylim([0.45, 0.8])
                    elif args.outcome == 'efs':
                        ax.set_ylim([0.45, 0.72])

                if 'C_index_' in metric:
                    metric_label = 'C-index'
                else:
                    metric_label = metric

                ax.set_ylabel(metric_label, fontsize=14, fontname="Arial")
                plt.xticks(fontsize=8, fontname="Arial")
                plt.yticks(fontsize=10, fontname="Arial")

                for ax in fig.axes:
                    mpl.pyplot.sca(ax)
                    plt.xticks(rotation=45)
                    # ax.axes.get_xaxis().set_visible(False)
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)
                    ax.spines['bottom'].set_visible(False)
                    ax.spines['left'].set_color('black')
                    ax.spines['left'].set_linewidth(1)

                he_ind = 0

                '''Model comparisons for significance'''
                for a, b in it.combinations(median_df[item].unique().tolist(), 2):
                    print('before', a, b)

                    if group == '2~10~100' and not args.ext_val:
                        if (comp_ncnph in a and comp_rsf in b) or (comp_rsf in a and comp_ncnph in b) or \
                            (comp_ncnph in a and comp_cph in b) or (comp_cph in a and comp_ncnph in b) or \
                                (comp_cph in a and comp_iss in b) or (comp_iss in a and comp_cph in b) \
                                    or (comp_gep in a and comp_iss in b) or (comp_iss in a and comp_gep in b):
                            '''pairwise rank-sums test'''
                            calc_pval(median_df, a, b, df_metr_conc, item, metric, bars, he_ind, heights, ax)
                            he_ind +=1
                        else:
                            print('else', a, b)
                    elif group == '0~1~2~3~4~5~6~7~9~10~100' and not args.ext_val:
                        if (comp_gep in a and comp_riss in b) or (comp_riss in a and comp_gep in b) \
                            or (comp_riss in a and comp_iss in b) or (comp_iss in a and comp_riss in b) \
                            or (comp_r2iss in a and comp_iss in b) or (comp_iss in a and comp_r2iss in b) \
                            or (comp_iss in a and comp_ncnph_3 in b) or (comp_ncnph_3 in a and comp_iss in b) \
                            or (comp_ncnph_4 in a and comp_ncnph_3 in b) or (comp_ncnph_3 in a and comp_ncnph_4 in b) \
                            or (comp_ncnph_5 in a and comp_ncnph_4 in b) or (comp_ncnph_4 in a and comp_ncnph_5 in b) \
                            or (comp_ncnph_6 in a and comp_ncnph_5 in b) or (comp_ncnph_5 in a and comp_ncnph_6 in b) \
                            or (comp_ncnph_7 in a and comp_ncnph_6 in b) or (comp_ncnph_6 in a and comp_ncnph_7 in b) \
                            or (comp_ncnph_9 in a and comp_ncnph_7 in b) or (comp_ncnph_7 in a and comp_ncnph_9 in b) \
                                or (comp_ncnph_9 in a and comp_ncnph_100 in b) or (comp_ncnph_100 in a and comp_ncnph_9 in b):

                            '''pairwise rank-sums test'''
                            calc_pval(median_df, a, b, df_metr_conc, item, metric, bars, he_ind, heights, ax)
                            he_ind +=1
                        else:
                            print('else', a, b)
                    elif group == '2~9' and args.ext_val:
                        if (comp_ncnph in a and comp_iss in b) or (comp_iss in a and comp_ncnph in b):
                            '''pairwise rank-sums test'''
                            calc_pval(median_df, a, b, df_metr_conc, item, metric, bars, he_ind, heights, ax)
                            he_ind +=1
                        else:
                            print('else', a, b)
                    if group == '2~9~100' and args.ext_val:
                        if (comp_9 in a and comp_rsf in b) or (comp_rsf in a and comp_9 in b) or \
                            (comp_ncnph in a and comp_rsf in b) or (comp_rsf in a and comp_ncnph in b) or \
                                (comp_9 in a and comp_cph in b) or (comp_cph in a and comp_9 in b) or \
                                 (comp_cph in a and comp_iss in b) or (comp_iss in a and comp_cph in b) \
                                    or (comp_9 in a and comp_ncnph in b) or (comp_ncnph in a and comp_9 in b):
                            '''pairwise rank-sums test'''
                            calc_pval(median_df, a, b, df_metr_conc, item, metric, bars, he_ind, heights, ax)
                            he_ind +=1
                        else:
                            print('else', a, b)
                plt.tight_layout()
                fig.savefig(args.path + '/plots/' + metric + '_' + args.outcome + '_ext_val_' + str(args.ext_val) + '_group_' + group + '.png')

    return None


def calc_pval(median_df, a, b, df_metr_conc, item, metric, bars, he_ind, heights, ax):
    ind_a = median_df.index[median_df[item] == a].tolist()[0]
    ind_b = median_df.index[median_df[item] == b].tolist()[0]

    stat, p_val = scipy.stats.ranksums(df_metr_conc[df_metr_conc[item] == a][metric], df_metr_conc[df_metr_conc[item] == b][metric])

    if p_val < 0.005:
        text = 'p<0.005'
    elif p_val < 0.001:
        text = 'p<0.001'
    elif p_val < 0.01:
        text = 'p<0.01'
    elif p_val < 0.05:
        text = 'p<0.05'
    elif p_val > 0.05:
        text = 'p=' + str(round(p_val, 2))

    barplot_annotate_brackets(ind_a, ind_b, text, bars, he_ind, heights, ax)
    return ax


def multistate(df):
    if args.c_index_group_addition and args.outcome != 'None':
        output_files = []
        for group in args.group_feat_add.split('~'):
            if args.ext_val:
                output_files.extend(
                    glob.glob(args.path + '/kfold/state_probs/' + args.outcome + '/testgroup-' + group + '-*'))
            else:
                output_files.extend(glob.glob(args.path + '/kfold/state_probs/' + args.outcome + '/group-' + group + '-*'))

    df_metr_conc = pd.DataFrame()
    for filename in output_files:
        df_metr = prep_ground_truth(filename, df, df_metr_conc, 'multistate')
        df_metr_conc = pd.concat([df_metr_conc, df_metr], ignore_index=True)

    return df_metr_conc


'''
Genomic feature selection (model)
'''
def wilcox_feat_signi(group, df_group, df_feat_imp):
    # for p-value;
    permut_files = glob.glob(args.path + '/permut_feat_imps/compose/fold_*_all_feat_weights_' + group[0] + '_*.csv')
    df_concat_ = pd.DataFrame()
    for filename in permut_files:
        single_df = pd.read_csv(filename, index_col=None, header=0)
        single_df['fold'] = filename.split('/')[-1].split('_')[1]
        df_concat_ = df_concat_.append(single_df, ignore_index=True)

    df_group['pval'], df_group['signfic_pval'], df_group['text_signf'] = '', '', ''

    print(list(df_feat_imp))

    for ind, feat in enumerate(df_group[group[0]].tolist()):
        print(feat, df_feat_imp[df_feat_imp['display_id'] == feat].size)
        if df_feat_imp[df_feat_imp['display_id'] == feat].size != 0:
            try:
                w, p_val = wilcoxon(df_feat_imp[df_feat_imp['display_id'] == feat]['weight'].tolist())
                print(w, p_val)

                if p_val < 0.005:
                    signif = 2 ; text = '****'
                elif p_val < 0.001:
                    signif = 2 ; text = '***'
                elif p_val < 0.01:
                    signif = 2 ; text = '**'
                elif p_val < 0.05:
                    signif = 2; text = '*'
                elif p_val > 0.05:
                    signif = 0 ; text = ''

                df_group.loc[ind, 'pval'] = p_val
                df_group.loc[ind, 'signfic_pval'] = signif
                df_group.loc[ind, 'text_signf'] = text
            except Exception:
                print('exception')
    return df_group


def plot_c_index_model_bygroups(perf_type):
    # from baseline = clinical + treatment + ISS + SCT + cont. treat
    # add on genomics
    output_files = glob.glob(args.path + "/kfold/metrics/model_genomics/*_" + args.surv_model + "*.csv")
    output_files.sort()

    data_list = []
    for filename in output_files:
        single_df = pd.read_csv(filename, index_col=None, header=0)
        feats_group = int(filename.split('/')[-1].split('_')[2])

        gp_name = str(feats_group)

        single_df['Feat_group'] = gp_name
        data_list.append(single_df)

    df = pd.concat(data_list, axis=0, ignore_index=True)

    df.mask(df.eq('None')).dropna(inplace=True)

    top_gen_feats_by_group = {}
    # c index
    with plt.style.context("ggplot"):
        sn.set_theme(style="whitegrid")

        fig, ax = plt.subplots(2, 2, figsize=(26, 26), sharex=True)

        for i, item in enumerate(['Progress (P1)', 'Progress (P2)',
                                  'Progress & deceased (P1)', 'Progress & deceased (P2)',
                                  ]):

            if 0 <= i <= 1:
                j = 0
            elif 2 <= i <= 3:
                j = 1 ; i = i - 2
            subset_cdf = df[df['Group'] == item]
            subset_cdf['Feat_group'] = subset_cdf['Feat_group'].astype(int)

            c_df = subset_cdf.groupby(['Group', 'Feat_group'])['C_index' + '_' + perf_type].median().sort_values(ascending=False).reset_index()

            # get perc drop from model (with all features)
            c_df_copy = c_df.copy()
            baseline_c = c_df_copy[c_df_copy['Feat_group'] == 0]['C_index_' + perf_type].values.item()
            c_df_copy['actual_change'] = (c_df_copy['C_index_' + perf_type] - baseline_c)/baseline_c *100
            # c_df_copy.drop('C_index_val', axis=1, inplace=True)

            # statistically significant difference in cross validated c-index performance over baseline
            dict_p_vals = {}
            for a, b in map(lambda e: (e, 0), range(0,10)):
                try:
                    stat, p_val = scipy.stats.ranksums(subset_cdf[subset_cdf['Feat_group'] == a]['C_index_val'],
                                                   subset_cdf[subset_cdf['Feat_group'] == b]['C_index_val'])
                    dict_p_vals[int(a)] = p_val
                except Exception:
                    l = 'none'

            # choose only when improvement over baseline
            filter_keys = c_df_copy[c_df_copy['actual_change'] > 0]
            filter_keys['p-val'] = filter_keys.index.map(dict_p_vals)
            filter_keys = filter_keys[filter_keys['Feat_group'] <= 5]
            filter_keys.dropna(inplace=True)
            dict_p_vals = {your_key: dict_p_vals[your_key]
                           for your_key in filter_keys['Feat_group'].tolist()}

            # most significant set of genomic features from baselines
            ind_min_p_val = min(dict_p_vals, key=dict_p_vals.get)
            act_change_best_feat = c_df_copy[c_df_copy['Feat_group'] == ind_min_p_val]
            print('after filter', item, baseline_c, ind_min_p_val, dict_p_vals[ind_min_p_val], act_change_best_feat['actual_change'].to_numpy()[0])

            top_gen_feats_by_group[item] = ind_min_p_val # save

            ###### change in c-index from baseline with number of features
            subset_cdf = c_df_copy[c_df_copy['Group'] == item]
            c_df = subset_cdf.sort_values(by=['Feat_group'], ascending=True).reset_index(drop=True)
            c_df.plot(x="Feat_group", y=["actual_change"], ax=ax[j, i], linewidth=4)
            ax[j,i].tick_params(axis='both', which='major', labelsize=32)
            ax[j,i].set(xlabel=None)
            ax[j,i].axhline(0, color='black', linestyle="-", linewidth=3)
            legend = ax[j, i].legend(prop={'size': 19})
            legend.get_frame().set_alpha(0.1)
            legend.remove()
            ax[j, i].set_xlabel("# Features (ranked by model)", fontsize=36, fontname="Arial")
            ax[j, i].set_ylabel("Change in c-index from baseline (%)", fontsize=36, fontname="Arial")
            ax[j, i].set_title(item, fontsize=36, fontname="Arial")

            if p_val < 0.005:
                text = 'p<0.005'
            elif p_val < 0.001:
                text = 'p<0.001'
            elif p_val < 0.01:
                text = 'p<0.01'
            elif p_val < 0.05:
                text = 'p<0.05'
            elif p_val > 0.05:
                text = 'p=' + str(round(p_val, 4))

            if item == 'Progress (P1)':
                y_scal = 2
                x_ = 2; y_ = 0.05; x_text = 2; y_text = 0.05
            elif item == 'Progress & deceased (P1)':
                y_scal = 1.5
                x_ = 2; y_ = 0.05; x_text = 2; y_text = 0.05
            elif item == 'Progress (P2)':
                y_scal = 1
                x_ = 2; y_ = 0.05; x_text = 2; y_text = 0.05
            elif item == 'Progress & deceased (P2)':
                y_scal = 3.5
                x_ = 2; y_ = 0.05; x_text = 2; y_text = 0.05
            elif item == 'Non-progress & deceased (P2)':
                y_scal = 6
                x_ = 2; y_ = 0.05; x_text = 2; y_text = 0.05
            elif item == 'Move to phase 2':
                y_scal = 1.5
                x_ = 2; y_ = 0.05; x_text = 2; y_text = 0.05

            ax[j,i].annotate('d=' + str(ind_min_p_val) + '; c-index=' + str(format(act_change_best_feat['C_index_' + perf_type].to_numpy()[0], ".3f")) + ';\n' + text,
                        xy = (ind_min_p_val, act_change_best_feat['actual_change'].to_numpy()[0]), xytext = (ind_min_p_val + 2, y_scal),
                             horizontalalignment='left',
                             verticalalignment='bottom', arrowprops={"arrowstyle":"->", "color":"red", "lw": 4},
                             fontsize=28, fontname="Arial")

            ax[j, i].annotate('c-index=' + str(format(baseline_c, ".3f")),
                              xy=(x_, y_),
                              xytext=(x_text, y_text),
                              horizontalalignment='left',
                              verticalalignment='bottom',
                              fontsize=28, fontname="Arial")

        plt.suptitle(perf_type.upper(), fontsize=38, fontname="Arial")
        plt.tight_layout()
        fig.savefig(args.path + '/plots/change_in_c_index_state_' + perf_type + '.png')
        plt.close()

        return top_gen_feats_by_group


def heatmap_bygroups_imp_feats(top_n_gen, perf_type):
    union_feats = []
    df_ = pd.DataFrame()
    df_pval_concat = pd.DataFrame()

    # get unique list of features from all groups as provided by Neural cox-nph
    for gp in [('Progress (P1)', top_n_gen['Progress (P1)']), ('Progress & deceased (P1)', top_n_gen['Progress & deceased (P1)']),
               ('Progress (P2)', top_n_gen['Progress (P2)']), ('Progress & deceased (P2)', top_n_gen['Progress & deceased (P2)'])]:

        if args.re_rank_feats:
            df_feat_rank = pd.read_csv(
                args.path + '/permut_feat_imps/aggregate_reranked/' + gp[0] + '_neural_cox_non_prop_folds_feat_imp.csv')
            df_feat_rank = df_feat_rank[['column_id', 'color', 'weight', 'std']]
        else:
            df_feat_rank = pd.read_csv(args.path + '/permut_feat_imps/aggregate/' + gp[0] + '_neural_cox_non_prop_folds_feat_imp.csv')

        overall_fea_class = pd.read_excel(args.legend)
        df_feat_rank = pd.merge(df_feat_rank, overall_fea_class[['column_id', 'display_id', 'class_generic']], on='column_id', how='inner')

        # top n genomic features
        top_gen_feats = df_feat_rank[df_feat_rank['class_generic'] == 'Genomics'].iloc[:gp[1]] # genomics
        # top_gen_feats = df_feat_rank # any class of feature

        union_feats.extend(top_gen_feats['display_id'].tolist())

        df_group = pd.DataFrame({gp[0]: top_gen_feats['display_id'].tolist()})

        df_ = pd.concat([df_, df_group], axis=1)

        df_pval = wilcox_feat_signi(gp, df_group, df_feat_rank)

        # to verify
        df_pval_concat = pd.concat([df_pval_concat, df_pval[[gp[0], 'pval', 'signfic_pval', 'text_signf']]], axis=0)

    unique_elements_union = list(set(union_feats))

    df_heatmap = pd.DataFrame(index=unique_elements_union,
                        columns=['Progress (P1)', 'Progress & deceased (P1)', 'Progress (P2)', 'Progress & deceased (P2)'])

    for feat in unique_elements_union:
        for idx, row in df_.iterrows():
            if row['Progress (P1)'] == feat:
                df_heatmap.loc[feat, 'Progress (P1)'] = 1
            if row['Progress & deceased (P1)'] == feat:
                df_heatmap.loc[feat, 'Progress & deceased (P1)'] = 1
            if row['Progress (P2)'] == feat:
                df_heatmap.loc[feat, 'Progress (P2)'] = 1
            if row['Progress & deceased (P2)'] == feat:
                df_heatmap.loc[feat, 'Progress & deceased (P2)'] = 1

    df_heatmap.fillna(0, inplace=True)
    df_heatmap.sort_index(axis = 0, inplace=True)
    df_heatmap_pval = df_heatmap.copy()
    df_heatmap_pval_symbol = df_heatmap_pval.copy()

    for item in ['Progress (P1)', 'Progress & deceased (P1)',
                 'Progress (P2)', 'Progress & deceased (P2)']:
        subset_ = df_pval_concat.dropna(subset=[item])[[item, 'signfic_pval', 'text_signf']]
        for idx, row in subset_.iterrows():
            if row['signfic_pval'] == 2: # only write if significant
                df_heatmap_pval.loc[row[item], item] = row['signfic_pval']
            df_heatmap_pval_symbol.loc[row[item], item] = row['text_signf']

    # to francesco
    #df_heatmap_pval_symbol = df_heatmap_pval_symbol.replace('', 0)
    #df_heatmap_pval_symbol = df_heatmap_pval_symbol.replace('*', 1)
    #df_heatmap_pval_symbol = df_heatmap_pval_symbol.replace('**', 2)
    #df_heatmap_pval_symbol = df_heatmap_pval_symbol.replace('***', 3)
    #df_heatmap_pval_symbol = df_heatmap_pval_symbol.replace('****', 4)
    # df_heatmap_pval_symbol = df_heatmap_pval_symbol.replace('', 1)

    df_heatmap_pval_symbol.to_csv(args.path + '/plots/df_feat_imp_p_val_' + perf_type + '.csv')

    df_heatmap_pval_symbol[df_heatmap_pval_symbol == 0] = ''
    fig, ax = plt.subplots(figsize=(10, 10))
    sn.heatmap(df_heatmap_pval, fmt="", cbar=True, linewidth=4, ax=ax)
    sn.heatmap(df_heatmap_pval, annot=df_heatmap_pval_symbol, annot_kws={'va': 'top', 'size': '18'}, fmt="", cbar=False)
    plt.tick_params(axis='both', which='major', labelsize=14, labelrotation=45,
                    labelbottom=False, bottom=False, top=False, labeltop=True)
    plt.yticks(rotation=0)
    plt.tight_layout()
    fig.savefig(args.path + '/plots/heatmap_importance_gen_across_states_' + perf_type + '.png')
    plt.close()

    return None


'''
Feature ranking
'''
def rank_perm_feature_imp():
    output_files = glob.glob(args.path + "/permut_feat_imps/compose/*_all_feat_*" + args.surv_model + "*.csv")
    data_list = []
    for filename in output_files:
        df = pd.read_csv(filename, index_col=None, header=0)
        data_list.append(df)
        df['fold'] = filename.split('/')[-1].split('_')[1]

    concat_feats = pd.concat(data_list, axis=0, ignore_index=True)

    concat_feats = concat_feats[['group', 'class_generic', 'weight', 'std', 'fold', 'column_id']]

    overall_fea_class = pd.read_excel(args.legend)
    concat_feats = pd.merge(concat_feats, overall_fea_class[['column_id', 'display_id']],
                            on='column_id', how='inner')

    unique_feat_groups = concat_feats['class_generic'].unique().tolist()
    color = ['deepskyblue', 'pink', 'darkorange', 'lightgreen', 'red', 'magenta']
    concat_feats['color'] = ''
    for i, item in enumerate(unique_feat_groups):
        print(item, color[i])
        concat_feats['color'].mask(concat_feats['class_generic'] == item, color[i], inplace=True)

    concat_feats['display_id'] = concat_feats['display_id'].str.replace(r'(CNV_)|(CNV\.SNV_)|(SNV_)', '', regex=True)

    unique_states = concat_feats['group'].unique()
    for state in unique_states:
        print('state is', state)
        df_inst = concat_feats[(concat_feats['group'] == state)]

        df_inst['std'] = df_inst['std'] / (df_inst[['fold']].max(axis=1) + 1) # np.sqrt(df_inst[['fold']].max(axis=1) + 1)

        df_inst = df_inst.groupby(['column_id', 'display_id', 'color'])['weight', 'std'].mean().reset_index()
        df_inst = df_inst.sort_values(by='weight', ascending=False)

        df_inst.to_csv(args.path + '/permut_feat_imps/aggregate/' + state + '_' + args.surv_model + '_' + 'folds_feat_imp.csv', index=False)

        df_inst = df_inst.iloc[:50] # top 50 feats

        ### supplementary
        fig, ax = plt.subplots(figsize=(10, 34))
        sn.barplot(y="display_id", x="weight", orient='h', xerr= df_inst['std'] * 1, palette=df_inst['color'], data=df_inst)
        ax.set_ylabel("Feature", fontsize=30, fontname="Arial")
        ax.set_xlabel("Score", fontsize=30, fontname="Arial")
        ax.tick_params(axis='both', which='major', labelsize=28)
        fig.suptitle(state + '\n', fontsize=30, fontname="Arial")
        fig.tight_layout()
        fig.savefig(args.path + '/plots/agg_perm_imp_' + state + '_' + args.surv_model + '_' + 'folds_feat_imp.png')


def get_reranked_perm_features(gp, df_feat_uni):
    df_feat_rank = pd.read_csv(args.path + '/permut_feat_imps/aggregate/' + gp +
                                    '_' + args.surv_model + '_folds_feat_imp.csv')

    keep_feats = df_feat_uni[~ ((df_feat_uni['PFS.pvalue'] <= 0.01) & (df_feat_uni['OS.pvalue'] <= 0.01)) ]

    df_feat_rank = df_feat_rank[df_feat_rank['column_id'].isin(keep_feats['Genomic.Feature'].tolist())]

    df_feat_rank = df_feat_rank[df_feat_rank['weight'] > 0]

    # df_feat_rank = df_feat_rank[~df_feat_rank['column_id'].isin(df_feat_uni['Genomic.Feature'].tolist())]

    df_feat_rank.to_csv(args.path + '/permut_feat_imps/aggregate_reranked/' + gp +
                                    '_' + args.surv_model + '_folds_feat_imp.csv', index=False)

    return None


def ABS_SHAP(df_shap, df, group):
    shap_v = pd.DataFrame(df_shap)
    feature_list = df.columns
    shap_v.columns = feature_list
    df_v = df.copy().reset_index().drop('index', axis=1)

    print('first', df_shap.shape, df.shape)

    # Determine the correlation in order to plot with different colors
    corr_list = list()
    for i in feature_list:
        b = np.corrcoef(shap_v[i], df_v[i])[1][0]
        corr_list.append(b)
    corr_df = pd.concat([pd.Series(feature_list), pd.Series(corr_list)], axis=1).fillna(0)
    # Make a data frame. Column 1 is the feature, and Column 2 is the correlation coefficient
    corr_df.columns = ['Variable', 'Corr']

    colors_to_use = ['lightsalmon', 'skyblue']
    corr_df['Sign'] = np.where(corr_df['Corr'] > 0, colors_to_use[0], colors_to_use[1])

    shap_abs = np.abs(shap_v)
    k = pd.DataFrame(shap_abs.mean()).reset_index()
    k.columns = ['Variable', 'SHAP_abs']

    k2 = k.merge(corr_df, left_on='Variable', right_on='Variable', how='inner')
    k2 = k2.sort_values(by='SHAP_abs', ascending=True)

    # clean up any duplicate variables
    # k2 = k2[~k2['Variable'].str.contains('.1$')]

    # calculate 95% CI
    k2_std = pd.DataFrame(shap_abs.std()).reset_index()
    k2_std.columns = ['Variable', 'SHAP_std']

    k2 = k2.merge(k2_std, left_on='Variable', right_on='Variable', how='inner')

    var_actual_count = df.astype(bool).sum(axis=0).reset_index().rename(
        columns={'index': 'Variable', 0: 'count_actual_observed'})

    k2 = k2.merge(var_actual_count, left_on='Variable', right_on='Variable', how='inner')

    k2 = k2[k2['count_actual_observed'] != 0]

    stats = k2[['SHAP_abs', 'count_actual_observed', 'SHAP_std']]

    ci95_hi = []
    ci95_lo = []

    for i in stats.index:
        m, c, s = stats.loc[i]
        # print(i, m,c,s)
        ci95_hi.append(m + 1.96 * s / math.sqrt(c))
        ci95_lo.append(m - 1.96 * s / math.sqrt(c))

    k2['ci95_hi'] = ci95_hi
    k2['ci95_lo'] = ci95_lo

    k2 = k2[(k2['ci95_lo'] > 0)].sort_values(by='SHAP_abs', ascending=False)
    # k2.to_csv(args.path + '/ci_shap/' + group + '_ci.csv', index=False)

    k2 = k2.reset_index().drop(columns='index')
    k2_cp2 = k2.copy()

    # read feats and their attributes
    overall_fea_class = pd.read_excel(args.legend)

    # k2[::-1][['Variable', 'SHAP_abs']].to_csv(args.path + '/ci_shap/' + group + '_shap_features.csv', index=False)

    #### top 45
    k2 = k2.reset_index().drop(columns='index').iloc[:45]
    k2 = k2.reset_index().drop(columns='index')

    k2 = pd.merge(k2, overall_fea_class[['column_id', 'display_id']], left_on='Variable', right_on='column_id')

    '''
    fig, ax = plt.subplots(figsize=(8, 12))
    k2.plot.barh(x='feature_renamed', y='SHAP_abs', color=k2['Sign'].tolist(), ax=ax, legend=False)
    ax.set_xlabel('Importance Value', fontsize=11, fontname='Arial')
    ax.set_ylabel('Feature', fontsize=11)
    ax.set_title(group)
    plt.yticks(fontsize=10)
    plt.xticks(fontsize=10)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    fig.savefig(args.path + '/plots/' + group + '_bar_top.png', dpi=600)
    '''

    # grant / paper
    k2.loc[k2['Variable'] == 'gender', 'display_id'] = 'Female'

    print(k2)
    fig, ax = plt.subplots(figsize=(12, 8))
    k2.plot.bar(x='display_id', y='SHAP_abs', color=k2['Sign'].tolist(), ax=ax, legend=False)
    ax.set_ylabel('Importance Value', fontsize=14, fontname='Arial')
    ax.set_xlabel('Feature', fontsize=14, fontname='Arial')

    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    plt.tight_layout()
    fig.savefig(args.path + '/plots/' + group + '_bar_top.pdf', dpi=600)

    df_feat_imp = pd.merge(k2_cp2, overall_fea_class[['column_id', 'display_id', 'class_generic', 'class_granular']],
                           left_on='Variable', right_on='column_id', how='left')

    df_feat_imp_orig = df_feat_imp.copy()

    #### pie
    df_feat_imp = pd.merge(k2_cp2, overall_fea_class[['column_id', 'display_id', 'class_generic']],
                           left_on='Variable', right_on='column_id', how='left')

    size_df = df_feat_imp.groupby(['class_generic']).size()
    size_df = pd.DataFrame(size_df, columns=['count'])
    feat_df = df_feat_imp.groupby(['class_generic']).mean('SHAP_abs')
    feat_size_df = pd.merge(feat_df, size_df, on=['class_generic'])

    # weight
    size_df['reweighted'] = 1 - (size_df['count'] / size_df['count'].sum())
    size_df = size_df.sort_values(by='reweighted')

    feat_size_df = pd.merge(df_feat_imp, size_df, on=['class_generic'])
    feat_size_df = feat_size_df.sort_values(by=['SHAP_abs'], ascending=False).reset_index(drop=True)
    feat_size_df['updated_weight'] = (
                size_df['count'].sum() - feat_size_df.reset_index()['index']) # * feat_size_df['reweighted']

    feat_size_df['SHAP_abs'] = feat_size_df['updated_weight'] * feat_size_df['SHAP_abs']
    feat_size_df = feat_size_df.groupby(['class_generic']).sum('SHAP_abs')

    fig, ax = plt.subplots()
    pie_plot_feat_df = feat_size_df.reset_index()
    pie_plot_feat_df['name_var'] = pie_plot_feat_df.reset_index()['Class']
    pie_plot_feat_df['color'] = np.select(
        [(pie_plot_feat_df['class_generic'] == 'Clinical'), (pie_plot_feat_df['class_generic'] == 'SCT'),
         (pie_plot_feat_df['class_generic'] == 'Genomics'), (pie_plot_feat_df['class_generic'] == 'Treatment'),
         (pie_plot_feat_df['class_generic'] == 'Continuous treatment')
         ],
        ['lightseagreen', 'pink', 'orange', 'khaki', 'lightcoral'])

    pie_plot_feat_df.plot(kind='pie', y='SHAP_abs', normalize=True, autopct='%.1f%%',
                          labels=pie_plot_feat_df['name_var'], colors=pie_plot_feat_df['color'],
                          ax=ax, textprops={'fontsize': 15, 'fontname': "Arial"}, labeldistance=1.1, radius=0.8)

    ax.set_ylabel('')
    legend = ax.legend(loc='lower left')
    legend.remove()
    plt.tight_layout()
    fig.savefig(args.path + '/plots/' + group + '_pie.pdf')

    #### pie granular
    size_df = df_feat_imp_orig.groupby(['class_granular']).size()
    size_df = pd.DataFrame(size_df, columns=['count'])

    # weight
    size_df['reweighted'] = 1 - (size_df['count'] / size_df['count'].sum())
    size_df = size_df.sort_values(by='reweighted')

    feat_size_df = pd.merge(df_feat_imp_orig, size_df, on=['class_granular'])
    feat_size_df = feat_size_df.sort_values(by=['SHAP_abs'], ascending=False).reset_index(drop=True)
    feat_size_df['updated_weight'] = (
                size_df['count'].sum() - feat_size_df.reset_index()['index']) # * feat_size_df['reweighted']

    feat_size_df['SHAP_abs'] = feat_size_df['updated_weight'] * feat_size_df['SHAP_abs']
    feat_size_df = feat_size_df.groupby(['class_granular']).sum('SHAP_abs')

    fig, ax = plt.subplots()
    pie_plot_feat_df = feat_size_df.reset_index()
    pie_plot_feat_df['name_var'] = pie_plot_feat_df.reset_index()['Class_granular']
    pie_plot_feat_df['color'] = np.select(
        [(pie_plot_feat_df['class_granular'] == 'Clinical'), (pie_plot_feat_df['class_granular'] == 'SCT'),
         (pie_plot_feat_df['class_granular'] == 'CNV'), (pie_plot_feat_df['class_granular'] == 'Translocation'),
         (pie_plot_feat_df['class_granular'] == 'Mutations'), (pie_plot_feat_df['class_granular'] == 'CNV_Mutations'),
         (pie_plot_feat_df['class_granular'] == 'Treatment'),
         (pie_plot_feat_df['class_granular'] == 'Continuous treatment')
         ],
        ['lightseagreen', 'pink', 'rosybrown', 'gray', 'lightgray', 'magenta', 'khaki', 'lightcoral'])

    pie_plot_feat_df.plot(kind='pie', y='SHAP_abs', normalize=True, autopct='%.1f%%',
                          labels=pie_plot_feat_df['name_var'], colors=pie_plot_feat_df['color'],
                          ax=ax, textprops={'fontsize': 15, 'fontname': "Arial"}, labeldistance=1.1, radius=0.8)

    ax.set_ylabel('')
    legend = ax.legend(loc='lower left')
    legend.remove()
    plt.tight_layout()
    fig.savefig(args.path + '/plots/' + group + '_pie_granular.pdf')

    return None


def plot_shap_values_per_group(group):
    shap_output_files = glob.glob(args.path + '/shap_feat_imps/fold_*_' + group + '_*.csv')

    shap_output_files.sort()

    df_concat_shap_values = pd.DataFrame()
    for filename in shap_output_files:
        single_df = pd.read_csv(filename, index_col=None, header=0)
        df_concat_shap_values = df_concat_shap_values.append(single_df, ignore_index=True)

    test_files = glob.glob(args.path + '/shap_feat_imps/test_data_fold_*_' + group + '_*.csv')

    test_files.sort()
    df_concat_shap_test = pd.DataFrame()
    for filename in test_files:
        single_df = pd.read_csv(filename, index_col=None, header=0)
        df_concat_shap_test = df_concat_shap_test.append(single_df, ignore_index=True)

    print(df_concat_shap_values.shape, df_concat_shap_test.shape)
    fig = plt.figure()
    shap.summary_plot(df_concat_shap_values.values, df_concat_shap_test.values, list(df_concat_shap_test), show=False, max_display=45)
    plt.tight_layout()
    fig.savefig(args.path + '/plots/' + group + '_scatter_summary_plot.png', dpi=600)

    # correlation of feat imp in bar plot
    ABS_SHAP(df_concat_shap_values, df_concat_shap_test, group)

    # df_concat_shap_values.mean(axis=0)

    fig = plt.figure()
    shap.summary_plot(df_concat_shap_values.values, df_concat_shap_test.values, list(df_concat_shap_test), plot_type='violin', show=False, max_display=45)
    plt.tight_layout()
    fig.savefig(args.path + '/plots/' + group + '_shap_violin_summary_plot.png', dpi=600)

    return None


def plot_shap_values():
    for group in ['Progress (P1)', 'Progress & deceased (P1)',
              'Progress (P2)', 'Progress & deceased (P2)']:
        plot_shap_values_per_group(group)


'''
survival/event curves by model
'''
def get_frac_of_pats():
    if args.ext_val:
        output_files = glob.glob(args.path + "/kfold/all_states/testGroup*~fold~0~*")
    else:
        output_files = glob.glob(args.path + "/loo/state_probs/*,orig,*")
        print(output_files)

    output_files = [item for item in output_files if 'P2' in item]

    data_list = []
    for filename in output_files:
        single_df = pd.read_csv(filename, index_col=None, header=0)
        data_list.append(single_df.head(args.preds_time_in_days))

    df = pd.concat(data_list, axis=1)
    # time_df = df['Time (Days)']
    df = df.drop(['sample', 'Time (Days)'], axis=1)

    df_accum = pd.DataFrame()
    for ind, group in enumerate(df.columns.unique()):
        print('group', group)
        spec_cols = [col for col in df.columns if group == col]
        df_spec = df[spec_cols]
        df_mean = df_spec.mean(axis=1).abs()
        df_mean = df_mean.reset_index()
        df_mean = df_mean.drop('index', axis=1)

        df_mean.columns = [group]

        df_accum = pd.concat((df_accum, df_mean), axis=1)

    # print(time_df, type(time_df))
    df_accum.index = df.index/365
    df_accum.columns = ['Alive in Induction', 'Alive after progression (Induction)', 'Alive in Phase 2',
            'Alive after Progression (Induction)', 'Death after progression (Induction)',
                        'Death in Phase 2', 'Death after Progression (Phase 2)']

    df_accum.to_csv(args.path + '/plots/mean_frac_of_patients_ext_val_' + str(args.ext_val) + '_at' + str(args.preds_time_in_days) + '.csv', index=False)


def plot_frac_of_pats():
    df = pd.read_csv(args.path + '/plots/mean_frac_of_patients_ext_val_' + str(args.ext_val) + '_at' + str(args.preds_time_in_days) + '.csv')
    df.columns = ['Non-POD \nin induction', 'POD in induction and alive',
                  'Non-POD and alive', 'POD after induction and alive', 'POD in induction and death',
                  'Non-POD and death', 'POD after induction and death']

    df = df[['Non-POD \nin induction', 'POD in induction and alive', 'POD in induction and death',
               'Non-POD and alive', 'POD after induction and alive',
               'Non-POD and death', 'POD after induction and death']]

    df.index = df.index / 365

    fig, ax = plt.subplots(figsize=(6,6))
    df.plot.area(stacked=True, color=('lightgrey', 'lightseagreen', 'beige', 'pink', 'lightblue', 'gold', 'lightcoral'), ax=ax)

    xy_annotate = [(0.04, 0.02), (1.8, 0.03), (2, 0.12), (2.5, 0.4), (3.15, 0.7), (3.34, 0.9), (3.1, 0.96)]

    for ind, item in enumerate(list(df)):
        if ind == 7:
            ax.annotate(item, xy=(xy_annotate[ind][2], xy_annotate[ind][3]),
                        xytext=(xy_annotate[ind][0], xy_annotate[ind][1]), va='center',
                        arrowprops=dict(facecolor='black', arrowstyle="->"), fontsize='small')
        else:
            ax.annotate(item, xy=(xy_annotate[ind][0], xy_annotate[ind][1]), fontsize=8, fontname='Arial')

    ax.set_xlabel('Time from start of induction (years)', fontsize=12, fontname="Arial")
    ax.set_ylabel('Fraction of patients', fontsize=12, fontname="Arial")
    ax.tick_params(axis='both', which='major', labelsize=10)
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    legend = ax.legend(loc='lower left')
    legend.remove()
    plt.xlim(xmin=0, xmax=5)
    plt.ylim(ymax=1)
    fig.tight_layout()
    fig.savefig(args.path + '/plots/mean_frac_of_patients_ext_val_' + str(args.ext_val) + '_at' + str(args.preds_time_in_days) + '.png')


'''
Outputs
'''
def get_risks_kfold_model(val, agg_type='mean'):
    if val == 'True':
        path = args.path + "/kfold/model_probs/test_group_9*neural_cox_non_prop"
    else:
        path = args.path + "/kfold/model_probs/group_9*neural_cox_non_prop"

    output_files_prp1 = glob.glob(path + '*_Progress (P1).csv')
    output_files_pr_dec_p1 = glob.glob(path + '*_Progress & deceased (P1).csv')
    output_files_prp2 = glob.glob(path + '*_Progress (P2).csv')
    output_files_pr_dec_p2 = glob.glob(path + '*_Progress & deceased (P2).csv')

    df_concat = pd.DataFrame()
    for file_ in [output_files_prp1, output_files_pr_dec_p1, output_files_prp2, output_files_pr_dec_p2]:
        df_concat_per_fold = pd.DataFrame()
        for filename in file_:
            df_file = pd.read_csv(filename, index_col=None, header=0)
            df_file = df_file.rename(columns=lambda x: re.sub('\.', '-', x))
            single_df = df_file.iloc[args.preds_time_in_days-1].rename(filename.split('/')[-1].split('_')[-1].split('.')[0]).to_frame()
            df_concat_per_fold = pd.concat([df_concat_per_fold, single_df])

        if agg_type == 'mean':
            df_concat_per_fold = df_concat_per_fold.groupby(level=0).mean()
        elif agg_type == 'median':
            df_concat_per_fold = df_concat_per_fold.groupby(level=0).median()

        df_concat = pd.concat([df_concat, df_concat_per_fold], axis=1)

    df_concat = df_concat.sort_index().reset_index().rename(columns={'index': 'sample'})
    df_concat.to_csv(args.path + '/preds_treatment/kfold_' + agg_type + '_' +
                    args.surv_model + '_model_risks_ext_val_' + str(val) + '_at_' + str(
        args.preds_time_in_days) + '.csv', index=False)

    return None


def get_risks_kfold_multistate(val, df, agg_type='mean'):
    if val == 'True':
        output_files = glob.glob(args.path + "/kfold/state_probs/all_states/testGroup~hd~*.csv")
    else:
        output_files = glob.glob(args.path + "/kfold/state_probs/all_states/group~9*neural_cox_non_prop.csv")

    df_years = pd.DataFrame([])
    for sampl in df['sample']:
        if ('I-H-' in sampl) or ('K08K-' in sampl):
            sampl_orig = sampl
            sampl = re.sub('-', '.', sampl)
        else:
            sampl_orig = sampl

        if val == 'True':
            output_files_cp = [item for item in output_files if sampl == item.split('/')[-1].split('~')[2]]
        else:
            output_files_cp = [item for item in output_files if sampl == item.split('/')[-1].split('~')[3]]

        if len(output_files_cp) > 0:
            df_concat_per_sample = pd.DataFrame()
            for filename in output_files_cp:
                single_df = pd.read_csv(filename, index_col=None, header=0)
                single_df_years = single_df[single_df['Time (Days)'] == args.preds_time_in_days-1].drop('Time (Days)', axis=1)
                single_df_years['sample'] = sampl_orig
                single_df_years.set_index('sample', inplace=True)
                df_concat_per_sample = pd.concat([df_concat_per_sample, single_df_years])

            if agg_type == 'mean':
                df_concat_per_sample = df_concat_per_sample.groupby(level=0).mean()
            elif agg_type == 'median':
                df_concat_per_sample = df_concat_per_sample.groupby(level=0).median()

            df_years = pd.concat([df_years, df_concat_per_sample])

    print(df_years)
    df_years = df_years.sort_index().reset_index().rename(columns={'index': 'sample'})
    df_years.to_csv(args.path + '/preds_treatment/kfold_' + agg_type + '_' +
                    args.surv_model + '_multistate_risks_ext_val_' + str(val) + '_at_' + str(args.preds_time_in_days) + '.csv', index=False)

    return None


def get_risks_loo_model():
    output_files = glob.glob(args.path + "/loo/model_probs/*orig*neural_cox_non_prop.csv")

    df_years = pd.DataFrame([])

    for filename in output_files:
        single_df = pd.read_csv(filename, index_col=None, header=0)
        single_df_years = single_df.iloc[args.preds_time_in_days-1]
        df_years = df_years.append(single_df_years, ignore_index=True)

    df_years = df_years.drop_duplicates(subset=['sample'], keep='first')

    df_years = df_years[['sample', 'Progress (P1)', 'Progress & deceased (P1)', 'Progress (P2)', 'Progress & deceased (P2)']]

    df_years.to_csv(args.path + '/preds_treatment/loo_' + args.surv_model + '_model_risks_at' + str(args.preds_time_in_days) + '.csv', index=False)

    return None


def get_risks_loo_multistate():
    output_files = glob.glob(args.path + "/loo/state_probs/*orig*neural_cox_non_prop.csv")

    df_years = pd.DataFrame([])

    for filename in output_files:
        single_df = pd.read_csv(filename, index_col=None, header=0)
        single_df_years = single_df.iloc[args.preds_time_in_days-1]
        df_years = df_years.append(single_df_years, ignore_index=True)

    df_years = df_years.drop_duplicates(subset=['sample'], keep='first')

    df_years = df_years[['sample', 'Alive in Induction', 'Alive in POD (induction)',
                           'Alive in Phase 2', 'Alive in POD', 'Death in POD (induction)',
                           'Death in Phase 2', 'Death in POD (Phase 2)']]

    df_years.to_csv(args.path + '/preds_treatment/loo_' + args.surv_model + '_multistate_risks_at' + str(args.preds_time_in_days) + '.csv', index=False)

    return None


'''
Scorer to understand overall feature group importance
'''
def scorer(base, comp):
    score = comp['median'] - (base['median'])
    score = score /base['median']
    score = score * ( 1 +
            (comp['Group Feat. size']/comp['Accum. Feat. size']) - (base['Group Feat. size']/base['Accum. Feat. size'])
            )
    return score * 100


def relative_importance(df_sub_c_ind_vals):
    df_median = df_sub_c_ind_vals.median().sort_values(ascending=True)

    # [orig size, previous feature set + current size]
    feat_size = {
        'Plus_Top genomics': [17,34],
        'Plus_Cont. treat': [1,11], 'Plus_Therapy': [11,17], 'Plus_SCT': [1,10],
        'All features': [153,153], 'Plus_Rec. Transloc': [3,9],
        'Plus_CLIN': [5, 6], 'ISS': [1,1], 'R2-ISS': [1,1], 'R-ISS': [1,1],
         'ISS+SCT': [1, 2],'ISS+Cont. treat': [1, 2], 'ISS+Therapy': [1, 7],
        'ISS+TopGen': [1, 18], 'ISS+Rec. Transloc': [1, 4]
         }

    df_median = pd.DataFrame(df_median).rename(columns={0: 'median'})
    rename_df = pd.DataFrame.from_dict(feat_size, orient='index')\
                        .reset_index().rename(columns={'index': 'group',
                              0: 'Group Feat. size', 1: 'Accum. Feat. size'})

    df_median = pd.merge(df_median, rename_df, left_index=True, right_on='group')
    df_median = df_median[(df_median['group'] != 'R-ISS') & (df_median['group'] != 'All features') &
                                (df_median['group'] != 'R2-ISS')]
    df_median = df_median.set_index('group')

    df_median_cp = df_median.copy()

    df_score = pd.DataFrame(columns=['score'])

    # R-ISS -> ISS
    #df_score.loc['R-ISS => ISS'] = scorer(df_median_cp.loc['R-ISS'], df_median_cp.loc['ISS'])

    # R2-ISS -> ISS
    #df_score.loc['R2-ISS => ISS'] = scorer(df_median_cp.loc['R2-ISS'], df_median_cp.loc['ISS'])

    '''
    # ISS -> CLIN
    df_score.loc['ISS => ISS+CLIN'] = scorer(df_median_cp.loc['ISS'], df_median_cp.loc['Plus_CLIN'])
    
    df_score.loc['ISS => ISS+CLIN+Rec.Transloc.'] = scorer(df_median_cp.loc['ISS'], df_median_cp.loc['Plus_Rec. Transloc'])

    df_score.loc['ISS => ISS+CLIN+Rec.Transloc.+SCT'] = scorer(df_median_cp.loc['ISS'], df_median_cp.loc['Plus_SCT'])

    df_score.loc['ISS => ISS+CLIN+Rec.Transloc.+SCT+Maintenance'] = scorer(df_median_cp.loc['ISS'], df_median_cp.loc['Plus_Cont. treat'])

    df_score.loc['ISS => ISS+CLIN+Rec.Transloc.+SCT+Maintenance+Therapy'] = scorer(df_median_cp.loc['ISS'], df_median_cp.loc['Plus_Therapy'])

    df_score.loc['ISS => ISS+CLIN+Rec.Transloc.+SCT+Maintenance+Therapy+Top genomics'] = scorer(df_median_cp.loc['ISS'], df_median_cp.loc['Plus_Top genomics'])

    df_score.loc['ISS => ISS+SCT'] = scorer(df_median_cp.loc['ISS'], df_median_cp.loc['ISS+SCT'])

    df_score.loc['ISS => ISS+Maintenance'] = scorer(df_median_cp.loc['ISS'], df_median_cp.loc['ISS+Cont. treat'])

    df_score.loc['ISS => ISS+Therapy'] = scorer(df_median_cp.loc['ISS'], df_median_cp.loc['ISS+Therapy'])

    df_score.loc['ISS => ISS+TopGen'] = scorer(df_median_cp.loc['ISS'], df_median_cp.loc['ISS+TopGen'])

    df_score.loc['ISS => ISS+Rec.Transloc'] = scorer(df_median_cp.loc['ISS'], df_median_cp.loc['ISS+Rec. Transloc'])
    '''

    # CLIN -> REC. transl
    df_score.loc['ISS+CLIN => ISS+CLIN+Rec.transloc.'] = scorer(df_median_cp.loc['Plus_CLIN'], df_median_cp.loc['Plus_Rec. Transloc'])

    # REC. transl -> SCT
    df_score.loc['ISS+CLIN+Rec.transloc. => ISS+CLIN+Rec.Transloc.+SCT'] = scorer(df_median_cp.loc['Plus_Rec. Transloc'], df_median_cp.loc['Plus_SCT'])

    # SCT -> Therapy
    df_score.loc['ISS+CLIN+Rec.Transloc.+SCT => ISS+CLIN+Rec.Transloc.+SCT+Maintenance'] = scorer(df_median_cp.loc['Plus_SCT'], df_median_cp.loc['Plus_Cont. treat'])

    # Therapy -> Cont. treat
    df_score.loc['ISS+CLIN+Rec.Transloc.+SCT+Maintenance => ISS+CLIN+Rec.Transloc.+SCT+Maintenance+Therapy'] = scorer(df_median_cp.loc['Plus_Cont. treat'], df_median_cp.loc['Plus_Therapy'])

    # Cont. treat -> Genomics
    df_score.loc['ISS+CLIN+Rec.Transloc.+SCT+Maintenance+Therapy => ISS+CLIN+Rec.Transloc.+SCT+Maintenance+Therapy+Top genomics'] = scorer(df_median_cp.loc['Plus_Cont. treat'], df_median_cp.loc['Plus_Top genomics'])

    df_score.to_csv(args.path + '/plots/' + args.outcome + '_relative_imp.csv')

    fig, ax = plt.subplots(figsize=(12, 3))
    sn.set_theme(style="white")
    df_score = df_score.sort_values(by='score', ascending=False).reset_index()
    sn.barplot(x='score', y='index', data=df_score, palette='coolwarm', orient='h', ax=ax)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.tick_params(which='both', width=2)
    ax.tick_params(which='major', length=7)
    ax.tick_params(which='minor', length=4)

    ax.tick_params(axis='both', which='major', labelsize=6)

    plt.xlabel('Relative Feature Group Impact Score (%)', fontsize=10, fontname='Arial')
    plt.ylabel('')
    ax.tick_params(axis='both', which='major', labelsize=8)

    plt.axvline(x=0)

    plt.tight_layout()
    fig.savefig(args.path + '/plots/' + args.outcome + '_relative_imp.png')

    # plot 2
    fig, ax = plt.subplots(figsize=(6, 3))
    df_median_cp[['Group Feat. size', 'Accum. Feat. size']].sort_values(by='Accum. Feat. size').plot(kind='barh', stacked=True, color=['lightblue', 'pink'], ax=ax)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.grid(True, which='major', axis='x')
    plt.ylabel('')
    ax.tick_params(axis='both', which='major', labelsize=8)
    plt.xlabel('# Features', fontsize=10, fontname='Arial')
    plt.tight_layout()
    fig.savefig(args.path + '/plots/' + args.outcome + '_num_features.png')
    return None


def test_graph():
    df_loo_ours = pd.read_csv(args.path + '/preds_treatment/loo_neural_cox_non_prop_multistate_risks_at1825.csv')
    df_cv_ours = pd.read_csv(args.path + '/preds_treatment/kfold_median_neural_cox_non_prop_multistate_risks_ext_val_False_at_1825.csv')
    df_hd = pd.read_csv(args.path + '/preds_treatment/kfold_median_neural_cox_non_prop_multistate_risks_ext_val_True_at_1825.csv')

    #
    df_cv_ours['risk_os'] = 1 - (
                df_cv_ours['Alive in Induction'] + df_cv_ours['Alive in POD (induction)'] +
                df_cv_ours['Alive in Phase 2'] + df_cv_ours['Alive in POD'])
    df_cv_ours['risk_pfs'] = 1 - df_cv_ours['Alive in Phase 2']

    #
    df_loo_ours['risk_os'] = 1 - (
            df_loo_ours['Alive in Induction'] + df_loo_ours['Alive in POD (induction)'] +
            df_loo_ours['Alive in Phase 2'] + df_loo_ours['Alive in POD'])
    df_loo_ours['risk_pfs'] = 1 - df_loo_ours['Alive in Phase 2']

    #
    df_hd['risk_os'] = 1 - (df_hd['Alive in Induction'] + df_hd['Alive in POD (induction)'] +
                            df_hd['Alive in Phase 2'] + df_hd['Alive in POD'])
    df_hd['risk_pfs'] = 1 - df_hd['Alive in Phase 2']

    return None


def main():
    global args
    args = parser.parse_args()

    df_val = load_dataset(val=True)
    df = load_dataset(val=False)

    # Fig. state c-index
    if args.c_index_state:
        plot_c_index_state('val')

    # Fig. multistate c-index
    if args.c_index_group_addition:
        print('multistate c-index')
        if args.ext_val is not True:
            df_mul = multistate(df)
        else:
            df_mul = multistate(df_val)
        plot_metrics_multistate_models(df_mul)

    # Fig. state rank across folds
    if args.rank_folds and args.c_index_state:
        rank_folds_first_level() # state

    if args.rank_folds:
        df_sub_c_ind_vals = rank_folds_multistate() # multistate

    # Save aggregate rank scores
    # Fig. Rank scores by state
    if args.agg_perm_imp:
        rank_perm_feature_imp()

    # Fig. Incremental gen; train/val
    if args.get_top_gen_feats:
        top_gen_feats_train = plot_c_index_model_bygroups('train')
        # top_gen_feats_val = plot_c_index_model_bygroups('val')

    # Fig. heatmap of importance of the top gen feats
    if args.heatmap_bygroup:
        heatmap_bygroups_imp_feats(top_gen_feats_train, 'train')
        # heatmap_bygroups_imp_feats(top_gen_feats_val, 'val')

    # Fig. SHAP scores
    if args.shap:
        plot_shap_values()

    # Fig. Sediment plot
    if args.get_frac_plot:
        get_frac_of_pats()
        plot_frac_of_pats()

    if args.re_rank_feats and not args.get_top_gen_feats:
        print('re-rank-feats in metrics file')
        df_feat_uni_ = pd.read_csv(args.in_data_file + '/Genomic_features_PFS_OS.txt', sep='\t')
        df_feat_uni_ = df_feat_uni_[(df_feat_uni_['PFS.pvalue'] < 0.05) & (df_feat_uni_['OS.pvalue'] < 0.05)]
        df_feat_uni_.to_csv(args.path + '/plots/df_feat_uni_p_val_05.csv', index=False)

        df_feat_uni_p_05 = df_feat_uni_.copy()

        df_feat_uni_ = df_feat_uni_[(df_feat_uni_['PFS.pvalue'] < 0.01) & (df_feat_uni_['OS.pvalue'] < 0.01)]
        df_feat_uni_.to_csv(args.path + '/plots/df_feat_uni.csv', index=False)

        for gp in ['Progress (P1)', 'Progress & deceased (P1)', 'Move to phase 2',
                   'Progress (P2)', 'Progress & deceased (P2)', 'Non-progress & deceased (P2)']:
            get_reranked_perm_features(gp, df_feat_uni_p_05)

    if args.get_risk_scores:
        get_risks_kfold_model('False', 'mean')  # cv ours
        # get_risks_kfold_model('True', 'mean')  # hd
        # get_risks_kfold_multistate('False', df, 'mean')  # cv ours
        # get_risks_kfold_multistate('True', df_val, 'mean')  # hd
        # get_risks_loo_multistate()  # loo ours

    if args.relative_imp:
        relative_importance(df_sub_c_ind_vals)

    if args.test_graph:
        test_graph()

    return None


if __name__ == '__main__':
    main()
