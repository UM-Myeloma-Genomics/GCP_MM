import pandas as pd
import argparse
import glob as glob
import re
import seaborn as sn
import matplotlib.pyplot as plt
import numpy as np
from itertools import combinations
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.mixture import GaussianMixture
from sklearn import metrics
from sklearn.manifold import TSNE
import scipy.stats
import plotly.express as px
import matplotlib as mpl
from functools import reduce

parser = argparse.ArgumentParser(description='Prediction 1.0 Multiple Myeloma workflow')
parser.add_argument('--in_data_file', type=str, default='/Users/axr2376/Desktop/pred_1_0_paper/data_in/', help='input dataset path')
parser.add_argument('--path', type=str, default='/Users/axr2376/Desktop/pred_1_0_paper/data_out/expt_1', help='input dataset path')
parser.add_argument('--legend', type=str, default='/Users/axr2376/Desktop/pred_1_0_paper/data_in/legend_prediction_may_2020.xlsx', help='input dataset path')
parser.add_argument('--surv_model', type=str, default='neural_cox_non_prop', help='model to train to survival model | rsf, cph, svms')
parser.add_argument('--preds_time_in_days', type=int, default=1825, help='time in days')
parser.add_argument('--thresh_m_to_p2', type=int, default=365, help='days')

parser.add_argument('--outcome', type=str, default='risk_pfs', help='none in the case of first level comparisons, os, pfs')
parser.add_argument('--pred_risk_at', type=int, default=1825, help='time in days')
parser.add_argument('--get_treat_stats', action='store_true', default=False, help='get treat stats from dataset')
parser.add_argument('--sample', type=str, default='PD5886a', help='SAMPLE ids: MMRF_1332, PD5886a, MRC_1018')
parser.add_argument('--plot_treat_combo_indiv_pat', action='store_true', default=False, help='indiv. predictions')
parser.add_argument('--top_display', type=int, default=8, help='time in days')
parser.add_argument('--plot_treat_combo_all_pats', action='store_true', default=False, help='summarize all predictions')
parser.add_argument('--ext_val', action='store_true', default=False, help='external test set')


def data_filter(df):
    df = df[df['duration'] >= 0]
    df = df.dropna(subset=['os_time', 'pfs_time', 'pfs_code', 'os_code'])

    df['time_SCT'].fillna(0, inplace=True)

    if args.ext_val is not True: # not for independent Heidelberg cohort
        df_r_iss = pd.read_csv(args.in_data_file + '/ISS_RISS_R2ISS_LDH_all_cohort.txt', sep='\t')
        df_r_iss['sample'] = np.where(df_r_iss['study'] == 'MRC_XI', 'MRC_' + df_r_iss['sample'], df_r_iss['sample'])
        df_r_iss['sample'] = np.where(df_r_iss['study'] == 'UAMS', 'UAMS_' + df_r_iss['sample'], df_r_iss['sample'])
        df = pd.merge(df, df_r_iss[['sample', 'R_ISS', 'R2_ISS']], on='sample', how='inner')

    if args.ext_val:
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

    y_df = x_df[['pfs_code', 'pfs_time', 'os_code', 'os_time']]

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

    # write patients to predict
    x_df['sample'].to_csv(args.path + '/preds_treatment/patients.txt', index=False, header=False)

    # print('unique', df_x_entire.groupby('study').size())
    x_df = x_df.drop(['DARA', 'ELO'], axis=1)
    #x_df = x_df.drop(
    #    ['os_code', 'os_time', 'pfs_code', 'pfs_time', 'time_SCT', 'phases', 'stages', 'sub_stages', 'duration', 'study',
    #     'sample'], axis=1)

    df_x_entire = x_df.copy()

    return df_x_entire


def load_dataset():
    if args.ext_val:
        df = pd.read_csv(args.in_data_file + '/heidelberg_matrix.txt', sep='\t')
    else:
        df = pd.read_csv(args.in_data_file + '/PMMM_matrix_12052022.txt', sep='\t')

    df = df.drop('SCT_line', axis=1)
    df = df[df['duration'] >= 0]
    df = data_filter(df)

    return df


def plot_density(df_density, ax):
    sn.kdeplot(data=df_density, x="Risk", fill=True, ax=ax)
    ax.set_xlabel('Risk', fontsize=12)
    ax.xaxis.set_tick_params(labelsize=10)
    return ax


"""
def plot_state_risks(df_entire, single_df, ax):
    df = df_entire[df_entire['sample'] == single_df['sample'].unique()[0]].reset_index(drop=True)

    if df.size != 0:
        single_df = single_df.head(args.preds_time_in_days)
        single_df = single_df.drop(['sample', 'Time (Days)'], axis=1).abs()
        single_df.index = single_df.index / 365

        times_df = pd.DataFrame({'induction': df['duration'].values[0],
                                 'pfs_time': (df['pfs_time'].values[0] - df['duration'].values[0]),
                                 'os_time': df['os_time'].values[0] - df['pfs_time'].values[0]}, index=[0])
        times_df = times_df / 365

        single_df = single_df[['Alive in Induction', 'Alive in POD (induction)', 'Death in POD (induction)',
                               'Alive in Phase 2', 'Alive in POD', 'Death in Phase 2',
                               'Death in POD (Phase 2)']]

        single_df.columns = ['Non-POD \nin induction', 'POD in induction and alive',
                             'POD in induction and death',
                             'Non-POD and alive', 'POD after induction and alive',
                             'Non-POD and death', 'POD after induction and death']

        single_df.plot.area(stacked=True, ax=ax,
                            color=(
                                'lightgrey', 'lightseagreen', 'beige', 'pink', 'lightblue', 'gold', 'lightcoral'))

        ax2 = ax.twinx()

        # pfs
        c = ''
        if (df['pfs_time'].values[0] <= df['duration'].values[0]) and (df['pfs_code'].values[0] == 1):
            c = 'black'
            plt.scatter(df['pfs_time'].values[0] / 365, 0.07, color=c, s=60)
            print(df['pfs_code'].values[0], df['pfs_time'].values[0] / 365)
        elif (df['duration'].values <= df['pfs_time'].values) and (df['pfs_code'].values[0] == 1):
            c = 'blue'
            plt.scatter(df['pfs_time'].values[0] / 365, 0.07, color=c, s=60)

        # os
        c = ''
        if (df['os_time'].values[0] <= df['duration'].values[0]) and (df['os_code'].values[0] == 1):
            c = 'beige'
            plt.scatter(df['os_time'].values[0] / 365, 0.07, color=c, marker="x", s=50)

        elif (df['duration'].values[0] < df['os_time'].values[0] <= df['pfs_time'].values[0]) and (
                df['os_code'].values[0] == 1):
            c = 'goldenrod'
            plt.scatter(df['os_time'].values[0] / 365, 0.07, color=c, marker="x", s=50)

        elif (df['os_time'].values[0] > df['pfs_time'].values[0]) and (df['os_code'].values[0] == 1):
            c = 'red'
            plt.scatter(df['os_time'].values[0] / 365, 0.07, color=c, marker="x", s=50)

        '''

        times_df.plot.barh(stacked=True, rot=0, width=0.07, legend=False, ax=ax2,
                           color=('black', 'deeppink', 'dodgerblue'))
        '''

        ax.set_xlim(0, 5)
        ax2.set_ylim(0, 1)
        # ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax2.set_yticklabels([])
        # ax2.set_xticklabels([])
        # ax.axis('off')
        ax.tick_params(axis='x', labelsize=14)
        ax2.axis('off')
        legend = ax.legend(loc='lower left')
        legend.remove()
        return ax


def plot_stacked_indiv(df):
    '''
    Sediment plot of individual patients (args.sample)
    '''
    fig, ax = plt.subplots(args.top_display, 5, figsize=(16, 30), sharex=True)
    df = df[df['sample'] == args.sample].reset_index(drop=True)

    duration_time = df['duration'].values[0]
    os_time = df['os_time'].values[0]
    pfs_time = df['pfs_time'].values[0]

    # print(duration_time, pfs_time, os_time, df['os_code'].values[0])

    count_col = 0
    for ind_phase, phase in enumerate(['P2_SCT_only', 'P2_cont_only', 'P2_SCT_cont', 'P2_no_SCT_no_cont']):
        output_files_state_preds = glob.glob(args.path + '/preds_treatment/state_probs/' + args.sample +
                                       '~' + phase + '*' + args.surv_model + '.csv')

        # re-rank from least to highest risk (Progress p2)
        sorted_list =[]
        for ind, item in enumerate(output_files_state_preds):
            df_probs = pd.read_csv(output_files_state_preds[ind], index_col=None, header=0)
            sorted_list.append((output_files_state_preds[ind],
                                df_probs.iloc[duration_time - 1][['sample', 'Alive in POD (induction)']].values[1],
                                df_probs.iloc[args.pred_risk_at - 1][['sample', 'Alive in Phase 2']].values[1]
                                ))

        # sort by P1, P2
        if phase == 'P1':
            sorted_list = sorted(sorted_list, key=lambda x: x[1], reverse=False)
        else:
            sorted_list = sorted(sorted_list, key=lambda x: x[2], reverse=True)

        output_files_probs = [item[0] for item in sorted_list]

        # top 8
        output_files_probs = output_files_probs[:args.top_display]
        for ind, filename in enumerate(output_files_probs):
            ax_ = ax[ind, count_col]
            file_ = output_files_probs[ind].split('/')[-1]
            # df_probs = pd.read_csv(output_files_probs[ind], index_col=None, header=0)
            df_state_probs = pd.read_csv(args.path + '/loo/state_probs/' + file_, index_col=None, header=0)

            df_state_probs = df_state_probs.head(args.preds_time_in_days)

            plot_state_risks(df, df_state_probs, ax_)

            file_ = re.sub(',SCT_first_line|,continuos_treat', '', file_.split('~')[1])

            if phase == 'P1':
                # risk to prog in p1
                risk_ = round(df_state_probs.iloc[duration_time][['sample', 'Alive in POD (induction)']].values[1] * 100, 2)
                file_ = re.sub(phase + ',', '', file_)
            else:
                # risk of either progr. or dec. = 1 - prob. of alive in remission
                risk_ = round(100 - df_state_probs.iloc[args.pred_risk_at - 1][['sample', 'Alive in Phase 2']].values[1] * 100, 2)
                file_ = re.sub(phase + ',', '', file_)

            if ind == 0:
                ax_.set_title(phase + '\n' + file_ + '\nRisk is ' + str(risk_) + '%', fontsize=18, fontweight="bold")
            else:
                ax_.set_title(file_ + '\nRisk is ' + str(risk_) + '%', fontsize=18)

        count_col += 1

    plt.tight_layout()
    fig.savefig(args.path + '/plots/' + args.sample + '~' + args.surv_model + '.png')

    return None
"""


def plot_indiv_preds(sample_list, df_entire):
    """
    sample list provided
    """
    if args.ext_val:
        path_indiv = args.path + '/kfold/state_probs/all_states/'
        output_files = glob.glob(path_indiv + '*~hd~*fold~0~*' + args.surv_model + '.csv') # *~hd~*fold~0~* + args.surv_model
    else:
        path_indiv = args.path + '/loo/state_probs/'
        output_files = glob.glob(path_indiv + '*orig*neural_cox_non_prop.csv')

    df_ = pd.DataFrame(output_files)

    if args.ext_val:
        df_[1] = df_[0].apply(lambda x : x.split('/')[-1].split('~')[2])
    else:
        df_[1] = df_[0].apply(lambda x: x.split('/')[-1].split('~')[0])
    # df_ = df_.drop_duplicates([1])
    df_[1] = df_[1].apply(lambda x: re.sub('\.', '-', x))

    df_ = df_[df_[1].isin(sample_list)]
    df_ = df_.set_index(1).reindex(sample_list)
    samples_leave_one_out = df_[0].tolist()

    if args.plot_treat_combo_indiv_pat and args.ext_val:
        n_samples_rows = 16  ; n_samples_cols = 16
        fig, axes = plt.subplots(n_samples_rows, n_samples_cols, figsize=(36, 30), sharex=True)
    elif args.plot_treat_combo_indiv_pat and not args.ext_val:
        n_samples_rows = 4  ; n_samples_cols = 4
        fig, axes = plt.subplots(n_samples_rows, n_samples_cols, figsize=(16, 20), sharex=True)

    ax_list = [plt.subplot(n_samples_rows, n_samples_cols, i + 1) for i in range(len(samples_leave_one_out))]
    plt.subplots_adjust(wspace=0.06, hspace=0)

    for ind, (filename, ax) in enumerate(zip(samples_leave_one_out, ax_list)):
        if args.ext_val:
            sample = filename.split('/')[-1].split('~')[2]
            sample = re.sub('\.', '-', sample)
        else:
            sample = filename.split('/')[-1].split('~')[0]
        single_df = pd.read_csv(filename, index_col=None, header=0)
        single_df = single_df.head(args.preds_time_in_days)

        if single_df.size == 0:
            print('sample df', sample, single_df.shape)

        if args.ext_val:
            df = df_entire[df_entire['sample'] == sample].reset_index(drop=True)
        else:
            df = df_entire[df_entire['sample'] == single_df['sample'].unique()[0]].reset_index(drop=True)

        if df.size == 0:
            print('df', sample, df.shape)

        if df.size != 0:
            single_df = single_df.head(args.preds_time_in_days)
            if not args.ext_val:
                single_df = single_df.drop(['sample', 'Time (Days)'], axis=1).abs()
            else:
                single_df = single_df.drop(['Time (Days)'], axis=1).abs()
            single_df.index = single_df.index/ 365

            times_df = pd.DataFrame({'induction': df['duration'].values[0],
                                     'pfs_time': (df['pfs_time'].values[0] - df['duration'].values[0]),
                                     'os_time': df['os_time'].values[0] - df['pfs_time'].values[0]}, index=[0])
            times_df = times_df/365

            single_df = single_df[['Alive in Induction', 'Alive in POD (induction)', 'Death in POD (induction)',
                                'Alive in Phase 2', 'Alive in POD', 'Death in Phase 2', 'Death in POD (Phase 2)']]

            single_df.columns = ['Non-POD \nin induction', 'POD in induction and alive', 'POD in induction and death',
                     'Non-POD and alive', 'POD after induction and alive',
                     'Non-POD and death', 'POD after induction and death']

            single_df.plot.area(stacked=True, ax=ax,
                         color=('lightgrey', 'lightseagreen', 'beige', 'pink', 'lightblue', 'gold', 'lightcoral'))

            ax2 = ax.twinx()

            # efs
            c = ''
            if (df['pfs_time'].values <= df['duration'].values) and (df['pfs_code'].values[0] == 1) :
                c = 'lightgrey'
                plt.scatter(df['pfs_time'].values[0] / 365, 0.07, color=c, s=60)
            elif (df['pfs_time'].values >= df['os_time'].values) and (df['pfs_code'].values[0] == 1) :
                c = 'blue'
                plt.scatter(df['pfs_time'].values[0]/365, 0.07, color=c, s=60)
            elif (df['pfs_time'].values < df['os_time'].values) and (df['pfs_code'].values[0] == 1) :
                c = 'magenta'
                plt.scatter(df['pfs_time'].values[0]/365, 0.07, color=c, s=60)

            # os
            c = ''
            if (df['os_time'].values <= df['duration'].values) and (df['os_code'].values[0] == 1) :
                c = 'beige'
                plt.scatter(df['os_time'].values[0]/365, 0.07, color=c, marker="x", s=120)

            elif (df['duration'].values < df['os_time'].values <= df['pfs_time'].values) and (df['pfs_code'].values[0] == 0) and (df['os_code'].values[0] == 1):
                c = 'cyan'
                plt.scatter(df['os_time'].values[0]/365, 0.07, color=c, marker="x", s=120)

            elif (df['os_time'].values > df['pfs_time'].values) and (df['pfs_code'].values[0] == 1) and (df['os_code'].values[0] == 1):
                c = 'red'
                plt.scatter(df['os_time'].values[0]/365, 0.07, color=c, marker="x", s=120)

            times_df.plot.barh(stacked=True, rot=0, width=0.07, legend=False, ax=ax2, color=('black', 'deeppink', 'dodgerblue'))

            # ax.set_title(sample)
            ax.set_xlim(0, 5)
            ax2.set_ylim(0,1)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax2.set_yticklabels([])
            ax2.set_xticklabels([])
            ax.axis('off')
            ax2.axis('off')
            legend = ax.legend(loc='lower left')
            legend.remove()
            print('done')

    if args.plot_treat_combo_indiv_pat:
        fig.savefig(args.path + '/plots/pats_preds_ext_val_' + str(args.ext_val) + '_' + args.surv_model + '_at_' + str(args.preds_time_in_days) + '_days.png')

    return None


def plot_treat_cohort(df):
    """
    Treatment plot of individual patients
    """


    df_copy = df.copy()
    list_of_lists_all_pats = []
    '''
    # run only once ; computationally expensive
    for ind_phase, phase in enumerate(['P2_SCT_only', 'P2_cont_only', 'P2_SCT_cont', 'P2_no_SCT_no_cont']):
        if args.ext_val:
            output_files_state_preds = glob.glob(args.path + '/kfold/state_probs/all_states/' + '*~hd~*fold~0~*' + phase + '*' + args.surv_model + '.csv')
        else:
            output_files_state_preds = glob.glob(args.path + '/loo/state_probs/' +
                                       '*~' + phase + '*' + args.surv_model + '.csv')

        for filename in output_files_state_preds:
            df_state_probs = pd.read_csv(filename, index_col=None, header=0)
            df_state_probs = df_state_probs.head(args.preds_time_in_days)

            sample_id = filename.split('/')[-1].split('~')[0]
            filename_ = re.sub(',orig', '', filename.split('/')[-1].split('~')[1])
            state_ = filename.split('/')[-1].split('~')[1].split(',')[0]
            treat_ = re.sub(state_ + ',', '', filename_)

            df = df_copy[df_copy['sample'] == sample_id.split('/')[-1].split('~')[0]].reset_index(drop=True)
            duration_time = df['duration'].values[0]
            os_time = df['os_time'].values[0]
            pfs_time = df['pfs_time'].values[0]

            if phase == 'P1':
                # risk to prog in p1
                risk_pfs = round(df_state_probs.iloc[duration_time][['sample', 'Alive in POD (induction)']].values[1] * 100, 2)
                risk_os = round(df_state_probs.iloc[duration_time][['sample', 'Alive in Induction']].values[1] * 100, 2)
            else:
                # risk of either progr. or dec. = 1 - prob. of alive in remission
                risk_pfs = round(100 - df_state_probs.iloc[args.pred_risk_at-1][['sample', 'Alive in Phase 2']].values[1] * 100, 2)
                risk_os = round(
                    100 -
                    (df_state_probs.iloc[args.pred_risk_at-1][['sample', 'Alive in Phase 2']].values[1] * 100
                    + df_state_probs.iloc[args.pred_risk_at-1][['sample', 'Death in Phase 2']].values[1] * 100),
                    2)

            list_of_lists_all_pats.append([sample_id, state_, treat_,
                                           duration_time, pfs_time, df['pfs_code'].values[0],
                                           os_time, df['os_code'].values[0], risk_pfs, risk_os])

    df_table = pd.DataFrame(list_of_lists_all_pats,
                            columns=['sample_id', 'state', 'treatment', 'duration', 'pfs_time', 'pfs_code', 'os_time',
                                                'os_code', 'risk_pfs', 'risk_os'])

    df_table.to_csv(args.path + '/preds_treatment/df_table_risk_ext_val' + str(args.ext_val) + '.csv', index=False)
    '''
    ### comment till here

    # comment till here

    # view 1

    df_table = pd.read_csv(args.path +'/preds_treatment/df_table_risk_ext_val' + str(args.ext_val) + '.csv')

    # conditions
    # df_table = df_table[df_table['state'] != 'P1']
    # df_table = df_table[(df_table['os_time'] <= args.pred_risk_at-1) & (df_table['os_time'] >= 365)]

    df_treat_combo = pd.read_csv(args.path + '/preds_treatment/combo_freq.csv')
    df_treat_combo = df_treat_combo[df_treat_combo['Freq'] > 10]

    df_table = df_table[df_table['treatment'].isin(df_treat_combo['Combo'].tolist())]

    # view 2
    treat_df = pd.pivot_table(df_table, index=['sample_id', 'treatment'],
                                       columns=['state'], values=[args.outcome]).reset_index()
    treat_df.columns = [''.join(col).strip() for col in treat_df.columns.values]
    treat_df.columns = [re.sub(args.outcome, '', item) for item in treat_df.columns]

    ### reorder
    state_order = ["P2_no_SCT_no_cont", "P2_cont_only", "P2_SCT_only", "P2_SCT_cont"]
    df_table = df_table.set_index('state').loc[state_order].reset_index()

    df_table['state'] = np.select([(df_table['state'] == 'P2_SCT_only'),
                                   (df_table['state'] == 'P2_cont_only'),
                                   (df_table['state'] == 'P2_no_SCT_no_cont'),
                                   (df_table['state'] == 'P2_SCT_cont')],
                                  ['SCT', 'Cont. treatment', 'No SCT No cont. treatment', 'SCT & cont. treatment'])

    ## variance
    df_table_withoutP1 = df_table.copy()
    df_for_var = pd.pivot_table(df_table_withoutP1, index=['sample_id', 'state'],
                              columns=['treatment'], values=[args.outcome]).reset_index()
    df_for_var.columns = [''.join(col).strip() for col in df_for_var.columns.values]
    df_for_var.columns = [re.sub(args.outcome, '', item) for item in df_for_var.columns]

    ## for variance plots
    sct = df_for_var[(df_for_var['state'] == 'SCT')].drop('state', axis=1).set_index('sample_id')
    sct.columns = ['SCT-' + str(col) for col in sct.columns]

    cont_treat = df_for_var[(df_for_var['state'] == 'Cont. treatment')].drop('state', axis=1).set_index('sample_id')
    cont_treat.columns = ['Cont. treatment-' + str(col) for col in cont_treat.columns]

    no_sct_no_cont = df_for_var[(df_for_var['state'] == 'No SCT No cont. treatment')].drop('state',axis=1).set_index('sample_id')
    no_sct_no_cont.columns = ['No SCT No cont. treatment-' + str(col) for col in no_sct_no_cont.columns]

    sct_cont = df_for_var[(df_for_var['state'] == 'SCT & cont. treatment')].drop('state',axis=1).set_index('sample_id')
    sct_cont.columns = ['SCT & cont. treatment-' + str(col) for col in sct_cont.columns]

    dfs = [sct, cont_treat, sct_cont, no_sct_no_cont]
    df_for_var_matrix = reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True), dfs)

    df_for_var_matrix.to_csv(args.path + '/preds_treatment/df_treatment_ext_val_' + str(args.ext_val) + '_' + args.outcome + '.csv')
    df = df.sort_values(by='sample')
    df_for_var_matrix = df_for_var_matrix.sort_values(by='sample_id')

    # std, var across all treatment combinations
    var_mat_treats = df_for_var_matrix.T.var()
    var_mat_treats = pd.DataFrame(var_mat_treats, index=var_mat_treats.index).rename(columns={0: 'variance'})

    # for clust
    var_mat_treats = pd.merge(var_mat_treats, df_for_var_matrix.min(axis=1).rename('best_risk'),
                                           left_index=True, right_index=True, how='inner')
    var_mat_treats = pd.merge(var_mat_treats, df_for_var_matrix.max(axis=1).rename('worst_risk'),
                              left_index=True, right_index=True, how='inner')

    var_mat_treats.to_csv(args.path + '/preds_treatment/df_var_ext_val_' + str(args.ext_val) + '_' + args.outcome + '.csv')

    '''
    std_clf = make_pipeline(StandardScaler())
    scaler = std_clf.named_steps["standardscaler"]
    normalized_for_clust = scaler.fit_transform(var_mat_treats.to_numpy())

    num_clust = 11
    temp_predict_labels = ''
    in_ =  normalized_for_clust
    score_list = []
    for i in range(2, 33, 1):
        gm = GaussianMixture(n_components=i, random_state=688).fit(in_)
        predict_labels = gm.predict(in_)
        score = metrics.silhouette_score(in_, predict_labels, metric='euclidean')
        # sample_silhouette_values = metrics.silhouette_samples(var_mat_treats.to_numpy(), predict_labels)
        score_list.append((i, score)) # gm.means_,
        if i == num_clust:
            temp_predict_labels = predict_labels
        print('done')

    fig = plt.figure(figsize=(13,15))
    ax = plt.subplot(1, 3, 1)
    n_components = np.arange(1, 33)
    models = [GaussianMixture(n, covariance_type='full', random_state=0).fit(in_) for n in n_components]
    ax.plot(n_components, [m.bic(in_) for m in models], label='BIC')
    ax.plot(n_components, [m.aic(in_) for m in models], label='AIC')
    ax.legend(loc='best')
    ax.set_xlabel('Clusters')
    ax.set_title('A. AIC, BIC', fontsize=14)

    ax = plt.subplot(1, 3, 2)
    df_sil_score = pd.DataFrame(score_list, columns=['clusters', 'silhouette_score'])
    df_sil_score.sort_values(by='silhouette_score').plot(kind='barh', x='clusters', y='silhouette_score', ax=ax)
    ax.set_title('B. Silhouette score', fontsize=14)

    temp_predict_labels_cp = temp_predict_labels.copy()
    to_save = pd.DataFrame(temp_predict_labels_cp, index=var_mat_treats.index).rename(columns={0: str(num_clust) + '_clust'})
    # to_save.to_csv(args.path + '/treat_preds/' + str(num_clust) + '_clust.csv')

    ax = plt.subplot(2, 3, 3)
    temp_predict_labels = pd.DataFrame(temp_predict_labels, columns=['cluster_num'])\
        .reset_index(drop=True).value_counts().reset_index().rename(columns={0: 'count_patients'})
    temp_predict_labels.sort_values(by='count_patients').plot(kind='barh', x='cluster_num', y='count_patients', ax=ax)
    ax.set_title('C. Cluster distribution (n=' + str(num_clust) + ')', fontsize=14)
    plt.tight_layout()
    fig.savefig(args.path + '/plots/' + str(num_clust) + '_cluster_distribution_' + args.outcome + '.png')
    plt.close()

    # 2d tsne
    per_ = 20
    tsne = TSNE(n_components=3, init = 'random', perplexity=per_, random_state=20)
    var_transformed_tsne = tsne.fit_transform(var_mat_treats.to_numpy())

    fig, ax = plt.subplots()
    sn.scatterplot(var_transformed_tsne[:, 0], var_transformed_tsne[:, 1],
                   hue=temp_predict_labels_cp, legend='full', ax=ax, palette="tab20")
    ax.set_title('Tsne', fontsize=14)
    plt.tight_layout()
    fig.savefig(args.path + '/plots/' + str(num_clust) + '_tsne_2d_' + args.outcome + '.png')
    plt.close()

    # 3d interactive tsne
    tsne = TSNE(n_components=3, init = 'random', perplexity=per_, random_state=20)
    var_transformed_tsne = tsne.fit_transform(var_mat_treats.to_numpy())

    temp_predict_labels_cp_3d = pd.DataFrame(temp_predict_labels_cp, index=var_mat_treats.index).rename(columns={0: 'clust'})
    temp_predict_labels_cp_3d = pd.merge(temp_predict_labels_cp_3d, var_mat_treats, on='sample_id')
    temp_predict_labels_cp_3d = pd.merge(temp_predict_labels_cp_3d, df, left_on='sample_id', right_on='sample')

    temp_predict_labels_cp_3d = pd.merge(temp_predict_labels_cp_3d,
                                         pd.DataFrame(var_transformed_tsne), how='inner',
                                         left_index=True, right_index=True)

    temp_predict_labels_cp_3d['variance'] = temp_predict_labels_cp_3d.variance.round(2)
    temp_predict_labels_cp_3d = temp_predict_labels_cp_3d.fillna('NA')

    fig_3d = px.scatter_3d(data_frame=temp_predict_labels_cp_3d, x=0, y=1, z=2,
                           custom_data=['ISS', 'age', 'variance', 'gender', 'SCT_first_line', 'continuos_treat',
                                        ],
                           color=temp_predict_labels_cp_3d.clust.astype(str), labels={'color': 'clust'},
                           title="Tsne to map therapies", opacity=0.5
            )

    fig_3d.update_traces(
        hovertemplate="<br>".join([
            # "0: %{x}",
            # "1: %{y}",
            # "2: %{z}",
            "ISS: %{customdata[0]}",
            "Age: %{customdata[1]}",
            "Variance of treatment: %{customdata[2]}",
            "Gender: %{customdata[3]}",
            "SCT: %{customdata[4]}",
            "Cont. treat: %{customdata[5]}",
            #"CNV_1q23.3_AMP: %{customdata[6]}"
        ])
    )

    fig_3d.show()
    '''

    # boxplot spread stats
    box_plot_stats = mpl.cbook.boxplot_stats(var_mat_treats)
    box_plot_stats = box_plot_stats[0]
    box_plot_stats['fliers'].sort()

    range_of_box_plot_qar = [(box_plot_stats['whislo'], box_plot_stats['q1']),
                             (box_plot_stats['q1'], box_plot_stats['q3']),
                             (box_plot_stats['q3'], box_plot_stats['whishi']),
                             (box_plot_stats['whishi'], box_plot_stats['fliers'][-1])]

    # 2. plot to get pats where treatment matters, dont or slightly
    var_mat_treats_cp = var_mat_treats.copy()
    n_choose_samples = 4 # from each quartile of boxplot
    samples_to_plot = []
    for ind, (low_val, high_val) in enumerate(range_of_box_plot_qar):
        df_sub = var_mat_treats_cp[(var_mat_treats_cp['variance'] >= low_val) & (var_mat_treats_cp['variance'] < high_val)]
        df_sub = df_sub.sort_values(by='variance', ascending=True)
        samples_ = df_sub.index.tolist()
        if ind == 0:
            samples_ = samples_[:n_choose_samples]
        elif ind == 1:
            samples_ = samples_[int(len(samples_)/4):int(len(samples_)/2)][:n_choose_samples]
        elif ind == 2 or ind == 3:
            samples_ = samples_[-n_choose_samples:]

        samples_to_plot.extend(samples_)

    list_of_cases = []
    fig, axes = plt.subplots(n_choose_samples,n_choose_samples, figsize=(14, 14), sharex=True)

    for ind, sampl in enumerate(samples_to_plot):
        if 0 <= ind <= 3:
            j = 0
            i = ind
        elif 4 <= ind <= 7:
            j = 1
            i = ind - 4
        elif 8 <= ind <= 11:
            j = 2
            i = ind - 8
        elif 12 <= ind <= 15:
            j = 3
            i = ind - 12

        ax = axes[j, i]

        indiv = df_table_withoutP1[df_table_withoutP1['sample_id'] == sampl]
        sn.pointplot(x="state", y=args.outcome, hue='treatment', palette="Paired", data=indiv, ax=ax)
        list_of_cases.append(sampl)
        # ax.set_title(sampl)
        ax.set_ylim(20, 100)
        ax.tick_params(axis='x', labelsize=14, rotation=75)
        ax.tick_params(axis='y', labelsize=14)
        plt.xticks(fontname="Arial")
        plt.yticks(fontname="Arial")
        if ind != 15:
            legend = ax.legend(loc='lower left')
            legend.remove()
        else:
            legend = ax.legend(loc='center left', bbox_to_anchor=(1, 2), prop={'size': 12, 'family': 'Arial'})
            legend.get_frame().set_alpha(0.1)

        ax.set_xlabel('')
        if ind % 4 == 0:
            ax.set_ylabel('Risk', fontsize=14, fontname='Arial')
        else:
            ax.set_ylabel('')

    with open(args.path + '/preds_treatment/indiv_cases_treat_preds.txt', 'w') as f:
        for line in list_of_cases:
            f.write(f"{line}\n")

    fig.savefig(args.path + '/plots/indiv_treat_preds_ext_val_' + str(args.ext_val) + '_' + args.outcome + '.png', bbox_inches='tight')

    '''
    ## plot treatment variance
    fig, ax = plt.subplots(figsize=(2,4))
    var_mat_treats = pd.DataFrame(var_mat_treats.variance).boxplot(ax=ax)
    ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.grid(visible=True, which='major', axis='both', color='w', linewidth=1.0)
    ax.grid(visible=True, which='minor', axis='both', color='w', linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_color('black')
    ax.spines['left'].set_linewidth(1)

    plt.tight_layout()
    fig.savefig(args.path + '/plots/var_mat_ext_val_' + str(args.ext_val) + '_' + args.outcome + '.png')
    '''

    return None


def barplot_combo_freq(list_of_freq):
    fig,ax = plt.subplots(figsize=(6, 3))
    df_ = pd.DataFrame.from_records(list_of_freq, columns=['Combo', 'Freq'])
    df_ = df_.drop_duplicates()
    print('occurrences', df_['Freq'].sum())

    df_.to_csv(args.path + '/preds_treatment/combo_freq.csv', index=False)

    df_[::-1].plot.barh(x='Combo', y='Freq', legend=False, ax=ax)

    plt.xlabel('Frequency', fontname='Arial', fontsize=10)
    plt.ylabel('Therapy', fontname='Arial', fontsize=10)
    ax.tick_params(axis='both', which='major', labelsize=8)
    plt.tight_layout()
    fig.savefig(args.path + '/plots/combo_freq.png')

    return None


def get_treat_combo_stats(df):
    # df = df[df['phases'] == 'phase_2']
    set_treat_two = [",".join(map(str, comb)) for comb in
                     combinations(df[['chemo', 'KAR', 'LEN', 'BORT', 'THAL', ]], 2)]

    set_treat_three = [",".join(map(str, comb)) for comb in
                     combinations(df[['chemo', 'KAR', 'LEN', 'BORT', 'THAL', ]], 3)]

    set_treat_quad = [",".join(map(str, comb)) for comb in
                     combinations(df[['chemo', 'KAR', 'LEN', 'BORT', 'THAL', ]], 4)]

    df_sub = df.copy()
    list_all_treats = ['chemo', 'KAR', 'LEN', 'BORT', 'THAL', 'SCT_first_line']
    ranked_combos = []
    for item in set_treat_two + set_treat_three + set_treat_quad:
        combo_id = item.split(',')
        if len(combo_id) == 2:
            df_sub = df.loc[(df[combo_id[0]] == 1) & (df[combo_id[1]] == 1)]
        if len(combo_id) == 3:
            df_sub = df.loc[(df[combo_id[0]] == 1) & (df[combo_id[1]] == 1) & (df[combo_id[2]] == 1)]
        if len(combo_id) == 4:
            df_sub = df.loc[(df[combo_id[0]] == 1) & (df[combo_id[1]] == 1) & (df[combo_id[2]] == 1)
                                & (df[combo_id[3]] == 1)]

        rest = list(set(list_all_treats) - set(combo_id))
        df_sub = df_sub.loc[(df_sub[rest] == 0).all(axis=1)]

        if df_sub.shape[0] != 0:
            # print(item, df_sub.shape)
            ranked_combos.append((item, df_sub.shape[0]))

    ranked_combos = sorted(ranked_combos, key=lambda x: x[1], reverse=True)
    barplot_combo_freq(ranked_combos)

    return None


def main():
    global args
    args = parser.parse_args()

    df = load_dataset()

    if args.get_treat_stats:
        get_treat_combo_stats(df)

    if args.plot_treat_combo_indiv_pat:
        # plot_stacked_indiv(df)
        plot_indiv_preds(df['sample'].tolist(), df)

    if args.plot_treat_combo_all_pats:
        plot_treat_cohort(df)

    return None


if __name__ == '__main__':
    main()