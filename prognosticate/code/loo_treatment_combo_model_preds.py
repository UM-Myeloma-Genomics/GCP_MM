import pandas as pd
import argparse
import warnings
import torch
import torchtuples as tt
from pycox.models import CoxTime
from pycox.models.cox_time import MLPVanillaCoxTime
from sklearn_pandas import DataFrameMapper
import re
from modules_preprocess import imputer_knn
from modules_mstate import *

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description='Leave one out cross validation based on different therapy strategies')
parser.add_argument('--seed', type=int, default=123, help='random seed (default: 42)')
parser.add_argument('--surv_model', type=str, default='neural_cox_non_prop', help='model to train to survival model | rsf, cph, svms')
parser.add_argument('--max_time_days', type=int, default=1825, help='predict curves till time t (days)')
parser.add_argument('--in_data_file', type=str, default='/Users/axr2376/Desktop/pred_1_0_paper/data_in/', help='input dataset path')
parser.add_argument('--path', type=str, default='/Users/axr2376/Desktop/pred_1_0_paper/data_out/expt_1', help='input dataset path')
parser.add_argument('--thresh_m_to_p2', type=int, default=365, help='days')
parser.add_argument('--sample', type=str, default='PD5886a', help='SAMPLE ids')
parser.add_argument('--legend', type=str, default='/Users/axr2376/Desktop/pred_1_0_paper/data_in/legend_PMMM_2022.xlsx', help='input dataset path')
parser.add_argument('--get_feats_group', type=str, default='9', help='if set None (100) then all feats')


def data_filter(df):
    df = df[df['duration'] >= 0]
    df = df.dropna(subset=['os_time', 'pfs_time', 'pfs_code', 'os_code'])

    df['time_SCT'].fillna(0, inplace=True)

    df_r_iss = pd.read_csv(args.in_data_file + '/ISS_RISS_R2ISS_LDH_all_cohort.txt', sep='\t')
    df_r_iss['sample'] = np.where(df_r_iss['study'] == 'MRC_XI', 'MRC_' + df_r_iss['sample'], df_r_iss['sample'])
    df_r_iss['sample'] = np.where(df_r_iss['study'] == 'UAMS', 'UAMS_' + df_r_iss['sample'], df_r_iss['sample'])
    df = pd.merge(df, df_r_iss[['sample', 'R_ISS', 'R2_ISS']], on='sample', how='inner')

    # map clinical vars
    map_gender = {'male': 0, 'female': 1}
    x_df = df.replace({'gender': map_gender})
    map_study = {'MMRF': 0, 'MPG': 1, 'Moffit': 2, 'MSKCC_292': 3, 'UAMS': 4}
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

    df_x_entire = x_df.copy()

    # print('unique', df_x_entire.groupby('study').size())

    x_df = x_df.drop(
        ['os_code', 'os_time', 'pfs_code', 'pfs_time', 'time_SCT', 'phases', 'stages', 'sub_stages', 'duration', 'study',
         'sample',  'DARA', 'ELO'], axis=1)

    # select groups of features for c-index
    x_df = select_features(x_df)
    # print('dataset shape is', x_df.shape)

    return x_df, y_df, df_x_entire


def select_features(df):
    df_fea_anno = pd.read_excel(args.legend)
    if args.get_feats_group == '10':
        df_fea_anno = df_fea_anno[
                            (df_fea_anno['column_id'] != 'time_SCT') & (df_fea_anno['column_id'] != 'study') &
                            (df_fea_anno['column_id'] != 'sample') & (df_fea_anno['column_id'] != 'phase') &
                            (df_fea_anno['column_id'] != 'os_time') & (df_fea_anno['column_id'] != 'pfs_time') &
                            (df_fea_anno['column_id'] != 'duration') & (df_fea_anno['column_id'] != 'pfs_code') &
                            (df_fea_anno['column_id'] != 'os_code') & (df_fea_anno['column_id'] != 'SCT_line')
                            & (df_fea_anno['column_id'] != 'DARA') & (df_fea_anno['column_id'] != 'ELO')
                            ]
    else:
        df_fea_anno = df_fea_anno[
            (df_fea_anno['column_id'] != 'time_SCT') & (df_fea_anno['column_id'] != 'study') &
            (df_fea_anno['column_id'] != 'sample') & (df_fea_anno['column_id'] != 'phase') &
            (df_fea_anno['column_id'] != 'os_time') & (df_fea_anno['column_id'] != 'pfs_time') &
            (df_fea_anno['column_id'] != 'duration') & (df_fea_anno['column_id'] != 'pfs_code') &
            (df_fea_anno['column_id'] != 'os_code') & (df_fea_anno['column_id'] != 'SCT_line')
            & (df_fea_anno['column_id'] != 'DARA') & (df_fea_anno['column_id'] != 'ELO') & (df_fea_anno['column_id'] != 'RNA_GEP70_dichotomous')
            ]

    df_ = df[df_fea_anno['column_id'].tolist()]

    # select group of feats
    master_df = pd.read_csv(args.path + '/feat_matrix/main_keeper.csv')
    df_ = df_[master_df.loc[int(args.get_feats_group), 'feat_combo'].split(' ')]

    # df_ = df_.drop(['R_ISS', 'R2_ISS'], axis=1)

    return df_


def raw_data_read_prepare():
    df = pd.read_csv(args.in_data_file + '/PMMM_matrix_12052022.txt', sep='\t')

    df = df.drop('SCT_line', axis=1)

    df = df.sort_values(by=['sample'])

    x_df, y_df, df_x_entire = data_filter(df)

    return x_df, y_df, df_x_entire


def cox_time_orig(train_x, train_y, test_x, cols_x, model_param):
    np.random.seed(args.seed)
    _ = torch.manual_seed(model_param[2])

    cols_leave = [x for x in cols_x ]
    leave = [(col, None) for col in cols_leave]
    x_mapper = DataFrameMapper(leave)

    train_x = pd.DataFrame(train_x, columns=cols_x)
    val_x = train_x.sample(frac=0.2)
    test_x = pd.DataFrame(test_x, columns=cols_x)

    train_y = train_y[['binary_event', 'event_time']]

    val_y = train_y.loc[val_x.index, :]

    # delete
    train_x = train_x.drop(val_x.index)
    train_y = train_y.drop(val_y.index)

    train_x = x_mapper.fit_transform(train_x).astype('float32')
    val_x = x_mapper.fit_transform(val_x).astype('float32')
    test_x = x_mapper.transform(test_x).astype('float32')

    labtrans = CoxTime.label_transform()
    get_target = lambda df: (df['event_time'].values, df['binary_event'].values)
    train_y = labtrans.fit_transform(*get_target(train_y))

    y_val = labtrans.transform(*get_target(val_y))
    val = tt.tuplefy(val_x, y_val)

    in_features = train_x.shape[1]
    num_nodes = model_param[0]
    batch_norm = False
    dropout = model_param[4]
    net = MLPVanillaCoxTime(in_features, num_nodes, batch_norm, dropout)

    optimizer = tt.optim.AdamWR(decoupled_weight_decay=0.025, cycle_eta_multiplier=0.5)
    model = CoxTime(net, optimizer, labtrans=labtrans)

    model.optimizer.set_lr(0.005)
    batch_size = model_param[1]
    epochs = 20
    callbacks = [tt.cb.EarlyStopping()]

    log = model.fit(train_x, train_y, batch_size, epochs, callbacks, val_data=val, val_batch_size=batch_size, verbose=False)

    _ = log.plot()

    _ = model.compute_baseline_hazards()

    net.eval()
    with torch.set_grad_enabled(False):
        surv_test = model.predict_surv_df(test_x)

    return surv_test


def cox_time(train_x, train_y, test_x, cols_x, model_param):
    np.random.seed(args.seed)
    _ = torch.manual_seed(model_param[2])

    cols_leave = [x for x in cols_x ]
    leave = [(col, None) for col in cols_leave]
    x_mapper = DataFrameMapper(leave)


    train_x = pd.DataFrame(train_x, columns=cols_x)
    val_x = train_x.sample(frac=0.2)
    test_x = pd.DataFrame(test_x, columns=cols_x)

    train_y = train_y[['binary_event', 'event_time']]

    val_y = train_y.loc[val_x.index, :]

    # delete
    train_x = train_x.drop(val_x.index)
    train_y = train_y.drop(val_y.index)

    train_x = x_mapper.fit_transform(train_x).astype('float32')
    val_x = x_mapper.fit_transform(val_x).astype('float32')
    test_x = x_mapper.transform(test_x).astype('float32')

    labtrans = CoxTime.label_transform()
    get_target = lambda df: (df['event_time'].values, df['binary_event'].values)
    train_y = labtrans.fit_transform(*get_target(train_y))

    y_val = labtrans.transform(*get_target(val_y))
    val = tt.tuplefy(val_x, y_val)

    in_features = train_x.shape[1]
    num_nodes = model_param[0]
    batch_norm = False
    dropout = model_param[4]
    net = MLPVanillaCoxTime(in_features, num_nodes, batch_norm, dropout)

    optimizer = tt.optim.AdamWR(decoupled_weight_decay=0.025, cycle_eta_multiplier=0.5)
    model = CoxTime(net, optimizer, labtrans=labtrans)

    model.optimizer.set_lr(0.005)
    batch_size = model_param[1]
    epochs = 20
    callbacks = [tt.cb.EarlyStopping()]

    log = model.fit(train_x, train_y, batch_size, epochs, callbacks, val_data=val, val_batch_size=batch_size, verbose=False)

    _ = log.plot()

    _ = model.compute_baseline_hazards()

    # net.eval()
    #   with torch.set_grad_enabled(False):
    surv_test = model.predict_surv_df(test_x)

    return surv_test


def progress_p1(X, Y, df, sample, model_param_list):
    # print('Progress P1', Y['os_time'].min(), Y['os_time'].max())
    X, Y = crit_prog_p1(X, Y, df)

    df_pred = train_survival_(X.drop(['SCT_first_line', 'continuos_treat'], axis=1), Y,
                              sample.drop(['SCT_first_line', 'continuos_treat'], axis=1),
                              'Progress (P1)', model_param_list[0])

    # missing probabilities
    idx = np.arange(0, args.max_time_days, 1)
    df_interp_pred = pd.merge(df_pred, pd.DataFrame(index=idx),
                                       left_index=True, right_index=True, how='outer')
    df_interp_pred = df_interp_pred.interpolate()

    return df_interp_pred


def progress_dec_p1(X, Y, df, sample, model_param_list):
    # print('Progress dec P1', Y['os_time'].min(), Y['os_time'].max())

    X, Y = crit_prog_dec_p1(X, Y, df)

    df_pred = train_survival_(X.drop(['SCT_first_line', 'continuos_treat'], axis=1), Y,
                              sample.drop(['SCT_first_line', 'continuos_treat'], axis=1),
                          'Progress & deceased (P1)', model_param_list[2])

    # missing probabilities
    idx = np.arange(0, args.max_time_days, 1)
    df_interp_pred = pd.merge(df_pred, pd.DataFrame(index=idx),
                                       left_index=True, right_index=True, how='outer')
    df_interp_pred = df_interp_pred.interpolate()

    return df_interp_pred


def move_to_phase2(X, Y, df, sample, model_param_list):
    # print('Move to phase 2', Y['os_time'].min(), Y['os_time'].max())

    X, Y = crit_move_to_p2(X, Y, df)

    df_pred = train_survival_(X.drop(['SCT_first_line', 'continuos_treat'], axis=1), Y,
                              sample.drop(['SCT_first_line', 'continuos_treat'], axis=1),
                              'Move to phase 2', model_param_list[0])

    # missing probabilities
    idx = np.arange(0, args.max_time_days, 1)
    df_interp_pred = pd.merge(df_pred, pd.DataFrame(index=idx),
                                       left_index=True, right_index=True, how='outer')
    df_interp_pred = df_interp_pred.interpolate()
    df_interp_pred['Move to phase 2'].iloc[args.thresh_m_to_p2+1:] = 0
    return df_interp_pred


def progress_p2(X, Y, df, sample, model_param_list):
    group = 'Progress (P2)'
    # print(group, Y['os_time'].min(), Y['os_time'].max())
    X, Y = crit_prog_p2(X, Y, df)

    df_pred = train_survival_(X, Y, sample, group, model_param_list[2])

    # missing probabilities
    # print('missing values')
    idx = np.arange(0, args.max_time_days, 1)
    df_interp_pred = pd.merge(df_pred, pd.DataFrame(index=idx),
                                       left_index=True, right_index=True, how='outer')
    df_interp_pred = df_interp_pred.interpolate()
    # print(df_interp_pred)
    return df_interp_pred


def progress_dec_p2(X, Y, df, sample, model_param_list):
    group = 'Progress & deceased (P2)'
    # print(group, Y.shape, Y['os_time'].min(), Y['os_time'].max())

    X, Y = crit_prog_dec_p2(X, Y, df)

    df_pred = train_survival_(X, Y, sample, group, model_param_list[3])

    # missing probabilities
    idx = np.arange(0, args.max_time_days, 1)
    df_interp_pred = pd.merge(df_pred, pd.DataFrame(index=idx),
                                       left_index=True, right_index=True, how='outer')
    df_interp_pred = df_interp_pred.interpolate()
    return df_interp_pred


def non_progress_dec_p2(X, Y, df, sample, model_param_list):
    group = 'Non-progress & deceased (P2)'
    # print(group, Y['os_time'].min(), Y['os_time'].max())

    X, Y = crit_non_prog_dec_p2(X, Y, df)
    df_pred = train_survival_(X, Y, sample, group, model_param_list[4])

    # missing probabilities
    idx = np.arange(0, args.max_time_days, 1)
    df_interp_pred = pd.merge(df_pred, pd.DataFrame(index=idx),
                                       left_index=True, right_index=True, how='outer')
    df_interp_pred = df_interp_pred.interpolate()
    return df_interp_pred


def data_prep(sample):
    X_not_norm, Y, df_X_entire = raw_data_read_prepare()

    X_not_norm = X_not_norm.reset_index(drop=True)
    Y = Y.reset_index(drop=True)
    df_X_entire = df_X_entire.reset_index(drop=True)

    print(list(X_not_norm))
    # print('data prep', Y['os_time'].min(), Y['os_time'].max())

    cols = list(X_not_norm)

    # impute entire dataset
    X_not_norm = imputer_knn(X_not_norm)
    X_not_norm = pd.DataFrame(X_not_norm, columns=cols)

    # get test sample
    test_sample = df_X_entire[(df_X_entire['sample'] == sample)].reset_index()
    test_x = X_not_norm[X_not_norm.index.isin(test_sample['index'])]
    test_y = Y[Y.index.isin(test_sample['index'])]

    X_not_norm = X_not_norm.drop(index=test_sample['index'])
    Y = Y.drop(index=test_sample['index'])
    df_X_entire = df_X_entire.drop(index=test_sample['index'])

    # reset indices
    X_not_norm = X_not_norm.reset_index(drop=True)
    Y = Y.reset_index(drop=True)
    test_x = test_x.reset_index(drop=True)
    test_y = test_y.reset_index(drop=True)
    df_X_entire = df_X_entire.reset_index(drop=True)

    # print(test_x.shape, X_not_norm.shape, test_y.shape, Y.shape)

    return X_not_norm, Y, df_X_entire, test_x


def phase_1(X_not_norm, Y, df_X_entire, test_x, model_param_list):
    # print('phase 1 prog', Y['os_time'].min(), Y['os_time'].max())

    df_p1_prog = progress_p1(X_not_norm, Y, df_X_entire, test_x, model_param_list)
    # print('phase 1 prog and dec', Y['os_time'].min(), Y['os_time'].max())

    # G2 Phase 1 induction progressors alive vs deceased
    df_p1_prog_alive = progress_dec_p1(X_not_norm, Y, df_X_entire, test_x, model_param_list)

    # G2 phase 2 vs 1
    df_move_to_phase2 = move_to_phase2(X_not_norm, Y, df_X_entire, test_x, model_param_list)
    df_move_to_phase2.iloc[366:] = 0

    return df_p1_prog, df_p1_prog_alive, df_move_to_phase2


def phase_2(X_not_norm, Y, df_X_entire, test_x, model_param_list):
    df_p2 = df_X_entire[(df_X_entire['phases'] == 'phase_2')]
    X_p2 = X_not_norm[X_not_norm.index.isin(df_p2.index)]
    Y_p2 = Y[Y.index.isin(df_p2.index)]

    # reset
    post_induction = df_p2[(df_p2['phases'] == 'phase_2')]
    Y_p2.loc[post_induction.index.tolist(), 'os_time'] = Y_p2.loc[post_induction.index.tolist(), 'os_time'] - post_induction['duration']

    X_p2 = X_p2.reset_index(drop=True)
    Y_p2 = Y_p2.reset_index(drop=True)
    df_p2 = df_p2.reset_index(drop=True)

    # print('phase 2', Y_p2['os_time'].min(), Y_p2['os_time'].max())

    # G3 Phase 2 progress vs non-progress
    df_p2_prog = progress_p2(X_p2, Y_p2, df_p2, test_x, model_param_list)

    # G4 Phase 2 progress alive vs progress dec
    df_p2_prog_alive = progress_dec_p2(X_p2, Y_p2, df_p2, test_x, model_param_list)

    # G5 Phase 2 non-prog dec vs rest in p2
    df_p2_non_prog_alive = non_progress_dec_p2(X_p2, Y_p2, df_p2, test_x, model_param_list)

    return df_p2_prog, df_p2_prog_alive, df_p2_non_prog_alive


def train_survival_(x, y, test_x, group, model_param):
    # print('group is', group)
    train_x = x.reset_index(drop=True)
    train_y = y.reset_index(drop=True)

    surv_curv = ''

    cols_x = list(train_x)

    df_train_y = train_y.copy()
    # train_y = Surv.from_dataframe('os_event', 'os_time', train_y)

    if args.surv_model == 'neural_cox_non_prop':
        surv_curv = cox_time(train_x, df_train_y, test_x, cols_x, model_param)

    surv_curv.columns = [group]

    return surv_curv


def cycle_single_sample(sample, treat_combo, X_not_norm, Y, df_X_entire, test_x, model_param_list):
    ########## Phase 1 groups ################
    df_p1_prog, df_p1_prog_alive, df_move_to_phase2 = phase_1(X_not_norm, Y, df_X_entire, test_x, model_param_list)

    ########## Phase 2 groups ##########
    df_p2_prog, df_p2_prog_alive, df_p2_non_prog_dec = phase_2(X_not_norm, Y, df_X_entire, test_x, model_param_list)

    df_p1_prog.fillna(1, inplace=True)
    df_p1_prog_alive.fillna(1, inplace=True)
    df_move_to_phase2.fillna(1, inplace=True)
    df_p2_prog.fillna(1, inplace=True)
    df_p2_prog_alive.fillna(1, inplace=True)
    df_p2_non_prog_dec.fillna(1, inplace=True)

    df_p1_prog.index = df_p1_prog.index.to_series() /365
    df_p1_prog_alive.index = df_p1_prog_alive.index.to_series() / 365
    df_move_to_phase2.index = df_move_to_phase2.index.to_series() / 365

    df_p2_prog.index = df_p2_prog.index.to_series() / 365
    df_p2_prog_alive.index = df_p2_prog_alive.index.to_series() / 365
    df_p2_non_prog_dec.index = df_p2_non_prog_dec.index.to_series() / 365

    join_1 = pd.merge(df_p1_prog, df_p1_prog_alive, left_index=True, right_index=True, how='outer')
    join_2 = pd.merge(join_1, df_p2_prog, left_index=True, right_index=True, how='outer')
    join_3 = pd.merge(join_2, df_p2_prog_alive, left_index=True, right_index=True, how='outer')
    join_4 = pd.merge(join_3, df_p2_non_prog_dec, left_index=True, right_index=True, how='outer')
    join_5 = pd.merge(join_4, df_move_to_phase2, left_index=True, right_index=True, how='outer')

    join_5 = join_5.reset_index().rename(columns={'index' :'Time (years)'})
    join_5['sample'] = sample

    # print('here', sample)
    join_5.head(args.max_time_days).to_csv(args.path + '/loo/model_probs/' + sample + '~' + str(treat_combo) +
                                           '~' + args.surv_model + '.csv', index=False)

    '''
    df_p1_prog.plot(drawstyle='steps-post', ax=ax)
    df_p1_prog_alive.plot(drawstyle='steps-post', ax=ax)
    df_move_to_phase2.plot(drawstyle='steps-post', ax=ax)
    df_p2_prog.plot(drawstyle='steps-post', ax=ax)
    df_p2_prog_alive.plot(drawstyle='steps-post', ax=ax)
    df_p2_non_prog_dec.plot(drawstyle='steps-post', ax=ax)
    '''

    return None


def data_prep_all_samples():
    X_not_norm, Y, df_X_entire = raw_data_read_prepare()
    df_X_entire = df_X_entire.reset_index(drop=True)
    df_X_entire = df_X_entire.sort_values(by=['sample'])
    return df_X_entire['sample'].tolist()


def treat_combo_predictions(model_param_list):
    """
    Sample predictions on different combinations of therapy in P2
    """
    sample = args.sample

    X_not_norm, Y, df_X_entire, test_x = data_prep(sample)
    print('orig treat\n', test_x[['chemo', 'KAR', 'LEN', 'BORT', 'THAL', 'SCT_first_line', 'continuos_treat']])

    treat_on_sample = test_x[['chemo', 'KAR', 'LEN', 'BORT', 'THAL', 'SCT_first_line', 'continuos_treat']]
    treat_on_sample = treat_on_sample[treat_on_sample[['chemo', 'KAR', 'LEN', 'BORT', 'THAL',
                                                'SCT_first_line', 'continuos_treat']].ne(0)].dropna(axis=1, how='any')

    try:
        sct_ = treat_on_sample['SCT_first_line'].values[0]
    except Exception:
        sct_ = 0

    try:
        cont_treat = treat_on_sample['continuos_treat'].values[0]
    except Exception:
        cont_treat = 0

    orig_treats = ",".join(treat_on_sample)
    orig_treats = re.sub(',SCT_first_line|,continuos_treat', '', orig_treats)
    # orig_treats_list = orig_treats.split(',')

    df_treat_combo = pd.read_csv(args.path + '/preds_treatment/combo_freq.csv')
    df_treat_combo = df_treat_combo[df_treat_combo['Freq'] > 10] # change threshold of frequency (with dataset)
    treat_combos = df_treat_combo['Combo']

    for phase in ['P2_no_SCT_no_cont', 'P2_SCT_only', 'P2_cont_only', 'P2_SCT_cont']:
        #if phase == 'P1':
        #    test_x[['SCT_first_line', 'continuos_treat']] = 0
        if phase == 'P2_SCT_only':
            test_x['SCT_first_line'] = 1
            test_x['continuos_treat'] = 0
        elif phase == 'P2_cont_only':
            test_x['continuos_treat'] = 1
            test_x['SCT_first_line'] = 0
        elif phase == 'P2_SCT_cont':
            test_x[['SCT_first_line', 'continuos_treat']] = 1
        elif phase == 'P2_no_SCT_no_cont':
            test_x[['SCT_first_line', 'continuos_treat']] = 0

        # orig treatment
        if (orig_treats not in treat_combos) & (sct_ == test_x['SCT_first_line'].values) & (cont_treat == test_x['continuos_treat'].values):
            cycle_single_sample(sample, phase + ',orig,' + orig_treats, X_not_norm, Y, df_X_entire, test_x, model_param_list)

        if (orig_treats not in treat_combos) | (phase == 'P2_no_SCT_no_cont') | (sct_ == test_x['SCT_first_line'].values) | (cont_treat == test_x['continuos_treat'].values):
            cycle_single_sample(sample, phase + ',' + orig_treats, X_not_norm, Y, df_X_entire, test_x, model_param_list)

        # other combinations
        for ind, treat_combo in enumerate(treat_combos):
            if orig_treats != treat_combo:
                test_x_copy = test_x.copy()
                test_x_copy[['chemo', 'KAR', 'LEN', 'BORT', 'THAL']] = 0
                test_x_copy[treat_combo.split(',')] = 1
                cycle_single_sample(sample, phase + ',' + str(treat_combo), X_not_norm, Y, df_X_entire, test_x_copy, model_param_list)


def main():
    """
    Script to build a model to predict progression or no progression
    :return: None
    """
    global args, f
    args = parser.parse_args()

    model_param = {}; model_param_list = []

    if args.surv_model == 'neural_cox_non_prop':
        model_param['num_nodes'] = [[128, 128], [128, 96], [64, 64], [128, 96], [128, 128], [98, 36]]
        model_param['batch_size'] = [128, 64, 128, 128, 128, 128]
        model_param['seed'] = [666, 444, 666, 787, 888, 888]
        model_param['rate_optim'] = [0.01, 0.001, 0.01, 0.01, 0.01, 0.01]
        model_param['dropout'] = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
        model_param_list = [(i, j, k, l, m) for i, j, k, l, m in zip(model_param['num_nodes'],
                                        model_param['batch_size'], model_param['seed'],
                                        model_param['rate_optim'], model_param['dropout'] )]

    treat_combo_predictions(model_param_list)


if __name__ == '__main__':
    main()